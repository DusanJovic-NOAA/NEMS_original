#include "../../ESMFVersionDefine.h"

#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)
#undef ESMF_520r
#else
#define ESMF_520r
#endif

!-----------------------------------------------------------------------
!
      MODULE MODULE_PARENT_CHILD_CPL_COMP
!
!-----------------------------------------------------------------------
!
!***  This module contains the coupler that exchanges data between
!***  NMM-B parent domains and their children.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!
!   2008-06-12  Black - Module created.
!   2009-02-19  Black - Hydrostatic update of nest boundaries.
!   2010-01-20  Black - Enable parent tasks to update associations
!                       with nest boundary tasks throughout the
!                       integration.
!   2011-02     Yang  - Updated to use both the ESMF 4.0.0rp2 library,
!                       ESMF 5 series library and the the
!                       ESMF 3.1.0rp2 library.
!   2011-05-12  Yang  - Modified for using the ESMF 5.2.0r_beta_snapshot_07.
!   2011-07-16  Black - Add moving nest capability.
!   2011-09-27  Yang  - Modified for using the ESMF 5.2.0r library.
!
!-----------------------------------------------------------------------
!
! USAGE: 
!
!-----------------------------------------------------------------------
!
#ifdef ESMF_520r
      USE esmf
#else
      USE esmf_mod
#endif
!
      USE module_INCLUDE
!
      USE module_CONTROL,ONLY: TIMEF
!
      USE module_DM_PARALLEL,ONLY: MPI_COMM_COMP,MYPE_SHARE
!
      USE module_VARS,ONLY: VAR
!
      USE module_NESTING,ONLY: BNDS_2D                                  &
                              ,BUNDLE_X                                 &
                              ,CHECK_REAL                               &
                              ,INTERIOR_DATA_FROM_PARENT                &
                              ,MIXED_DATA                               &
                              ,REAL_DATA                                &
                              ,MIXED_DATA_TASKS                         &
                              ,REAL_DATA_TASKS                          &
                              ,CHILD_UPDATE_LINK                        &
                              ,MOVING_NEST_BOOKKEEPING                  &
                              ,MOVING_NEST_RECV_DATA                    &
                              ,PARENT_BOOKKEEPING_MOVING                &
                              ,PARENT_READS_MOVING_CHILD_TOPO           &
                              ,PARENT_UPDATES_HALOS                     &
                              ,PARENT_UPDATES_MOVING                    &
                              ,REAL_DATA_2D
!
!!!   USE module_CLOCKTIMES
!
      USE module_CONSTANTS,ONLY: G,P608,R_D
!
      USE module_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
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
               ,PARENT_CHILD_COUPLER_SETUP                              &
               ,NSTEP_CHILD_RECV
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
        REAL,DIMENSION(:),ALLOCATABLE :: STRING
        INTEGER :: LENGTH
        INTEGER :: ID_SOURCE
        INTEGER :: INDX_START
        INTEGER :: INDX_END
        INTEGER :: INDX_END_EXP
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
      INTEGER(kind=KINT),PARAMETER :: INDX_SW=1                         &
                                     ,INDX_SE=2                         &
                                     ,INDX_NW=3                         &
                                     ,INDX_NE=4
!
      INTEGER(kind=KINT),SAVE :: MOVE_TAG=1111                          &  !<-- Arbitrary tag used for child's move
                                ,MOVING_BC_TAG=1112                     &  !<-- Arbitrary tag used for moving nests' BC updates
                                ,PARENT_SHIFT_TAG=1113                  &  !<-- Arbitrary tag used for parent's move
                                ,TIME_TAG=1114                             !<-- Arbitrary tag used for child's move timestep
!
      INTEGER(kind=KINT),SAVE :: NEXT_MOVE_TIMESTEP=-999999                !<-- For all moving nests (a child or a parent)
!
      INTEGER(kind=KINT),SAVE :: COMM_MY_DOMAIN                         &
                                ,COMM_TO_MY_PARENT                      &
                                ,I_SHIFT_CHILD                          &  !<-- The shift of a moving nest in its own space.
                                ,J_SHIFT_CHILD                          &  !
                                ,I_SW_PARENT_CURRENT                    &  !<-- These are scalars used by the nests for their own 
                                ,J_SW_PARENT_CURRENT                    &  !     values of their parent I,J on their SW corner.
                                ,ITS,ITE,JTS,JTE,LM                     &
                                ,IMS,IME,JMS,JME                        &
                                ,IDS,IDE,JDS,JDE                        &
                                ,IM_1,INDX_CW,INDX_Q                    &
                                ,INPES,INPES_PARENT                     &
                                ,JM_1,JNPES,JNPES_PARENT                &
                                ,KOUNT_RATIOS_MN                        &
                                ,MY_DOMAIN_ID                           &
                                ,MY_LOCAL_RANK_CHILD                    &
                                ,MY_LOCAL_RANK_PARENT                   &
                                ,MYPE                                   &
                                ,N_BLEND_H,N_BLEND_V                    &
                                ,NHALO                                  &
                                ,NPHS                                   &
                                ,NROWS_P_UPD_E                          &
                                ,NROWS_P_UPD_N                          &
                                ,NROWS_P_UPD_S                          &
                                ,NROWS_P_UPD_W                          &
                                ,NUM_CHILDREN                           &
                                ,NUM_FIELDS_MOVE                        &
                                ,NUM_FIELDS_MOVE_2D_H_I                 &
                                ,NUM_FIELDS_MOVE_2D_X_I                 &
                                ,NUM_FIELDS_MOVE_2D_H_R                 &
                                ,NUM_FIELDS_MOVE_2D_X_R                 &
                                ,NUM_FIELDS_MOVE_3D_H                   &
                                ,NUM_FIELDS_MOVE_2D_V                   &
                                ,NUM_FIELDS_MOVE_3D_V                   &
                                ,NUM_LEVELS_MOVE_3D_H                   &
                                ,NUM_LEVELS_MOVE_3D_V                   &
                                ,NUM_MOVING_CHILDREN                    &
                                ,NUM_PES_FCST                           &
                                ,NUM_TASKS_PARENT                       &
                                ,SPACE_RATIO_MY_PARENT                  &
                                ,TIME_RATIO_MY_PARENT
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,SAVE :: LOCAL_ISTART      &
                                                     ,LOCAL_IEND        &
                                                     ,LOCAL_JSTART      &
                                                     ,LOCAL_JEND
!
      INTEGER(kind=KINT),DIMENSION(1:2) :: IJ_SHIFT_CHILD               &
                                          ,PARENT_SHIFT
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,SAVE :: HANDLE_BC_UPDATE    &
                                                         ,HANDLE_TIMESTEP     &
                                                         ,IM_CHILD            &
                                                         ,JM_CHILD            &
                                                         ,I_PARENT_SW         &  !<-- These are arrays used by parents to hold their
                                                         ,J_PARENT_SW         &  !     values of I,J for all their nests' SW corners.
                                                         ,ITE_PARENT          &
                                                         ,ITS_PARENT          &
                                                         ,JTE_PARENT          &
                                                         ,JTS_PARENT          &
                                                         ,LINK_MRANK_RATIO    &
                                                         ,LIST_OF_RATIOS      &
                                                         ,M_NEST_RATIO        &
                                                         ,N_BLEND_H_CHILD     &
                                                         ,N_BLEND_V_CHILD     &
                                                         ,NSTEP_CHILD_RECV    &
!
                                                         ,NUM_TASKS_SEND_H_S  &
                                                         ,NUM_TASKS_SEND_H_N  &
                                                         ,NUM_TASKS_SEND_H_W  &
                                                         ,NUM_TASKS_SEND_H_E  &
                                                         ,NUM_TASKS_SEND_V_S  &
                                                         ,NUM_TASKS_SEND_V_N  &
                                                         ,NUM_TASKS_SEND_V_W  &
                                                         ,NUM_TASKS_SEND_V_E 
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,SAVE :: COMM_TO_MY_CHILDREN       &
                                                     ,FTASKS_DOMAIN             &
                                                     ,ID_PARENTS                &
                                                     ,INC_FIX                   &
                                                     ,MY_CHILDREN_ID            &
                                                     ,PARENT_CHILD_SPACE_RATIO  &
                                                     ,RANK_MOVING_CHILD         &
                                                     ,TIME_RATIO_MY_CHILDREN
!
      REAL(kind=KFPT) :: EPS=1.E-4
!
      REAL(kind=KFPT),SAVE :: DYH,PDTOP,PT                                   &
                             ,RECIP_DPH_1,RECIP_DLM_1                        &
                             ,SB_1,WB_1                                      &
                             ,TPH0_1,TLM0_1
!
      REAL(kind=KFPT),DIMENSION(:),POINTER,SAVE :: BOUND_1D_SOUTH_H          &
                                                  ,BOUND_1D_SOUTH_V          &
                                                  ,BOUND_1D_NORTH_H          &
                                                  ,BOUND_1D_NORTH_V          &
                                                  ,BOUND_1D_WEST_H           &
                                                  ,BOUND_1D_WEST_V           &
                                                  ,BOUND_1D_EAST_H           &
                                                  ,BOUND_1D_EAST_V           &
!
                                                  ,CHILD_PARENT_SPACE_RATIO  &
!
                                                  ,DT_DOMAIN                 &
!
                                                  ,DSG2                      &
                                                  ,DXH                       &
                                                  ,PDSG1                     &
                                                  ,PSGML1                    &
                                                  ,SG1                       &
                                                  ,SG2                       &
                                                  ,SGML2
!
      REAL(kind=KFPT),DIMENSION(:,:),POINTER,SAVE :: FIS                     &
                                                    ,FIS_CHILD_ON_PARENT     &
                                                    ,GLAT                    &
                                                    ,GLON                    &
                                                    ,PD                      &
                                                    ,SM
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER,SAVE :: CW,PINT,Q,T,U,V
!
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE,SAVE :: PDB_S,PDB_N    &
                                                        ,PDB_W,PDB_E
!
#ifdef ESMF_3
      TYPE(ESMF_Logical),SAVE :: I_AM_A_FCST_TASK                       &
                                ,I_AM_A_PARENT                          &
                                ,MY_DOMAIN_MOVES
#else
      LOGICAL(kind=KLOG),SAVE :: I_AM_A_FCST_TASK                       &
                                ,I_AM_A_PARENT                          &
                                ,MY_DOMAIN_MOVES
#endif
!
      REAL(kind=KFPT),DIMENSION(:,:,:),ALLOCATABLE,SAVE :: CWB_S,CWB_N  &
                                                          ,CWB_W,CWB_E  &
                                                          ,QB_S,QB_N    &
                                                          ,QB_W,QB_E    &
                                                          ,TB_S,TB_N    &
                                                          ,TB_W,TB_E    &
                                                          ,UB_S,UB_N    &
                                                          ,UB_W,UB_E    &
                                                          ,VB_S,VB_N    &
                                                          ,VB_W,VB_E
!
      REAL(kind=KFPT),DIMENSION(:,:,:,:),POINTER,SAVE :: TRACERS
!
      LOGICAL(kind=KLOG),SAVE :: I_WANT_TO_MOVE                         &
                                ,MOVE_FLAG_SENT                         &
                                ,MY_PARENT_MOVES                        &
                                ,RECVD_MOVE_TIMESTEP
!
      LOGICAL(kind=KLOG),DIMENSION(:),ALLOCATABLE,SAVE :: MOVE_FLAG     &
                                                         ,SEND_CHILD_DATA
!
      CHARACTER(23) :: FILENAME_MOVING_NEST_RES='nest_res_topo_on_parent'
!
      CHARACTER(6),DIMENSION(:),ALLOCATABLE,SAVE :: STATIC_OR_MOVING
!
      TYPE(BNDS_2D),DIMENSION(:),ALLOCATABLE,SAVE ::                    &
                                               NEST_FIS_ON_PARENT_BNDS
!
      TYPE(REAL_DATA_2D),DIMENSION(:),ALLOCATABLE,SAVE ::               &
                                                   NEST_FIS_ON_PARENT   &
                                                  ,NEST_FIS_V_ON_PARENT
!
      TYPE(MIXED_DATA_TASKS),DIMENSION(:),ALLOCATABLE,SAVE ::           &
                                                   MOVING_CHILD_UPDATE
!
      TYPE(REAL_DATA_TASKS),DIMENSION(:,:),POINTER,SAVE ::              &
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
      TYPE(HANDLE_SEND),DIMENSION(:),POINTER,SAVE :: HANDLE_MOVE_DATA 
!
      TYPE(HANDLE_SEND),DIMENSION(:,:),POINTER,SAVE :: HANDLE_H_SOUTH   &
                                                      ,HANDLE_H_NORTH   &
                                                      ,HANDLE_H_WEST    &
                                                      ,HANDLE_H_EAST    &
                                                      ,HANDLE_V_SOUTH   &
                                                      ,HANDLE_V_NORTH   &
                                                      ,HANDLE_V_WEST    &
                                                      ,HANDLE_V_EAST
!
!!!   TYPE(CHILD_UPDATE_LINK),DIMENSION(:),ALLOCATABLE,target,SAVE ::          &
!!!                                                TASK_UPDATE_SPECS_H  &
!!!                                               ,TASK_UPDATE_SPECS_V
      TYPE(CHILD_UPDATE_LINK),DIMENSION(:),POINTER,SAVE ::              &
                                                   TASK_UPDATE_SPECS
!
      TYPE(ESMF_FieldBundle),SAVE :: MOVE_BUNDLE_H                      &
                                    ,MOVE_BUNDLE_V
!
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: btim,btim0
!
      REAL(kind=KDBL),SAVE :: cpl1_prelim_tim                           &
                             ,cpl1_south_h_tim,cpl1_south_v_tim         &
                             ,cpl1_north_h_tim,cpl1_north_v_tim         &
                             ,cpl1_west_h_tim, cpl1_west_v_tim          &
                             ,cpl1_east_h_tim, cpl1_east_v_tim          &
                             ,cpl1_recv_tim
!
      REAL(kind=KDBL),SAVE :: cpl1_south_h_recv_tim                     &
                             ,cpl1_south_h_undo_tim                     &
                             ,cpl1_south_h_exp_tim                      &
                             ,cpl1_south_v_recv_tim                     &
                             ,cpl1_south_v_undo_tim                     &
                             ,cpl1_south_v_exp_tim
!
      REAL(kind=KDBL),SAVE :: cpl2_comp_tim                             &
                             ,cpl2_send_tim                             &
                             ,cpl2_wait_tim
!
      REAL(kind=KDBL),SAVE :: moving_nest_bookkeep_tim                  &
                             ,moving_nest_update_tim
!
      REAL(kind=KDBL),SAVE :: parent_bookkeep_moving_tim                &
                             ,parent_update_moving_tim                  &
                             ,t0_recv_move_tim
!
      character(len=15) :: hhmmssffff
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
!***  Register the nesting coupler component's Initialize, Run, and 
!***  Finalize routines.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_CplComp) :: CPL_COMP                                      !<-- Coupler component
!
      INTEGER(kind=KINT),INTENT(OUT) :: RC_NEST_REG                                   !<-- Return code for register
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
#ifdef ESMF_3
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETINIT                       &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_INITIALIZE        &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
#else
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETINIT                       &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_INITIALIZE        &  !<-- User's subroutineName
                                    ,phase=ESMF_SINGLEPHASE             &  !<-- Phase
                                    ,rc=RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NEST_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Register the Parent-Child coupler's Run subroutines.
!
!     The Parent-Child Run step of the coupler has three distinct parts:
!     (1) The relevant parent tasks compute interpolation information
!         needed for generating data to be sent to the nests.
!     (2) At the top of certain timesteps a child receives data 
!         from its parent.
!     (3) At the end of every timestep a parent sends data to its
!         children.  For those children that moved, the parent must
!         regenerate new interpolation information before computing
!         the new data for those nests.
!
!***  Thus register the coupler's Run step with those two phases.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NEST_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Entry Point for Phase 1 of Nesting Coupler Run"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_RUN_RECV          &  !<-- User's subroutineName
                                    ,1                                  &  !<-- Phase
                                    ,RC)
#else
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_RUN_RECV          &  !<-- User's subroutineName
                                    ,phase=1                            &  !<-- Phase
                                    ,rc=RC)
#endif
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
#ifdef ESMF_3
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_RUN_SEND          &  !<-- User's subroutineName
                                    ,2                                  &  !<-- Phase
                                    ,RC)
#else
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETRUN                        &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_RUN_SEND          &  !<-- User's subroutineName
                                    ,phase=2                            &  !<-- Phase
                                    ,rc=RC)
#endif
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
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETFINAL                      &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_FINALIZE          &  !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE                   &  !<-- Phase
                                    ,RC)
#else
      CALL ESMF_CplCompSetEntryPoint(CPL_COMP                           &  !<-- The Parent-Child Coupler Component
                                    ,ESMF_SETFINAL                      &  !<-- subroutineType
                                    ,PARENT_CHILD_CPL_FINALIZE          &  !<-- User's subroutineName
                                    ,phase=ESMF_SINGLEPHASE             &  !<-- Phase
                                    ,rc=RC)
#endif
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
!***  Perform initial work needed by the Parent-Child coupler.
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
      INTEGER,INTENT(OUT) :: RC_FINAL
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J,L
!
      INTEGER(kind=KINT) :: CHILD_ID,CONFIG_ID                          &
                           ,I_RATIO,ICORNER                             &
                           ,ID_DOM,ID_MY_PARENT                         &
                           ,IDIM,IEND,ISTART,IUNIT_FIS_NEST             &
                           ,JCORNER,JDIM,JEND,JSTART,JSTOP              &
                           ,KOUNT,KR                                    &
                           ,LENGTH,LIM1_H,LIM1_V,LIM2_H,LIM2_V,LMP1,LOR &
                           ,MAX_DOMAINS                                 &
                           ,N,NN,N_CHILD,N_FIELD,N_START,N_END          &
                           ,N_H_EAST_WEST,N_H_NORTH_SOUTH               &
                           ,N_V_EAST_WEST,N_V_NORTH_SOUTH               &
                           ,N_MOVING,NROWS_P_UPD_X                      &
                           ,NKOUNT,NTIMESTEP                            &
                           ,NUM_BOUNDARY_WORDS,NUM_CHILD_TASKS          &
                           ,NUM_DIMS,NUM_DOMAINS                        &
                           ,SFC_FILE_RATIO                              &
                           ,UPDATE_TYPE_INT
!
      INTEGER(kind=KINT) :: IERR,RC,RC_CPL_INIT
!
      INTEGER(kind=KINT),DIMENSION(1:3) :: LBNDS_DYN,UBNDS_DYN          &
                                          ,LBNDS_PHY,UBNDS_PHY          &
                                          ,LIMITS_LO,LIMITS_HI
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER(kind=KINT),DIMENSION(1:3) :: INFO_EXT_DATA
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER :: DOMAIN_ID_TO_RANK
!
      REAL(kind=KFPT) :: DIST_NESTV_SOUTH_TO_PARENTV_SOUTH              &
                        ,DT_PARENT                                      &
                        ,REAL_I                                         &
                        ,REAL_J
!
      REAL(kind=KFPT) :: DPH_1,DLM_1                                    &
                        ,SBD_1,WBD_1                                    &
                        ,TPH0D_1,TLM0D_1
!
      REAL(kind=DOUBLE) :: D2R,D_ONE,D_180,PI
!
!!!   REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: HOLD_FIS_1D
!
!!!   REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE :: HOLD_FIS_NEST       &
!!!                                                ,HOLD_FIS_V_NEST
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: DUMMY_3D=>NULL()
!
      CHARACTER(len=1) :: UPDATE_TYPE_CHAR
!
      CHARACTER(len=2) :: INT_TO_CHAR
      CHARACTER(len=6) :: FMT='(I2.2)'
!
      CHARACTER(len=99) :: CONFIG_FILE_NAME
!
      LOGICAL(kind=KLOG) :: DOMAIN_MOVES,FOUND,OPENED
!
#ifdef ESMF_3
      TYPE(ESMF_Logical) :: I_AM_A_NEST
#else
      LOGICAL(kind=KLOG) :: I_AM_A_NEST
#endif
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(ESMF_Config) :: CF_1,CF_MINE,CF_PARENT
      TYPE(ESMF_Config),DIMENSION(:),ALLOCATABLE :: CF
!
      TYPE(ESMF_TypeKind) :: DATATYPE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      MYPE=MYPE_SHARE
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
!-------------------------------
!***  Maximum number of domains
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Maximum Number of Domains"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='MAX_DOMAINS'                        &  !<-- Name of the attribute to extract
                            ,value=MAX_DOMAINS                          &  !<-- Maximum # of domains
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
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
                            ,value=NUM_CHILDREN                         &  !<-- How many children does this DOMAIN Component have?
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NUM_MOVING_CHILDREN=0                                                !<-- Initialize the number of children who will move
!
!----------------------------------
!***  Communicator for each domain
!----------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Domain Intercommunicator in P-C Coupler Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='Single Domain Comm'                 &  !<-- Name of the attribute to extract
                            ,value=COMM_MY_DOMAIN                       &  !<-- MPI communicator for this domain
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------
!***  Number of fcst tasks on this domain
!-----------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract # of Forecast Tasks in P-C Coupler Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='NUM_PES_FCST'                       &  !<-- Name of the attribute to extract
                            ,value=NUM_PES_FCST                         &  !<-- MPI communicator for this domain
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
      MESSAGE_CHECK="Extract Child-to-Parent Comm in P-C Coupler Init"
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
      MESSAGE_CHECK="Extract Total Number of Domains in P-C Coupler Init"
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
!----------------------------------------------------------
!***  The association of domains and their configure files
!----------------------------------------------------------
!
      ALLOCATE(DOMAIN_ID_TO_RANK(1:MAX_DOMAINS))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Association of Domains and Config Files"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='DOMAIN_ID_TO_RANK'              &  !<-- Name of the attribute to extract
                            ,count    =MAX_DOMAINS                      &  !<-- Number of elements in Attribute
                            ,valueList=DOMAIN_ID_TO_RANK                &  !<-- Array associating domains and config files
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='DOMAIN_ID_TO_RANK'              &  !<-- Name of the attribute to extract
                            ,itemCount=MAX_DOMAINS                      &  !<-- Name of elements in the Attribute
                            ,valueList=DOMAIN_ID_TO_RANK                &  !<-- Array associating domains and config files
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------
!***  Fundamental Timestep on Each Domain
!-----------------------------------------
!
      ALLOCATE(DT_DOMAIN(1:NUM_DOMAINS))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Timestep of Domains in P-C Coupler Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='DOMAIN_DTs'                     &  !<-- Name of the attribute to extract
                            ,count    =NUM_DOMAINS                      &  !<-- # of items in the Attribute
                            ,valueList=DT_DOMAIN                        &  !<-- Timestep on each domain
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='DOMAIN_DTs'                     &  !<-- Name of the attribute to extract
                            ,itemCount=NUM_DOMAINS                      &  !<-- # of items in the Attribute
                            ,valueList=DT_DOMAIN                        &  !<-- Timestep on each domain
                            ,rc       =RC)
#endif
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
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='FTASKS_DOMAIN'                  &  !<-- Name of the attribute to extract
                            ,count    =NUM_DOMAINS                      &  !<-- # of items in the Attribute
                            ,valueList=FTASKS_DOMAIN                    &  !<-- # of forecast tasks on each domain
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='FTASKS_DOMAIN'                  &  !<-- Name of the attribute to extract
                            ,itemCount=NUM_DOMAINS                      &  !<-- # of items in the Attribute
                            ,valueList=FTASKS_DOMAIN                    &  !<-- # of forecast tasks on each domain
                            ,rc       =RC)
#endif
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
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='ID_PARENTS'                     &  !<-- Name of the attribute to extract
                            ,count    =NUM_DOMAINS                      &  !<-- # of items in the Attribute
                            ,valueList=ID_PARENTS                       &  !<-- Domain IDs of parents
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='ID_PARENTS'                     &  !<-- Name of the attribute to extract
                            ,itemCount=NUM_DOMAINS                      &  !<-- # of items in the Attribute
                            ,valueList=ID_PARENTS                       &  !<-- Domain IDs of parents
                            ,rc       =RC)
#endif
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
!-----------------------------------------------------
!***  Index Limits of All Forecast Tasks on My Domain
!-----------------------------------------------------
!
      ALLOCATE(LOCAL_ISTART(1:FTASKS_DOMAIN(MY_DOMAIN_ID)))
      ALLOCATE(LOCAL_IEND  (1:FTASKS_DOMAIN(MY_DOMAIN_ID)))
      ALLOCATE(LOCAL_JSTART(1:FTASKS_DOMAIN(MY_DOMAIN_ID)))
      ALLOCATE(LOCAL_JEND  (1:FTASKS_DOMAIN(MY_DOMAIN_ID)))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Index Limits of Fcst Tasks on My Domain in Init Step Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='LOCAL ISTART'                   &  !<-- Name of the attribute to extract
                            ,count    =FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- # of items in the Attribute
                            ,valueList=LOCAL_ISTART                     &  !<-- Starting I's of fcst tasks on my domain
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='LOCAL IEND'                     &  !<-- Name of the attribute to extract
                            ,count    =FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- # of items in the Attribute
                            ,valueList=LOCAL_IEND                       &  !<-- Ending I's of fcst tasks on my domain
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='LOCAL JSTART'                   &  !<-- Name of the attribute to extract
                            ,count    =FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- # of items in the Attribute
                            ,valueList=LOCAL_JSTART                     &  !<-- Starting J's of fcst tasks on my domain
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='LOCAL JEND'                     &  !<-- Name of the attribute to extract
                            ,count    =FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- # of items in the Attribute
                            ,valueList=LOCAL_JEND                       &  !<-- Ending J's of fcst tasks on my domain
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='LOCAL ISTART'                   &  !<-- Name of the attribute to extract
                            ,itemCount=FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- # of items in the Attribute
                            ,valueList=LOCAL_ISTART                     &  !<-- Starting I's of fcst tasks on my domain
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='LOCAL IEND'                     &  !<-- Name of the attribute to extract
                            ,itemCount=FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- # of items in the Attribute
                            ,valueList=LOCAL_IEND                       &  !<-- Ending I's of fcst tasks on my domain
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='LOCAL JSTART'                   &  !<-- Name of the attribute to extract
                            ,itemCount=FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- # of items in the Attribute
                            ,valueList=LOCAL_JSTART                     &  !<-- Starting J's of fcst tasks on my domain
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='LOCAL JEND'                     &  !<-- Name of the attribute to extract
                            ,itemCount=FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- # of items in the Attribute
                            ,valueList=LOCAL_JEND                       &  !<-- Ending J's of fcst tasks on my domain
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
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
!***  Unload the domains' prognostic data.
!-----------------------------------------------------------------------
!
!--------
!***  PD
!--------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PD Field from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='PD'                                  &  !<-- Extract PD
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PD from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =PD                                     &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=PD                                   &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------
!***  Layer Interface Pressures
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PINT from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='PINT'                                &  !<-- Extract layer interface pressures
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PINT from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =PINT                                   &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=PINT                                 &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
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
      MESSAGE_CHECK="Extract T Field from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='T'                                   &  !<-- Extract temperature
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract T from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =T                                      &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=T                                    &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
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
      MESSAGE_CHECK="Extract U Field from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='U'                                   &  !<-- Extract U wind
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract U from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =U                                      &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=U                                    &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
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
      MESSAGE_CHECK="Extract V Field from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='V'                                   &  !<-- Extract V wind
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract V from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =V                                      &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=V                                    &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
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
      MESSAGE_CHECK="Extract Tracers Field from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='TRACERS'                             &  !<-- Extract tracers
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Tracers from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =TRACERS                                &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=TRACERS                              &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
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
      MESSAGE_CHECK="Extract FIS Field from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='FIS'                                 &  !<-- Extract sfc geopotential
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract FIS from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =FIS                                    &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=FIS                                  &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------
!***  Sea Mask
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Sea Mask from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='SM'                                  &  !<-- Extract sea mask
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract SM from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =SM                                     &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=SM                                   &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------
!***  Geographic latitude
!-------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract GLAT Field from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='GLAT'                                &  !<-- Extract geographic latitude
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract GLAT from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =GLAT                                   &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=GLAT                                 &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------
!***  Geographic longitude
!--------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract GLON Field from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='GLON'                                &  !<-- Extract geographic longitude
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract GLON from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =GLON                                   &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=GLON                                 &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------------
!***  PT,PDTOP,PSGML1,SG1,SG2,SGML2,DSG2,PDSG1
!----------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PT from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='PT'                                 &  !<-- Extract PT
                            ,value=PT                                   &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PDTOP from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='PDTOP'                              &  !<-- Extract PDTOP
                            ,value=PDTOP                                &  !<-- Put the extracted Attribute here
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
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='PSGML1'                         &  !<-- Extract PSGML1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=PSGML1                           &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='PSGML1'                         &  !<-- Name of Attribute to extract
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=PSGML1                           &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      LMP1=LM+1
      ALLOCATE(SG1(1:LMP1))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract SG1 from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='SG1'                            &  !<-- Extract SG1
                            ,count    =LMP1                             &  !<-- # of words in data list
                            ,valueList=SG1                              &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='SG1'                            &  !<-- Name of Attribute to extract
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=SG1                              &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALLOCATE(SG2(1:LMP1))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract SG2 from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='SG2'                            &  !<-- Extract SG2
                            ,count    =LMP1                             &  !<-- # of words in data list
                            ,valueList=SG2                              &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='SG2'                            &  !<-- Name of Attribute to extract
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=SG2                              &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALLOCATE(SGML2(1:LM))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract SGML2 from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='SGML2'                          &  !<-- Extract SGML2
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=SGML2                            &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='SGML2'                          &  !<-- Name of Attribute to extract
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=SGML2                            &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALLOCATE(DSG2(1:LM))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract DSG2 from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='DSG2'                           &  !<-- Extract DSG2
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=DSG2                             &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='DSG2'                           &  !<-- Extract DSG2
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=DSG2                             &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALLOCATE(PDSG1(1:LM))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PDSG1 from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='PDSG1'                          &  !<-- Extract PDSG1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=PDSG1                            &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='PDSG1'                          &  !<-- Extract PDSG1
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=PDSG1                            &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------
!***  DYH,DXH
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract DYH from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='DYH'                                &  !<-- Extract DYH
                            ,value=DYH                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ALLOCATE(DXH(JDS:JDE))
      NKOUNT=JDE-JDS+1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract DXH from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='DXH'                            &  !<-- Extract DXH
                            ,count    =NKOUNT                           &  !<-- # of words in data list
                            ,valueList=DXH                              &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#else
      CALL ESMF_AttributeGet(state    =IMP_STATE                        &  !<-- The parent-child coupler import state
                            ,name     ='DXH'                            &  !<-- Name of Attribute to extract
                            ,itemCount=NKOUNT                           &  !<-- # of words in data list
                            ,valueList=DXH                              &  !<-- Put the extracted Attribute here
                            ,rc       =RC)
#endif
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
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                            &  !<-- The parent-child coupler import state
                            ,name ='INDX_Q'                             &  !<-- Name of Attribute to extract
                            ,value=INDX_Q                               &  !<-- Put the extracted Attribute here
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
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=IMP_STATE                              &  !<-- The parent-child coupler import state
                            ,name ='INDX_CW'                              &  !<-- Name of Attribute to extract
                            ,value=INDX_CW                                &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!   Q=>TRACERS(:,:,:,INDX_Q)
!!!   CW=>TRACERS(:,:,:,INDX_CW)
      Q=>TRACERS(IMS:IME,JMS:JME,1:LM,INDX_Q)
      CW=>TRACERS(IMS:IME,JMS:JME,1:LM,INDX_CW)
!
!-----------------------------------------------------------------------
!***  Children need to save the ratio of their parent's timestep and
!***  grid increment to their own.  The timestep ratio MUST be an
!***  integer and for now so must the space ratio.
!***  Then obtain the parent I,J of the nest's SW corner.
!-----------------------------------------------------------------------
!
      child_block: IF(COMM_TO_MY_PARENT/=-999)THEN                         !<-- Select the children
!
!-----------------------------------------------------------------------
!
        DT_PARENT=DT_DOMAIN(ID_PARENTS(MY_DOMAIN_ID))
        TIME_RATIO_MY_PARENT=NINT(DT_PARENT/DT_DOMAIN(MY_DOMAIN_ID))       !<-- Ratio of my parent's timestep to mine
!
!-----------------------------------------------------------------------
!***  In order to allow moving nests to be updated their tasks need
!***  to know their domain's forecast task layout as well as that of
!***  their parents.  Likewise the parents of moving nests need to
!***  know the forecast task layout of their moving children.  For
!***  simplicity we will provide that information to all domain tasks
!***  now.
!-----------------------------------------------------------------------
!
        CONFIG_ID=DOMAIN_ID_TO_RANK(MY_DOMAIN_ID)
        WRITE(INT_TO_CHAR,FMT)CONFIG_ID
        CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                    !<-- Prepare the config file name
!
        CF_MINE=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent-Child Init: Nest Loads Its Configure File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigLoadFile(config  =CF_MINE                       &
                                ,filename=CONFIG_FILE_NAME              &
                                ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent-Child Init: Child Gets Space Ratio"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF_MINE                     &  !<-- The child's config object
                                    ,value =SPACE_RATIO_MY_PARENT       &  !<-- The variable filled (Parent-to-child space ratio)
                                    ,label ='parent_child_space_ratio:' &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent-Child Init: Child Gets SW Corner Point"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF_MINE                     &  !<-- The child's config object
                                    ,value =I_SW_PARENT_CURRENT         &  !<-- The variable filled (parent I of nest SW corner)
                                    ,label ='i_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF_MINE                     &  !<-- The child's config object
                                    ,value =J_SW_PARENT_CURRENT         &  !<-- The variable filled (parent J of nest SW corner)
                                    ,label ='j_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent-Child Init: Child Gets INPES,JNPES"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF_MINE                     &  !<-- The child's config object
                                    ,value =INPES                       &  !<-- The variable filled (fcst tasks in I direction)
                                    ,label ='inpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF_MINE                     &  !<-- The child's config object
                                    ,value =JNPES                       &  !<-- The variable filled (fcst tasks in J direction)
                                    ,label ='jnpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Move Flag from Nest's Configure File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
        CALL ESMF_ConfigGetAttribute(config=CF_MINE                     &  !<-- The child's config object
                                    ,value =DOMAIN_MOVES                &  !<-- The variable filled (Move flag)
                                    ,label ='my_domain_moves:'          &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        MY_DOMAIN_MOVES=ESMF_FALSE
!
        IF(DOMAIN_MOVES)THEN
          MY_DOMAIN_MOVES=ESMF_TRUE
        ENDIF
#else
        CALL ESMF_ConfigGetAttribute(config=CF_MINE                     &  !<-- The child's config object
                                    ,value =MY_DOMAIN_MOVES             &  !<-- The variable filled (Move flag)
                                    ,label ='my_domain_moves:'          &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#endif
!
!-----------------------------------------------------------------------
!***  Moving nests need their parents' INPES.
!-----------------------------------------------------------------------
!
#ifdef ESMF_3
        IF(MY_DOMAIN_MOVES==ESMF_TRUE)THEN
#else
        IF(MY_DOMAIN_MOVES)THEN
#endif
!
          ID_MY_PARENT=ID_PARENTS(MY_DOMAIN_ID)
          CONFIG_ID=DOMAIN_ID_TO_RANK(ID_MY_PARENT)
          WRITE(INT_TO_CHAR,FMT)CONFIG_ID
          CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                  !<-- Prepare the config file name
!
          CF_PARENT=ESMF_ConfigCreate(rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Nest Loads Parent Config File"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF_PARENT                   &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Child Gets Parent INPES"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_PARENT                 &  !<-- The parent's config object
                                      ,value =INPES_PARENT              &  !<-- The variable filled (fcst tasks in I direction)
                                      ,label ='inpes:'                  &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Moving nests also must know if their parents move.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Nest Checks If Parent Moves"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_PARENT                 &  !<-- The parent's config object
                                      ,value =MY_PARENT_MOVES           &  !<-- The variable filled (does the parent move?)
                                      ,label ='my_domain_moves:'        &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Since the nests can only move on parent timesteps and
!***  are allowed to move only on physics timesteps then 
!***  warn the user if the Parent timestep ratio does not
!***  divide evenly into the nest's physics frequency.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Child Gets NPHS"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_MINE                   &  !<-- The nest's config object
                                      ,value =NPHS                      &  !<-- The variable filled (frequency of physics calls)
                                      ,label ='nphs:'                   &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(MOD(NPHS,TIME_RATIO_MY_PARENT)/=0)THEN
            WRITE(0,*)' WARNING: Moving nest parent time ratio does'    &
                     ,' not divide into its NPHS!!!'
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF child_block
!
!-----------------------------------------------------------------------
!
      parent_block: IF(NUM_CHILDREN>0)THEN                                 !<-- Select parents for additional setup 
!
!-----------------------------------------------------------------------
!***  First the parent tasks extract relevant data about themselves
!***  from the Parent-Child Coupler Import State.
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
#ifdef ESMF_3
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='Parent-to-Child Comms'        &  !<-- Name of the attribute to extract
                              ,count    =NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=COMM_TO_MY_CHILDREN            &  !<-- MPI communicators to my children
                              ,rc       =RC)
#else
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='Parent-to-Child Comms'        &  !<-- Name of the attribute to extract
                              ,itemCount=NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=COMM_TO_MY_CHILDREN            &  !<-- MPI communicators to my children
                              ,rc       =RC)
#endif
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
#ifdef ESMF_3
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='CHILD_IDs'                    &  !<-- Name of the attribute to extract
                              ,count    =NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=MY_CHILDREN_ID                 &  !<-- The domain IDs of the current domain's children
                              ,rc       =RC)
#else
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='CHILD_IDs'                    &  !<-- Name of the attribute to extract
                              ,itemCount=NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=MY_CHILDREN_ID                 &  !<-- The domain IDs of the current domain's children
                              ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract integer ratios of parent-to-child timesteps.
!-----------------------------------------------------------------------
!
        ALLOCATE(TIME_RATIO_MY_CHILDREN(1:NUM_CHILDREN))                   !<-- Integer ratio of parent timestep to children's
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Parent-to-Child DT Ratio from Imp State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='Parent-Child Time Ratio'      &  !<-- Name of the attribute to extract
                              ,count    =NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=TIME_RATIO_MY_CHILDREN         &  !<-- Ratio of parent to child DTs 
                              ,rc       =RC)
#else
        CALL ESMF_AttributeGet(state    =IMP_STATE                      &  !<-- The parent-child coupler import state
                              ,name     ='Parent-Child Time Ratio'      &  !<-- Name of the attribute to extract
                              ,itemCount=NUM_CHILDREN                   &  !<-- # of items in the Attribute
                              ,valueList=TIME_RATIO_MY_CHILDREN         &  !<-- Ratio of parent to child DTs 
                              ,rc       =RC)
#endif
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
          NSTEP_CHILD_RECV(N)=(NTIMESTEP-1)*TIME_RATIO_MY_CHILDREN(N)
        ENDDO
!
!-----------------------------------------------------------------------
!***  Allocate more arrays needed by the parent to hold child
!***  information derived from the children's configure files.
!-----------------------------------------------------------------------
!
        ALLOCATE(CF(1:NUM_CHILDREN))                                       !<-- Configure objects of this parent's children
        ALLOCATE(IM_CHILD(1:NUM_CHILDREN))                                 !<-- I extent of children's domains
        ALLOCATE(JM_CHILD(1:NUM_CHILDREN))                                 !<-- J extent of children's domains
        ALLOCATE(PARENT_CHILD_SPACE_RATIO(1:NUM_CHILDREN))                 !<-- Integer ratio of parent grid increment to children's
        ALLOCATE(CHILD_PARENT_SPACE_RATIO(1:NUM_CHILDREN))                 !<-- Inverse of PARENT_CHILD_SPACE_RATIO
        ALLOCATE(N_BLEND_H_CHILD(1:NUM_CHILDREN))                          !<-- Boundary blending width for child H points
        ALLOCATE(N_BLEND_V_CHILD(1:NUM_CHILDREN))                          !<-- Boundary blending width for child V points
        ALLOCATE(INC_FIX(1:NUM_CHILDREN))                                  !<-- See below where INC_FIX is filled
        ALLOCATE(RANK_MOVING_CHILD(1:NUM_CHILDREN))                        !<-- Location of moving nests in list of all children
        ALLOCATE(STATIC_OR_MOVING(1:NUM_CHILDREN))                         !<-- Are the individual children static or moving?
!
!-----------------------------------------------------------------------
!
        child_info_loop: DO N=1,NUM_CHILDREN
!
!-----------------------------------------------------------------------
!***  Initialize to nonsense the newly allocated arrays.
!-----------------------------------------------------------------------
!
          IM_CHILD(N)                =-999
          JM_CHILD(N)                =-999
          PARENT_CHILD_SPACE_RATIO(N)=-999
          CHILD_PARENT_SPACE_RATIO(N)=-999.
          INC_FIX(N)                 =-999
          RANK_MOVING_CHILD(N)       =-999
!
          STATIC_OR_MOVING(N)        ='Static'
!
!-----------------------------------------------------------------------
!***  The parent loads each of its children's configure files.
!-----------------------------------------------------------------------
!
          CF(N)=ESMF_ConfigCreate(rc=RC)
!
          CHILD_ID=MY_CHILDREN_ID(N)
          CONFIG_ID=DOMAIN_ID_TO_RANK(CHILD_ID)
          WRITE(INT_TO_CHAR,FMT)CONFIG_ID
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
          IF(RC/=0)THEN
            WRITE(0,*)' Parent unable to load child configure file '    &
                     ,TRIM(CONFIG_FILE_NAME)                            &
                     ,' in PARENT_CHILD_CPL_INITIALIZE'
            WRITE(0,*)' ABORTING!'
            CALL ESMF_FINALIZE(terminationflag=ESMF_ABORT               &
                              ,rc             =RC)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Extract the children's domain sizes from the configure files.
!-----------------------------------------------------------------------
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
                                      ,value =JM_CHILD(N)               &  !<-- The variable filled (JM of child domain)
                                      ,label ='jm:'                     &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Extract the children's boundary blending widths.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Extract Child Bndry Blending Width"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(N)                     &  !<-- The child's config object
                                      ,value =N_BLEND_H_CHILD(N)        &  !<-- The variable filled (N_BLEND_H of child domain N)
                                      ,label ='lnsh:'                   &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF(N)                     &  !<-- The child's config object
                                      ,value =N_BLEND_V_CHILD(N)        &  !<-- The variable filled (N_BLEND_V of child domain N)
                                      ,label ='lnsv:'                   &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(N_BLEND_V_CHILD(N)>N_BLEND_H_CHILD(N))THEN
            WRITE(0,*)' N_BLEND_V CANNOT EXCEED N_BLEND_H DUE TO PD AVERAGING!!!'
            WRITE(0,*)' ABORTING in child N=',N
            CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
          ENDIF
!
!-----------------------------------------------------------------------
!***  Extract the integer ratio of parent-to-child grid increments.
!-----------------------------------------------------------------------
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
!***  Use this ratio to compute an increment that is needed for
!***  selecting the appropriate nest tasks as mass values on the
!***  nest boundaries are averaged to the V points.  Its values
!***  are based on the nest grid increment distance from the
!***  southernmost V point on a nest's southernmost tasks to
!***  the nearest parent V point to the north.  Fractional values
!***  are increased to the next integer.
!-----------------------------------------------------------------------
!
          DIST_NESTV_SOUTH_TO_PARENTV_SOUTH=                            &
                                     (PARENT_CHILD_SPACE_RATIO(N)-1)*0.5
          INC_FIX(N)=INT(DIST_NESTV_SOUTH_TO_PARENTV_SOUTH+0.9)
!
!-----------------------------------------------------------------------
!***  Which of the children will be moving?  Save their ranking  
!***  in the list of all children.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract the Child's Flag Indicating Movability"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(N)                       &  !<-- The child's config object
                                      ,value =DOMAIN_MOVES                &  !<-- The variable filled (will the child move?)
                                      ,label ='my_domain_moves:'          &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!  
          IF(DOMAIN_MOVES)THEN                                               !<-- If true then child N moves.
!
            NUM_MOVING_CHILDREN=NUM_MOVING_CHILDREN+1                        !<-- Add up the # of moving children.
            RANK_MOVING_CHILD(NUM_MOVING_CHILDREN)=N                         !<-- Location in list of children of those who move.
            STATIC_OR_MOVING(N)='Moving'                                     !<-- Child N moves
!
          ENDIF 
!
!-----------------------------------------------------------------------
!***  We do not allow moving parents to have static children for now.
!-----------------------------------------------------------------------
!
#ifdef ESMF_3
          IF(MY_DOMAIN_MOVES==ESMF_TRUE)THEN
#else
          IF(MY_DOMAIN_MOVES)THEN
#endif
            IF(NUM_MOVING_CHILDREN/=NUM_CHILDREN)THEN
              WRITE(0,*)' You have specified a moving parent with'      &
                       ,' static children.  This is not allowed. '
              WRITE(0,*)' Moving parents can have only moving children.'
              WRITE(0,*)' Aborting!'
              CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
            ENDIF
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO child_info_loop
!
!-----------------------------------------------------------------------
!***  Allocate arrays/pointers needed by the parents to compute
!***  boundary data for their children.  The routine is called only
!***  by parents.
!-----------------------------------------------------------------------
!
        CALL PARENT_CHILD_INTERP_SETUP(NUM_CHILDREN                     &
                                      ,MY_CHILDREN_ID                   &
                                      ,IM_CHILD                         &
                                      ,JM_CHILD                         &
                                      ,FTASKS_DOMAIN                    &
                                      ,CTASK_LIMITS                     &
                                      ,N_BLEND_H_CHILD                  &
                                      ,N_BLEND_V_CHILD                  &
                                      ,CF                               &
                                      ,ITS,ITE,JTS,JTE                  &
                                      ,IDS,IDE,JDS,JDE )
!
!-----------------------------------------------------------------------
!***  We now compute various indices and weights needed by the parents
!***  to compute boundary data for their children.  It is here that
!***  location-dependent interpolation information is determined 
!***  regarding the parent and nests.  Again only parents call this 
!***  routine.
!-----------------------------------------------------------------------
!
        DO N=1,NUM_CHILDREN
          CALL PREPARE_NEST_INTERP_FACTORS(N)
        ENDDO
!
!-----------------------------------------------------------------------
!***  Allocate the pointers that will hold all of the interpolated
!***  boundary data for the child tasks if the parent task contains
!***  child boundary points on the four sides.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The data pointer within the CHILD_BOUND_* arrays will
!***  hold the boundary data of each child tasks' boundary data in a 
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
!***          allocated and instead are simply pointed into the
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
        ALLOCATE(CHILD_BOUND_H_SOUTH(1:NUM_CHILDREN,1:2))                   !<-- 1-D bndry data string for child tasks with Sbndry H points
        ALLOCATE(CHILD_BOUND_V_SOUTH(1:NUM_CHILDREN,1:2))                   !<-- 1-D bndry data string for child tasks with Sbndry V points
        ALLOCATE(WORDS_BOUND_H_SOUTH(1:NUM_CHILDREN))                       !<-- # of words in Sbndry H point 1-D data string
        ALLOCATE(WORDS_BOUND_V_SOUTH(1:NUM_CHILDREN))                       !<-- # of words in Sbndry V point 1-D data string
!
        ALLOCATE(PD_B_SOUTH(1:NUM_CHILDREN))                                !<-- South boundary PD on H points
        ALLOCATE(PD_B_SOUTH_V(1:NUM_CHILDREN))                              !<-- South boundary PD on V points
        ALLOCATE(T_B_SOUTH (1:NUM_CHILDREN))                                !<-- South boundary temperature
        ALLOCATE(Q_B_SOUTH (1:NUM_CHILDREN))                                !<-- South boundary specific humidity
        ALLOCATE(CW_B_SOUTH(1:NUM_CHILDREN))                                !<-- South boundary cloud condensate
        ALLOCATE(U_B_SOUTH (1:NUM_CHILDREN))                                !<-- South boundary U wind component
        ALLOCATE(V_B_SOUTH (1:NUM_CHILDREN))                                !<-- South boundary V wind component
!
        ALLOCATE(CHILD_BOUND_H_NORTH(1:NUM_CHILDREN,1:2))                   !<-- 1-D bndry data string for child tasks with Nbndry H points
        ALLOCATE(CHILD_BOUND_V_NORTH(1:NUM_CHILDREN,1:2))                   !<-- 1-D bndry data string for child tasks with Nbndry V points
        ALLOCATE(WORDS_BOUND_H_NORTH(1:NUM_CHILDREN))                       !<-- # of words in Nbndry H point 1-D data string
        ALLOCATE(WORDS_BOUND_V_NORTH(1:NUM_CHILDREN))                       !<-- # of words in Nbndry V point 1-D data string
!
        ALLOCATE(PD_B_NORTH(1:NUM_CHILDREN))                                !<-- North boundary PD on H points
        ALLOCATE(PD_B_NORTH_V(1:NUM_CHILDREN))                              !<-- North boundary PD on V points
        ALLOCATE(T_B_NORTH (1:NUM_CHILDREN))                                !<-- North boundary temperature
        ALLOCATE(Q_B_NORTH (1:NUM_CHILDREN))                                !<-- North boundary specific humidity
        ALLOCATE(CW_B_NORTH(1:NUM_CHILDREN))                                !<-- North boundary cloud condensate
        ALLOCATE(U_B_NORTH (1:NUM_CHILDREN))                                !<-- North boundary U wind component
        ALLOCATE(V_B_NORTH (1:NUM_CHILDREN))                                !<-- North boundary V wind component
!
        ALLOCATE(CHILD_BOUND_H_WEST(1:NUM_CHILDREN,1:2))                    !<-- 1-D bndry data string for child tasks with Wbndry H points
        ALLOCATE(CHILD_BOUND_V_WEST(1:NUM_CHILDREN,1:2))                    !<-- 1-D bndry data string for child tasks with Wbndry V points
        ALLOCATE(WORDS_BOUND_H_WEST(1:NUM_CHILDREN))                        !<-- # of words in Wbndry H point 1-D data string
        ALLOCATE(WORDS_BOUND_V_WEST(1:NUM_CHILDREN))                        !<-- # of words in Wbndry V point 1-D data string
!
        ALLOCATE(PD_B_WEST(1:NUM_CHILDREN))                                 !<-- West boundary PD on H points
        ALLOCATE(PD_B_WEST_V(1:NUM_CHILDREN))                               !<-- West boundary PD on V points
        ALLOCATE(T_B_WEST (1:NUM_CHILDREN))                                 !<-- West boundary temperature
        ALLOCATE(Q_B_WEST (1:NUM_CHILDREN))                                 !<-- West boundary specific humidity
        ALLOCATE(CW_B_WEST(1:NUM_CHILDREN))                                 !<-- West boundary cloud condensate
        ALLOCATE(U_B_WEST (1:NUM_CHILDREN))                                 !<-- West boundary U wind component
        ALLOCATE(V_B_WEST (1:NUM_CHILDREN))                                 !<-- West boundary V wind component
!
        ALLOCATE(CHILD_BOUND_H_EAST(1:NUM_CHILDREN,1:2))                    !<-- 1-D bndry data string for child tasks with Ebndry H points
        ALLOCATE(CHILD_BOUND_V_EAST(1:NUM_CHILDREN,1:2))                    !<-- 1-D bndry data string for child tasks with Ebndry V points
        ALLOCATE(WORDS_BOUND_H_EAST(1:NUM_CHILDREN))                        !<-- # of words in Ebndry H point 1-D data string
        ALLOCATE(WORDS_BOUND_V_EAST(1:NUM_CHILDREN))                        !<-- # of words in Ebndry V point 1-D data string
!
        ALLOCATE(PD_B_EAST(1:NUM_CHILDREN))                                 !<-- East boundary PD on H points
        ALLOCATE(PD_B_EAST_V(1:NUM_CHILDREN))                               !<-- East boundary PD on V points
        ALLOCATE(T_B_EAST (1:NUM_CHILDREN))                                 !<-- East boundary temperature
        ALLOCATE(Q_B_EAST (1:NUM_CHILDREN))                                 !<-- East boundary specific humidity
        ALLOCATE(CW_B_EAST(1:NUM_CHILDREN))                                 !<-- East boundary cloud condensate
        ALLOCATE(U_B_EAST (1:NUM_CHILDREN))                                 !<-- East boundary U wind component
        ALLOCATE(V_B_EAST (1:NUM_CHILDREN))                                 !<-- East boundary V wind component
!
        DO NN=1,2
          DO N=1,NUM_CHILDREN
            CHILD_BOUND_H_SOUTH(N,NN)%TASKS=>NULL()
            CHILD_BOUND_H_NORTH(N,NN)%TASKS=>NULL()
            CHILD_BOUND_H_WEST(N,NN)%TASKS=>NULL()
            CHILD_BOUND_H_EAST(N,NN)%TASKS=>NULL()
          ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!***  Allocate dummy subcomponents for the working pointers for (N,2)
!***  which correspond to values used when parents must send BC data
!***  to a nest immediately after it has moved, i.e., for that nest's
!***  current time and not for its future time.  Pointers for (N,1)
!***  will always be allocated ahead of deallocation since they are
!***  continually needed for the parent's sending BC data to all the
!***  nests from their future.
!***  In the normal sequence these working pointers are deallocated
!***  each time a nest moves so that they can be reallocated properly
!***  for the given association of parent and nest tasks.  Therefore
!***  they must be allocated already for the deallocations that take
!***  place with each nest's first move.
!-----------------------------------------------------------------------
!
!     write(0,*)' P-C Init before allocate CHILD_BOUND_H_SOUTH(N,2)%TASKS'
!     if(associated(CHILD_BOUND_H_SOUTH(2,1)%TASKS))then
!       write(0,*)' P-C Init CHILD_BOUND_H_SOUTH(2,1) is associated here'
!     else
!       write(0,*)' P-C Init CHILD_BOUND_H_SOUTH(2,1) is not associated here'
!     endif
        DO N=1,NUM_CHILDREN
!
          ALLOCATE(CHILD_BOUND_H_SOUTH(N,2)%TASKS(1))      
          ALLOCATE(CHILD_BOUND_V_SOUTH(N,2)%TASKS(1))      
          CHILD_BOUND_H_SOUTH(N,2)%TASKS(1)%DATA=>NULL()
          CHILD_BOUND_V_SOUTH(N,2)%TASKS(1)%DATA=>NULL()
!
          ALLOCATE(CHILD_BOUND_H_NORTH(N,2)%TASKS(1))      
          ALLOCATE(CHILD_BOUND_V_NORTH(N,2)%TASKS(1))      
          CHILD_BOUND_H_NORTH(N,2)%TASKS(1)%DATA=>NULL()
          CHILD_BOUND_V_NORTH(N,2)%TASKS(1)%DATA=>NULL()
!
          ALLOCATE(CHILD_BOUND_H_WEST(N,2)%TASKS(1))      
          ALLOCATE(CHILD_BOUND_V_WEST(N,2)%TASKS(1))      
          CHILD_BOUND_H_WEST(N,2)%TASKS(1)%DATA=>NULL()
          CHILD_BOUND_V_WEST(N,2)%TASKS(1)%DATA=>NULL()
!
          ALLOCATE(CHILD_BOUND_H_EAST(N,2)%TASKS(1))      
          ALLOCATE(CHILD_BOUND_V_EAST(N,2)%TASKS(1))      
          CHILD_BOUND_H_EAST(N,2)%TASKS(1)%DATA=>NULL()
          CHILD_BOUND_V_EAST(N,2)%TASKS(1)%DATA=>NULL()
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Allocate logical flags indicating if parent task holds any
!***  child boundary points for the purpose of sending that data
!***  to pertinent child tasks.
!-----------------------------------------------------------------------
!
        ALLOCATE(SEND_CHILD_DATA(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  Allocate the handles to be used by parent tasks when they ISend
!***  data directly to the appropiate child boundary tasks.  The 2nd
!***  dimension is 2 because these handles are used in two different
!***  and essentially independent situations.  The first is when the
!***  parents send their children the usual boundary updates from the
!***  children's future so that the children can compute time tendencies
!***  for their integration through the next parent timestep.  The
!***  second is when parents send their moving children the same set
!***  of boundary values that they will receive at one of their
!***  current timesteps immediately after they move to a new location.
!***  Those values will serve as the time N values in the subsequent
!***  tendency computations for variable X: [X(N+1)-X(N)]/DT(parent)
!***  Note that while the 2nd dimension of all children is 2, the
!***  Handles' subcomponents associated with that index's value of 2
!***  will be allocated and used only for moving nests.
!-----------------------------------------------------------------------
!
        ALLOCATE(HANDLE_H_SOUTH(1:NUM_CHILDREN,1:2))
        ALLOCATE(HANDLE_H_NORTH(1:NUM_CHILDREN,1:2))
        ALLOCATE(HANDLE_H_WEST (1:NUM_CHILDREN,1:2))
        ALLOCATE(HANDLE_H_EAST (1:NUM_CHILDREN,1:2))
!
        ALLOCATE(HANDLE_V_SOUTH(1:NUM_CHILDREN,1:2))
        ALLOCATE(HANDLE_V_NORTH(1:NUM_CHILDREN,1:2))
        ALLOCATE(HANDLE_V_WEST (1:NUM_CHILDREN,1:2))
        ALLOCATE(HANDLE_V_EAST (1:NUM_CHILDREN,1:2))
!
        DO N=1,NUM_CHILDREN
!
          ALLOCATE(HANDLE_H_SOUTH(N,2)%NTASKS_TO_RECV(1))
          ALLOCATE(HANDLE_V_SOUTH(N,2)%NTASKS_TO_RECV(1))
          ALLOCATE(HANDLE_H_NORTH(N,2)%NTASKS_TO_RECV(1))
          ALLOCATE(HANDLE_V_NORTH(N,2)%NTASKS_TO_RECV(1))
          ALLOCATE(HANDLE_H_WEST(N,2)%NTASKS_TO_RECV(1))
          ALLOCATE(HANDLE_V_WEST(N,2)%NTASKS_TO_RECV(1))
          ALLOCATE(HANDLE_H_EAST(N,2)%NTASKS_TO_RECV(1))
          ALLOCATE(HANDLE_V_EAST(N,2)%NTASKS_TO_RECV(1))
!
          HANDLE_H_SOUTH(N,2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
          HANDLE_V_SOUTH(N,2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
          HANDLE_H_NORTH(N,2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
          HANDLE_V_NORTH(N,2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
          HANDLE_H_WEST(N,2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
          HANDLE_V_WEST(N,2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
          HANDLE_H_EAST(N,2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
          HANDLE_V_EAST(N,2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Point unallocated working pointers of parent's interpolated data
!***  into the allocated composite pointer holding all boundary data
!***  to be sent to each child from their future.
!-----------------------------------------------------------------------
!
        DO N=1,NUM_CHILDREN
          CALL POINT_INTERP_DATA_TO_MEMORY(N,'Future')
        ENDDO
!
!-----------------------------------------------------------------------
!
!!!     DEALLOCATE(CF)
!
!-----------------------------------------------------------------------
!***  Allocate the pointers that will hold the Surface Geopotential
!***  of child tasks on each side of the child boundaries.  The child
!***  tasks of static nests will send that data to the appropriate
!***  parent tasks.
!-----------------------------------------------------------------------
!
        ALLOCATE(FIS_CHILD_SOUTH(1:NUM_CHILDREN))
        ALLOCATE(FIS_CHILD_NORTH(1:NUM_CHILDREN))
        ALLOCATE(FIS_CHILD_WEST (1:NUM_CHILDREN))
        ALLOCATE(FIS_CHILD_EAST (1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  Allocate an array of logical flags the parent will use for its
!***  moving children to know when they want to move.
!***  Also allocate the composite data object that will hold all of
!***  the update data the parent sends to its moving children and
!***  the associated Handles for the ISends.
!-----------------------------------------------------------------------
!
        IF(NUM_MOVING_CHILDREN>0)THEN
          ALLOCATE(MOVE_FLAG(1:NUM_MOVING_CHILDREN))
          ALLOCATE(HANDLE_BC_UPDATE(1:NUM_MOVING_CHILDREN))
          ALLOCATE(HANDLE_TIMESTEP(1:NUM_MOVING_CHILDREN))
          ALLOCATE(HANDLE_MOVE_DATA(1:NUM_MOVING_CHILDREN))
!
          DO N=1,NUM_MOVING_CHILDREN
            MOVE_FLAG(N)       =.FALSE.
            HANDLE_BC_UPDATE(N)=MPI_REQUEST_NULL
            HANDLE_TIMESTEP(N) =MPI_REQUEST_NULL
!
            N_MOVING=RANK_MOVING_CHILD(N)
            NUM_CHILD_TASKS=FTASKS_DOMAIN(MY_CHILDREN_ID(N_MOVING))
            ALLOCATE(HANDLE_MOVE_DATA(N)%NTASKS_TO_RECV(1:NUM_CHILD_TASKS))
!
            DO NN=1,NUM_CHILD_TASKS
              HANDLE_MOVE_DATA(N)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
            ENDDO
          ENDDO
!
!!!       ALLOCATE(TASK_UPDATE_SPECS_H(1:NUM_MOVING_CHILDREN))
!!!       ALLOCATE(TASK_UPDATE_SPECS_V(1:NUM_MOVING_CHILDREN))
          ALLOCATE(TASK_UPDATE_SPECS(1:NUM_MOVING_CHILDREN))
          ALLOCATE(MOVING_CHILD_UPDATE(1:NUM_MOVING_CHILDREN))
!
          DO N=1,NUM_MOVING_CHILDREN
            TASK_UPDATE_SPECS(N)%TASK_ID=>NULL()
            TASK_UPDATE_SPECS(N)%NUM_PTS_UPDATE_HZ=>NULL()
            TASK_UPDATE_SPECS(N)%NEXT_LINK=>NULL()
!
            MOVING_CHILD_UPDATE(N)%TASKS=>NULL()
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Parents send their children the integration limits of all
!***  forecast tasks on the parent domain.
!-----------------------------------------------------------------------
!
        IF(MYPE==0)THEN
          DO N=1,NUM_CHILDREN
!
            CALL MPI_SEND(LOCAL_ISTART                                  &  !<-- Starting I's of fcst tasks on parent domain
                         ,FTASKS_DOMAIN(MY_DOMAIN_ID)                   &  !<-- # of fcst tasks on parent domain
                         ,MPI_INTEGER                                   &  !<-- Indices are integers
                         ,0                                             &  !<-- Send to each child's task 0
                         ,10001                                         &  !<-- MYPE tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator to child N
                         ,IERR )
!
            CALL MPI_SEND(LOCAL_JSTART                                  &  !<-- Starting J's of fcst tasks on parent domain
                         ,FTASKS_DOMAIN(MY_DOMAIN_ID)                   &  !<-- # of fcst tasks on parent domain
                         ,MPI_INTEGER                                   &  !<-- Indices are integers
                         ,0                                             &  !<-- Send to each child's task 0
                         ,10002                                         &  !<-- MYPE tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator to child N
                         ,IERR )
!
            CALL MPI_SEND(LOCAL_IEND                                    &  !<-- Ending I's of fcst tasks on parent domain
                         ,FTASKS_DOMAIN(MY_DOMAIN_ID)                   &  !<-- # of fcst tasks on parent domain
                         ,MPI_INTEGER                                   &  !<-- Indices are integers
                         ,0                                             &  !<-- Send to each child's task 0
                         ,10003                                         &  !<-- MYPE tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator to child N
                         ,IERR )
!
            CALL MPI_SEND(LOCAL_JEND                                    &  !<-- Ending J's of fcst tasks on parent domain
                         ,FTASKS_DOMAIN(MY_DOMAIN_ID)                   &  !<-- # of fcst tasks on parent domain
                         ,MPI_INTEGER                                   &  !<-- Indices are integers
                         ,0                                             &  !<-- Send to each child's task 0
                         ,10004                                         &  !<-- MYPE tag
                         ,COMM_TO_MY_CHILDREN(N)                        &  !<-- MPI communicator to child N
                         ,IERR )
!
          ENDDO
        ENDIF
!
!-----------------------------------------------------------------------
!***  The parents need to send their children some key information
!***  regarding the association of child boundary tasks with parent
!***  tasks so the children know how to receive boundary data from
!***  their parents.
!-----------------------------------------------------------------------
!
        DO N=1,NUM_CHILDREN
!
          CALL PARENT_SENDS_CHILD_DATA_LIMITS(N)
!
!-----------------------------------------------------------------------
!***  The parent receives the child's boundary topography.  This is
!***  needed to maintain hydrostatic balance when parent data is
!***  interpolated to child boundaries where the terrain is different.
!-----------------------------------------------------------------------
!
          CALL PARENT_RECVS_CHILD_TOPO(N)
!
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDIF parent_block
!
      DEALLOCATE(DOMAIN_ID_TO_RANK)
!
      IF(ASSOCIATED(LOCAL_ISTART))DEALLOCATE(LOCAL_ISTART)
      IF(ASSOCIATED(LOCAL_IEND  ))DEALLOCATE(LOCAL_IEND  )
      IF(ASSOCIATED(LOCAL_JSTART))DEALLOCATE(LOCAL_JSTART)
      IF(ASSOCIATED(LOCAL_JEND  ))DEALLOCATE(LOCAL_JEND  )
!
!-----------------------------------------------------------------------
!***  Task 0 on each child receives the integration limits of their
!***  parents' forecast tasks then broadcasts that information to
!***  the remaining child tasks.
!***  Also each child task needs to allocate the derived type that will
!***  hold: (i) Which parent task(s) will send boundary data to it;
!***  (ii) The grid index limits on the child boundary covered by 
!***  the parent task's data the child will receive.  Then it receives 
!***  those pieces of information from the parent.
!-----------------------------------------------------------------------
!
      IF(COMM_TO_MY_PARENT/=-999)THEN                                      !<-- Select the children.
!
        ID_DOM=ID_PARENTS(MY_DOMAIN_ID)                                    !<-- Domain ID of this child's parent
!
!-----------------------------------------------------------------------
!
        NUM_TASKS_PARENT=FTASKS_DOMAIN(ID_DOM)
!
        ALLOCATE(ITS_PARENT(0:NUM_TASKS_PARENT-1))
        ALLOCATE(ITE_PARENT(0:NUM_TASKS_PARENT-1))
        ALLOCATE(JTS_PARENT(0:NUM_TASKS_PARENT-1))
        ALLOCATE(JTE_PARENT(0:NUM_TASKS_PARENT-1))
!
        IF(MYPE==0)THEN
!
          CALL MPI_RECV(ITS_PARENT                                      &  !<-- Starting I on each parent forecast task's subdomain
                       ,NUM_TASKS_PARENT                                &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Indices are integers
                       ,0                                               &  !<-- Receive from this parent task
                       ,10001                                           &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(JTS_PARENT                                      &  !<-- Starting J on each parent forecast task's subdomain
                       ,NUM_TASKS_PARENT                                &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Indices are integers
                       ,0                                               &  !<-- Receive from this parent task
                       ,10002                                           &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(ITE_PARENT                                      &  !<-- Ending I on each parent forecast task's subdomain
                       ,NUM_TASKS_PARENT                                &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Indices are integers
                       ,0                                               &  !<-- Receive from this parent task
                       ,10003                                           &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(JTE_PARENT                                      &  !<-- Ending J on each parent forecast task's subdomain
                       ,NUM_TASKS_PARENT                                &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Indices are integers
                       ,0                                               &  !<-- Receive from this parent task
                       ,10004                                           &  !<-- MPI tag
                       ,COMM_TO_MY_PARENT                               &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
        ENDIF
!
        CALL MPI_BCAST(ITS_PARENT                                       &
                      ,NUM_TASKS_PARENT                                 &
                      ,MPI_INTEGER                                      &
                      ,0                                                &
                      ,MPI_COMM_COMP                                    &
                      ,IERR )
!
        CALL MPI_BCAST(JTS_PARENT                                       &
                      ,NUM_TASKS_PARENT                                 &
                      ,MPI_INTEGER                                      &
                      ,0                                                &
                      ,MPI_COMM_COMP                                    &
                      ,IERR )
!
        CALL MPI_BCAST(ITE_PARENT                                       &
                      ,NUM_TASKS_PARENT                                 &
                      ,MPI_INTEGER                                      &
                      ,0                                                &
                      ,MPI_COMM_COMP                                    &
                      ,IERR )
!
        CALL MPI_BCAST(JTE_PARENT                                       &
                      ,NUM_TASKS_PARENT                                 &
                      ,MPI_INTEGER                                      &
                      ,0                                                &
                      ,MPI_COMM_COMP                                    &
                      ,IERR )
!
!-----------------------------------------------------------------------
!
        ALLOCATE(PARENT_TASK(1:FTASKS_DOMAIN(ID_DOM)))                     !<-- Dimensioned as # of fcst tasks on domain of parent.
!
        CALL CHILD_RECVS_CHILD_DATA_LIMITS(EXP_STATE)                      !<-- Recv specs of new parent/child task associations
!
!-----------------------------------------------------------------------
!***  All the children send to their parents their boundary
!***  topography so that the parents can properly balance the data
!***  generated for the children's boundaries.  For moving nests
!***  these are only initial values that will change when the nests
!***  move.
!-----------------------------------------------------------------------
!
        CALL CHILD_SENDS_TOPO_TO_PARENT(IMP_STATE) 
!
!-----------------------------------------------------------------------
!
#ifdef ESMF_3
        I_AM_A_NEST=ESMF_TRUE
#else
        I_AM_A_NEST=.TRUE.
#endif
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(COMM_TO_MY_PARENT==-999)THEN
#ifdef ESMF_3
        I_AM_A_NEST=ESMF_FALSE
#else
        I_AM_A_NEST=.FALSE.
#endif
      ENDIF
!
!-----------------------------------------------------------------------
!***  Everyone loads the coupler export state with the flag indicating
!***  whether or not they are a nest.  Nests load the flag indicating
!***  if they move or not.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Nest Flag into the Coupler Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The Parent_Child coupler export state
                            ,name ='I-Am-A-Nest Flag'                   &  !<-- The name of the flag
                            ,value=I_AM_A_NEST                          &  !<-- The nest flag
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      IF(I_AM_A_NEST==ESMF_TRUE)THEN
#else
      IF(I_AM_A_NEST)THEN
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert Moving Nest Flag into the Coupler Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child coupler export state
                              ,name ='MY_DOMAIN_MOVES'                  &  !<-- The name of the flag
                              ,value=MY_DOMAIN_MOVES                    &  !<-- The moving nest flag
                              ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        I_WANT_TO_MOVE=.FALSE.                                             !<-- Initialize the nest 'move' flag
        MOVE_FLAG_SENT=.FALSE.                                             !<-- Initialize the flag for ISending the nest move flag
        RECVD_MOVE_TIMESTEP=.FALSE.                                        !<-- Initialize the flag for receiving move timestep
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Now take care of several issues related to moving nests.
!***  All moving nests and their parents must participate.
!-----------------------------------------------------------------------
!
      parents_and_moving: IF(NUM_MOVING_CHILDREN>0                      &  !<-- This is a parent of moving nests.
                                    .OR.                                &  !
#ifdef ESMF_3
                             MY_DOMAIN_MOVES==ESMF_TRUE)THEN               !<-- This is a moving nest.
#else
                             MY_DOMAIN_MOVES)THEN                          !<-- This is a moving nest.
#endif
!
!-----------------------------------------------------------------------
!***  Extract the Bundle with the 2-D and 3-D arrays of
!***  Dynamics and Physics internal state variables needed for 
!***  updating any nests that are moving.  Since the eventual update
!***  of moving nest data will be done via looping through the Fields
!***  in the Bundles we need to know how many Fields there are.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Bundles for Updates of Moving Nests in P-C Init"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE                        &  !<-- The Parent-Child coupler import state
                          ,itemname   ='Move_Bundle H'                  &  !<-- Name of Bundle of internal state arrays to update
                          ,fieldbundle=MOVE_BUNDLE_H                    &  !<-- The H-point ESMF Bundle 
                          ,rc         =RC)
!
        CALL ESMF_StateGet(state      =IMP_STATE                        &  !<-- The Parent-Child coupler import state
                          ,itemname   ='Move_Bundle V'                  &  !<-- Name of Bundle of internal state arrays to update
                          ,fieldbundle=MOVE_BUNDLE_V                    &  !<-- The V-point ESMF Bundle 
                          ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="How many Fields in the H Bundle?"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520r
          CALL ESMF_FieldBundleGet(fieldBundle=MOVE_BUNDLE_H            &  !<-- The ESMF Bundle of H update arrays for moving nests
                                  ,fieldCount =NUM_FIELDS_MOVE          &  !<-- # of Fields in the Bundle
                                  ,rc         =RC)
#else
        CALL ESMF_FieldBundleGet(bundle    =MOVE_BUNDLE_H               &  !<-- The ESMF Bundle of H update arrays for moving nests
                                ,fieldCount=NUM_FIELDS_MOVE             &  !<-- # of Fields in the Bundle
                                ,rc        =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Count the number of 2-D and 3-D Fields.  Those numbers will be
!***  needed to know how many points are updated on moving nest tasks.
!-----------------------------------------------------------------------
!
        NUM_FIELDS_MOVE_2D_H_I=0
        NUM_FIELDS_MOVE_2D_X_I=0
        NUM_FIELDS_MOVE_2D_H_R=0
        NUM_FIELDS_MOVE_2D_X_R=0
        NUM_FIELDS_MOVE_3D_H=0
        NUM_LEVELS_MOVE_3D_H=0
!
!-----------------------------------------------------------------------
!
        DO N_FIELD=1,NUM_FIELDS_MOVE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Fields from H Move Bundle for Counting"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520r
          CALL ESMF_FieldBundleGet(fieldBundle=MOVE_BUNDLE_H            &  !<-- Bundle holding the H arrays for move updates
                                  ,fieldIndex =N_FIELD                  &  !<-- Index of the Field in the Bundle
                                  ,field      =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                  ,rc         =RC)
#else
          CALL ESMF_FieldBundleGet(bundle    =MOVE_BUNDLE_H             &  !<-- Bundle holding the H arrays for move updates
                                  ,fieldIndex=N_FIELD                   &  !<-- Index of the Field in the Bundle
                                  ,field     =HOLD_FIELD                &  !<-- Field N_FIELD in the Bundle
                                  ,rc        =RC )
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="How Many Dims in H Move Bundle Field?"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field   =HOLD_FIELD                        &  !<-- Field N_FIELD in the Bundle
                            ,dimCount=NUM_DIMS                          &  !<-- Is this Field 2-D or 3-D?
                            ,rc      =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract UPDATE_TYPE from Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(field=HOLD_FIELD                       &  !<-- Get Attribute from this Field
                                ,name ='UPDATE_TYPE'                    &  !<-- Name of the attribute to extract
                                ,value=UPDATE_TYPE_INT                  &  !<-- Value of the Attribute
                                ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(UPDATE_TYPE_INT==1)THEN
            UPDATE_TYPE_CHAR='H'                                           !<-- Ordinary H-pt variable
          ELSEIF(UPDATE_TYPE_INT==2)THEN
            UPDATE_TYPE_CHAR='L'                                           !<-- H-pt land sfc variable
          ELSEIF(UPDATE_TYPE_INT==3)THEN
            UPDATE_TYPE_CHAR='S'                                           !<-- H-pt sea sfc variable
          ELSEIF(UPDATE_TYPE_INT==4)THEN
            UPDATE_TYPE_CHAR='F'                                           !<-- H-pt variable updated from external file
          ELSEIF(UPDATE_TYPE_INT==5)THEN
            UPDATE_TYPE_CHAR='V'                                           !<-- Ordinary V-pt variable
          ENDIF
!
!-----------------------------------------------------------------------
!
          IF(NUM_DIMS==2)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Does the Field Contain Integer or Real Data?"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_FieldGet(field   =HOLD_FIELD                      &  !<-- Field N_FIELD in the Bundle
                              ,typekind=DATATYPE                        &  !<-- Is the data Integer or Real?
                              ,rc      =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(DATATYPE==ESMF_TYPEKIND_I4)THEN
              NUM_FIELDS_MOVE_2D_H_I=NUM_FIELDS_MOVE_2D_H_I+1              !<-- Count ALL 2-D Integer Fields
              IF(UPDATE_TYPE_CHAR=='F')THEN
                NUM_FIELDS_MOVE_2D_X_I=NUM_FIELDS_MOVE_2D_X_I+1            !<-- Count the 2-D Integer variables updated from external files
              ENDIF
!
            ELSE
              NUM_FIELDS_MOVE_2D_H_R=NUM_FIELDS_MOVE_2D_H_R+1              !<-- Count ALL 2-D Real Fields
              IF(UPDATE_TYPE_CHAR=='F')THEN
                NUM_FIELDS_MOVE_2D_X_R=NUM_FIELDS_MOVE_2D_X_R+1            !<-- Count the 2-D Real variables updated from external files
              ENDIF
!
            ENDIF
!
!-----------------------------------------------------------------------
!
          ELSEIF(NUM_DIMS==3)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract 3rd Dimension Limits in 3-D H Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
            CALL ESMF_FieldGet(field      =HOLD_FIELD                   &  !<-- Field N in the Bundle
                              ,localDe    =0                            &
                              ,farray     =DUMMY_3D                     &  !<-- Dummy 3-D array with Field's data
                              ,totalLBound=LIMITS_LO                    &  !<-- Starting index in each dimension
                              ,totalUBound=LIMITS_HI                    &  !<-- Ending index in each dimension
                              ,rc         =RC )
#else
            CALL ESMF_FieldGet(field      =HOLD_FIELD                   &  !<-- Field N in the Bundle
                              ,localDe    =0                            &
                              ,farrayPtr  =DUMMY_3D                     &  !<-- Dummy 3-D array with Field's data
                              ,totalLBound=LIMITS_LO                    &  !<-- Starting index in each dimension
                              ,totalUBound=LIMITS_HI                    &  !<-- Ending index in each dimension
                              ,rc         =RC )
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NUM_FIELDS_MOVE_3D_H=NUM_FIELDS_MOVE_3D_H+1                    !<-- Count the 3-D Real H Fields
!
            NUM_LEVELS_MOVE_3D_H=LIMITS_HI(3)-LIMITS_LO(3)+1            &  !<-- Count the # of 2-D levels in the 3-D H Fields
                                +NUM_LEVELS_MOVE_3D_H
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="How many Fields in the V Bundle?"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520r
          CALL ESMF_FieldBundleGet(fieldBundle=MOVE_BUNDLE_V            &  !<-- The ESMF Bundle of V update arrays for moving nests
                                  ,fieldCount =NUM_FIELDS_MOVE          &  !<-- # of Fields in the Bundle
                                  ,rc         =RC)
#else
        CALL ESMF_FieldBundleGet(bundle    =MOVE_BUNDLE_V               &  !<-- The ESMF Bundle of V update arrays for moving nests
                                ,fieldCount=NUM_FIELDS_MOVE             &  !<-- # of Fields in the Bundle
                                ,rc        =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Count the number of 2-D and 3-D Fields.  Those numbers will be
!***  needed to know how many points are updated on moving nest tasks.
!-----------------------------------------------------------------------
!
        NUM_FIELDS_MOVE_2D_V=0
        NUM_FIELDS_MOVE_3D_V=0
        NUM_LEVELS_MOVE_3D_V=0
!
        DO N_FIELD=1,NUM_FIELDS_MOVE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Fields from V Move Bundle for Counting"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520r
          CALL ESMF_FieldBundleGet(fieldBundle=MOVE_BUNDLE_V            &  !<-- Bundle holding the H arrays for move updates
                                  ,fieldIndex =N_FIELD                  &  !<-- Index of the Field in the Bundle
                                  ,field      =HOLD_FIELD               &  !<-- Field N_FIELD in the Bundle
                                  ,rc         =RC)
#else
          CALL ESMF_FieldBundleGet(bundle    =MOVE_BUNDLE_V             &  !<-- Bundle holding the H arrays for move updates
                                  ,fieldIndex=N_FIELD                   &  !<-- Index of the Field in the Bundle
                                  ,field     =HOLD_FIELD                &  !<-- Field N_FIELD in the Bundle
                                  ,rc        =RC )
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="How Many Dims in V Move Bundle Field?"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field   =HOLD_FIELD                        &  !<-- Field N_FIELD in the Bundle
                            ,dimCount=NUM_DIMS                          &  !<-- Is this Field 2-D or 3-D?
                            ,rc      =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(NUM_DIMS==2)THEN
            NUM_FIELDS_MOVE_2D_V=NUM_FIELDS_MOVE_2D_V+1
!
          ELSEIF(NUM_DIMS==3)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Extract 3rd Dimension Limits in 3-D V Field"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
            CALL ESMF_FieldGet(field      =HOLD_FIELD                   &  !<-- Field N in the Bundle
                              ,localDe    =0                            &
                              ,farray     =DUMMY_3D                     &  !<-- Dummy 3-D array with Field's data
                              ,totalLBound=LIMITS_LO                    &  !<-- Starting index in each dimension
                              ,totalUBound=LIMITS_HI                    &  !<-- Ending index in each dimension
                              ,rc         =RC )
#else
            CALL ESMF_FieldGet(field      =HOLD_FIELD                   &  !<-- Field N in the Bundle
                              ,localDe    =0                            &
                              ,farrayPtr  =DUMMY_3D                     &  !<-- Dummy 3-D array with Field's data
                              ,totalLBound=LIMITS_LO                    &  !<-- Starting index in each dimension
                              ,totalUBound=LIMITS_HI                    &  !<-- Ending index in each dimension
                              ,rc         =RC )
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            NUM_FIELDS_MOVE_3D_V=NUM_FIELDS_MOVE_3D_V+1                    !<-- Count the 3-D V Fields
!
            NUM_LEVELS_MOVE_3D_V=LIMITS_HI(3)-LIMITS_LO(3)+1            &  !<-- Count the # of 2-D levels in the 3-D V Fields
                                +NUM_LEVELS_MOVE_3D_V
          ENDIF
!
        ENDDO
!
!       WRITE(0,*)' P-C Cpl Init '                                      &
!                ,' NUM_FIELDS_MOVE_2D_H_I=',NUM_FIELDS_MOVE_2D_H_I     &
!                ,' NUM_FIELDS_MOVE_2D_H_R=',NUM_FIELDS_MOVE_2D_H_R     &
!                ,' NUM_FIELDS_MOVE_3D_H=',NUM_FIELDS_MOVE_3D_H         &
!                ,' NUM_FIELDS_MOVE_2D_V=',NUM_FIELDS_MOVE_2D_V         &
!                ,' NUM_FIELDS_MOVE_3D_V=',NUM_FIELDS_MOVE_3D_V
!
!-----------------------------------------------------------------------
!***  The moving nests and their parents read in the four configure
!***  variables specifying the number of boundary rows on the nests'
!***  pre-move footprints whose locations will receive update data
!***  from the parent after the nests move.  All moving nests use
!***  the same values and the parent checks to be sure this is true.
!-----------------------------------------------------------------------
!
        parents_with_movers: IF(NUM_MOVING_CHILDREN>0)THEN                  !<-- Parents read moving nests' configure files
!
!-----------------------------------------------------------------------
!
          DO N=1,NUM_MOVING_CHILDREN
            NN=RANK_MOVING_CHILD(N)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Parent-Child Init: Parent Reads NROWS_P_UPD_W"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(N)                   &  !<-- The child's config object
                                        ,value =NROWS_P_UPD_X           &  !<-- # of footprint W bndry rows updated by parent
                                        ,label ='nrows_p_upd_w:'        &  !<-- The configure label
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(N==1)THEN
              NROWS_P_UPD_W=NROWS_P_UPD_X
            ELSE
              IF(NROWS_P_UPD_X/=NROWS_P_UPD_W)THEN
                WRITE(0,*)' Moving nests must have same values for NROWS_P_UPD_W!'
                WRITE(0,*)' Aborting!'
                CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
              ENDIF
            ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Parent-Child Init: Parent Reads NROWS_P_UPD_E"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(N)                   &  !<-- The child's config object
                                        ,value =NROWS_P_UPD_X           &  !<-- # of footprint E bndry rows updated by parent
                                        ,label ='nrows_p_upd_e:'        &  !<-- The configure label
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(N==1)THEN
              NROWS_P_UPD_E=NROWS_P_UPD_X
            ELSE
              IF(NROWS_P_UPD_X/=NROWS_P_UPD_E)THEN
                WRITE(0,*)' Moving nests must have same values for NROWS_P_UPD_E!'
                WRITE(0,*)' Aborting!'
                CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
              ENDIF
            ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Parent-Child Init: Parent Reads NROWS_P_UPD_S"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(N)                   &  !<-- The child's config object
                                        ,value =NROWS_P_UPD_X           &  !<-- # of footprint S bndry rows updated by parent
                                        ,label ='nrows_p_upd_s:'        &  !<-- The configure label
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(N==1)THEN
              NROWS_P_UPD_S=NROWS_P_UPD_X
            ELSE
              IF(NROWS_P_UPD_X/=NROWS_P_UPD_S)THEN
                WRITE(0,*)' Moving nests must have same values for NROWS_P_UPD_S!'
                WRITE(0,*)' Aborting!'
                CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
              ENDIF
            ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Parent-Child Init: Parent Reads NROWS_P_UPD_N"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_ConfigGetAttribute(config=CF(N)                   &  !<-- The child's config object
                                        ,value =NROWS_P_UPD_X           &  !<-- # of footprint N bndry rows updated by parent
                                        ,label ='nrows_p_upd_n:'        &  !<-- The configure label
                                        ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            IF(N==1)THEN
              NROWS_P_UPD_N=NROWS_P_UPD_X
            ELSE
              IF(NROWS_P_UPD_X/=NROWS_P_UPD_N)THEN
                WRITE(0,*)' Moving nests must have same values for NROWS_P_UPD_N!'
                WRITE(0,*)' Aborting!'
                CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
              ENDIF
            ENDIF
!
!-----------------------------------------------------------------------
!
          ENDDO
!
!-----------------------------------------------------------------------
!***  All parents of moving nests will be reading those children's
!***  full resolution topography files that span the entire uppermost
!***  parent.  This will require these parents to know the dimensions
!***  as well as other key aspects of the uppermost parent's grid.
!***  Read the pertinent data from the uppermost parent's configure
!***  file and save what will be needed later.
!-----------------------------------------------------------------------
!
          CF_1=ESMF_ConfigCreate(rc=RC)
!
          CONFIG_FILE_NAME='configure_file_01'                             !<-- Config file name of uppermost parent
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Load Configure Object of Upper Domain in P-C Cpl Init"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigLoadFile(config  =CF_1                        &
                                  ,filename=CONFIG_FILE_NAME            &
                                  ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Base Dimensions of Uppermost Domain"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_1                      &  !<-- The config object
                                      ,value =IM_1                      &  !<-- The variable filled
                                      ,label ='im:'                     &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF_1                      &  !<-- The config object
                                      ,value =JM_1                      &  !<-- The variable filled
                                      ,label ='jm:'                     &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Central Lat/Lon of Uppermost Domain"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_1                      &  !<-- The config object
                                      ,value =TPH0D_1                   &  !<-- The variable filled
                                      ,label ='tph0d:'                  &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF_1                      &  !<-- The config object
                                      ,value =TLM0D_1                   &  !<-- The variable filled
                                      ,label ='tlm0d:'                  &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Southern/Western Boundary of Uppermost Domain"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_1                      &  !<-- The config object
                                      ,value =SBD_1                     &  !<-- The variable filled
                                      ,label ='sbd:'                    &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
          CALL ESMF_ConfigGetAttribute(config=CF_1                      &  !<-- The config object
                                      ,value =WBD_1                     &  !<-- The variable filled
                                      ,label ='wbd:'                    &  !<-- Give this label's value to the previous variable
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
          D_ONE=1.
          D_180=180.
          PI=DACOS(-D_ONE)
          D2R=PI/D_180
!
          TPH0_1=TPH0D_1*D2R                                               !<-- Central geo lat of domain (radians, positive north)
          TLM0_1=TLM0D_1*D2R                                               !<-- Central geo lon of domain (radians, positive east)
          WB_1=WBD_1*D2R                                                   !<-- Rotated lon of west boundary (radians, positive east)
          SB_1=SBD_1*D2R                                                   !<-- Rotated lat of south boundary (radians, positive north)
!
          DPH_1=-2.*SB_1/(JM_1-1)                                          !<-- Uppermost parent's grid increment in J (radians)
          DLM_1=-2.*WB_1/(IM_1-1)                                          !<-- Uppermost parent's grid increment in I (radians)
!
          RECIP_DPH_1=1./DPH_1
          RECIP_DLM_1=1./DLM_1
!
!-----------------------------------------------------------------------
!
        ENDIF parents_with_movers
!
!-----------------------------------------------------------------------
!
#ifdef ESMF_3
        IF(MY_DOMAIN_MOVES==ESMF_TRUE)THEN                                 !<-- Moving nests read their configure files
#else
        IF(MY_DOMAIN_MOVES)THEN                                            !<-- Moving nests read their configure files
#endif
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Nest Reads NROWS_P_UPD_W"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_MINE                   &  !<-- The nest's config object
                                      ,value =NROWS_P_UPD_W             &  !<-- # of footprint W bndry rows updated by parent
                                      ,label ='nrows_p_upd_w:'          &  !<-- The configure label
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Nest Reads NROWS_P_UPD_E"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_MINE                   &  !<-- The nest's config object
                                      ,value =NROWS_P_UPD_E             &  !<-- # of footprint E bndry rows updated by parent
                                      ,label ='nrows_p_upd_e:'          &  !<-- The configure label
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Nest Reads NROWS_P_UPD_S"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_MINE                   &  !<-- The nest's config object
                                      ,value =NROWS_P_UPD_S             &  !<-- # of footprint S bndry rows updated by parent
                                      ,label ='nrows_p_upd_s:'          &  !<-- The configure label
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Parent-Child Init: Nest Reads NROWS_P_UPD_N"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF_MINE                   &  !<-- The nest's config object
                                      ,value =NROWS_P_UPD_N             &  !<-- # of footprint N bndry rows updated by parent
                                      ,label ='nrows_p_upd_n:'          &  !<-- The configure label
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF parents_and_moving
!
!-----------------------------------------------------------------------
!***  Additional setup for parents of moving nests and/or for moving
!***  parents with any children at all.
!-----------------------------------------------------------------------
!
      nests_move: IF(NUM_MOVING_CHILDREN>0)THEN
!
!-----------------------------------------------------------------------
!***  By 'moving nest' we mean any domain that moves across the earth
!***  and not those domains that move within their parent's domain.
!***  This would therefore include static children inside moving parents
!***  however that arrangement is not allowed at present.  That setup
!***  would require full updates of the child domain following the
!***  parent's shift since the child moved with resepct to the earth
!***  and atmosphere.  However a moving child in a moving parent will
!***  stay in place when its parent shifts and thus that child domain
!***  needs no updating at all following the parent shift.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Parents with moving nests need to know those nests' topography
!***  at the nests' own resolutions for the hydrostatic adjustment
!***  that must take place when the parents interpolate their data
!***  to moving nest grid points.  For the sake of generality all
!***  of those nest-resolution datasets must span the domain of the
!***  uppermost parent which is the true maximum range of any nest's
!***  motion.
!
!***  So each parent with moving nests must:
!***   (1) Know how many different space resolutions its moving 
!***       children use;
!***   (2) Associate each resolution with the appropriate moving
!***       child using the nest-to-uppermost parent space ratio
!***       that the user specified in each moving nest's configure 
!***       file;
!***   (3) Have each of its forecast tasks read in its own piece of
!***       each different resolution of topography data needed by
!***       all of its moving children.
!
!***  If a parent domain moves then it must refill its task subdomains
!***  with the topography of its moving children each time it (the
!***  parent) shifts its position.  That is handled in subroutine
!***  PARENT_CHILD_CPL_RUN_RECV.
!-----------------------------------------------------------------------
!
        ALLOCATE(M_NEST_RATIO(1:NUM_MOVING_CHILDREN))                      !<-- Associate moving nests with list of different space ratios
        ALLOCATE(LIST_OF_RATIOS(1:NUM_MOVING_CHILDREN))                    !<-- Keep a list of the different space ratios
        ALLOCATE(LINK_MRANK_RATIO(1:NUM_MOVING_CHILDREN))                  !<-- Which different space ratio for each moving child
!
        NN=0
        I_RATIO=0
        KOUNT_RATIOS_MN=0                                                  !<-- Count the different resolutions of moving nests
!
!-----------------------------------------------------------------------
!
!!!     DO N=1,NUM_MOVING_CHILDREN
        DO N=1,NUM_CHILDREN
!
          IF(STATIC_OR_MOVING(N)=='Static')CYCLE
          NN=NN+1
!
!-----------------------------------------------------------------------
!
          LIST_OF_RATIOS(NN)=0         
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Moving Child's Sfc File Ratio"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_ConfigGetAttribute(config=CF(N)                     &  !<-- Child N's config object
                                      ,value =SFC_FILE_RATIO            &  !<-- Save the configure value with the following label.
                                      ,label ='ratio_sfc_files:'        &  !<-- Ratio of upper parent's grid increment to this nest's
                                      ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!  
          M_NEST_RATIO(NN)=SFC_FILE_RATIO                                  !<-- Moving child NN uses topography file with this ratio/ID
!
          IF(NN==1)THEN
            KOUNT_RATIOS_MN=1                                              !<-- Begin counting the ratios
            LIST_OF_RATIOS(1)=SFC_FILE_RATIO                               !<-- The 1st sfc file ratio is that of the 1st moving nest
            LINK_MRANK_RATIO(1)=1                                          !<-- 1st moving nest uses 1st sfc file ratio
!
          ELSE
            FOUND=.FALSE.                                                  
            DO KR=1,KOUNT_RATIOS_MN
              IF(SFC_FILE_RATIO==LIST_OF_RATIOS(KR))THEN
                FOUND=.TRUE.
                LINK_MRANK_RATIO(NN)=KR                                    !<-- Moving nest NN uses existing KR'th sfc file ratio 
                EXIT
              ENDIF
            ENDDO
!
            IF(.NOT.FOUND)THEN
              KOUNT_RATIOS_MN=KOUNT_RATIOS_MN+1                            !<-- Increment the counter of different ratios
              LIST_OF_RATIOS(KOUNT_RATIOS_MN)=SFC_FILE_RATIO               !<-- Save the new ratio in the list of different ratios
              LINK_MRANK_RATIO(NN)=KOUNT_RATIOS_MN                         !<-- Moving child NN uses this rank in list of all different ratios
            ENDIF
!
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        ALLOCATE(NEST_FIS_ON_PARENT(1:KOUNT_RATIOS_MN))
        ALLOCATE(NEST_FIS_V_ON_PARENT(1:KOUNT_RATIOS_MN))
!
        DO N=1,KOUNT_RATIOS_MN
          NEST_FIS_ON_PARENT(N)%DATA=>NULL()
          NEST_FIS_V_ON_PARENT(N)%DATA=>NULL()
        ENDDO
!
        ALLOCATE(NEST_FIS_ON_PARENT_BNDS(1:KOUNT_RATIOS_MN))
!
!-----------------------------------------------------------------------
!***  Now fill the parent's data objects that hold the nest-resolution
!***  topography at child H points and child V points.
!-----------------------------------------------------------------------
!
        CALL PARENT_READS_MOVING_CHILD_TOPO(NUM_MOVING_CHILDREN         &
                                           ,LINK_MRANK_RATIO            &
                                           ,LIST_OF_RATIOS              &
                                           ,M_NEST_RATIO                &
                                           ,KOUNT_RATIOS_MN             &
                                           ,IM_1,JM_1                   &
                                           ,TPH0_1,TLM0_1               &
                                           ,SB_1,WB_1                   &
                                           ,RECIP_DPH_1,RECIP_DLM_1     &
                                           ,GLAT,GLON                   &
                                           ,NEST_FIS_ON_PARENT_BNDS     &
                                           ,NEST_FIS_ON_PARENT          &
                                           ,NEST_FIS_V_ON_PARENT        &
                                           ,IDS,IDE,IMS,IME,ITS,ITE     &
                                           ,JDS,JDE,JMS,JME,JTS,JTE)
!
!-----------------------------------------------------------------------
!
      ENDIF nests_move
!
!-----------------------------------------------------------------------
!
      IF(ALLOCATED(CF))DEALLOCATE(CF)
!
!-----------------------------------------------------------------------
!
!------------
!***  Timers
!------------
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
      moving_nest_bookkeep_tim  =0.
      moving_nest_update_tim    =0.
      parent_bookkeep_moving_tim=0.
      parent_update_moving_tim  =0.
      t0_recv_move_tim          =0.
!
!-----------------------------------------------------------------------
!
      IF(RC_CPL_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)"PARENT_CHILD_CPL INITIALIZE STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"PARENT_CHILD_CPL INITIALIZE STEP FAILED"
      ENDIF
!
      RC_FINAL=RC_CPL_INIT
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
!***  Run the coupler step where children receive data from parents.
!***  This is phase 1 of the coupler RUn step since it occurs at the 
!***  beginning of the timesteps.  The parents send the data in 
!***  phase 2 at the end of the timesteps.  Only child tasks enter
!***  this routine.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_CplComp) :: CPL_COMP                                       !<-- The Parent-Child Coupler Component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Coupler's Import State
                         ,EXP_STATE                                        !<-- The Coupler's Export State
!
      TYPE(ESMF_Clock) :: CLOCK                                            !<-- The ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_FINAL
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT),SAVE :: HANDLE_MOVE_FLAG=MPI_REQUEST_NULL      &
                                ,HANDLE_NEXT_MOVE=MPI_REQUEST_NULL      &
                                ,HANDLE_SW_CORNER=MPI_REQUEST_NULL
!
      INTEGER(kind=KINT),SAVE :: I_SW_PARENT_NEW                        &
                                ,J_SW_PARENT_NEW                        &
                                ,LAST_STEP_MOVED                        &
                                ,NEXT_MOVE_TIMESTEP_PARENT
!
      INTEGER(kind=KINT) :: BC_UPDATE_FLAG,N,NTIMESTEP,NTAG0
!
      INTEGER(kind=KINT) :: IERR,IRTN,RC,RC_CPL_RUN
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,SAVE ::               &
                                                   HANDLE_PARENT_SHIFT
!
      LOGICAL(kind=KLOG) :: MY_MOVE_FLAG
!
      TYPE(INTERIOR_DATA_FROM_PARENT),DIMENSION(1:4) :: SEND_TASK          !<-- Specifics about interior data from sending parent tasks
!
      LOGICAL(kind=KLOG) :: MOVE_TIME_IS_PRESENT                        &
                           ,PARENT_SHIFT_IS_PRESENT
!
#ifdef ESMF_3
      TYPE(ESMF_Logical) :: MOVE_NOW
#else
      LOGICAL(kind=KLOG) :: MOVE_NOW
#endif
!
      integer(kind=kint),save :: kount_moves=0
      integer(kind=kint),dimension(8) :: values
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
!***  Initialize the error signal variables.
!-----------------------------------------------------------------------
!
      RC        =ESMF_SUCCESS
      RC_FINAL  =ESMF_SUCCESS
      RC_CPL_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  What is the current timestep on this nest's Clock?
!-----------------------------------------------------------------------
!
      CALL ESMF_ClockGet(clock       =CLOCK                             &
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- The current timestep of this child (ESMF)
                        ,rc          =RC)
!
      NTIMESTEP=NTIMESTEP_ESMF                                             !<-- The current timestep of this child (integer)
!
!-----------------------------------------------------------------------
!***  The child is now at the beginning of a timestep that coincides
!***  with the beginning of a parent timestep.  If the child has
!***  determined that it wants its domain to move then it must send
!***  its parent a signal to that effect along with its new position
!***  on the parent's grid.  A parent must always run ahead of its
!***  children therefore the parent will receive 'move' messages at the
!***  end of a timestep in the future.  It is at that point in time
!***  that the parent will generate new internal data for the children
!***  who want to move as well as new starting boundary data for the
!***  child's new location and thus it is at that point in time that 
!***  the children will be able to execute their moves.  The parent
!***  must send back to the child that (the parent's) timestep so the
!***  child will know when to receive the new data.  The 'move' signal
!***  (true) is sent to the parent only when the child has determined
!***  it wants to move.
!-----------------------------------------------------------------------
!
#ifdef ESMF_3
      moving_children_a: IF(MY_DOMAIN_MOVES==ESMF_TRUE)THEN                !<-- Select the moving nests
#else
      moving_children_a: IF(MY_DOMAIN_MOVES)THEN                           !<-- Select the moving nests
#endif
!        
!----------------------------------------------------------
!***  If the child sent an affirmative 'move' signal at
!***  the beginning of the previous parent timestep then
!***  it now receives from the parent the parent's timestep
!***  at which the new data is to be received.  We wait to
!***  use the parent timestep for the move here rather than
!***  immediately after the child sent its 'move' signal
!***  to allow the child to proceed and not wait for a
!***  parent who might not be ready to send its response
!***  back to the child immediately.
!----------------------------------------------------------
!
        new_move_time: IF(I_WANT_TO_MOVE                   &  !<-- Does this nest want to move but doesn't know when to?
                               .AND.                       &  
                          .NOT.RECVD_MOVE_TIMESTEP)THEN      
!
          nest_task_0_a: IF(MYPE==0)THEN                                   !<-- Lead forecast task on this moving nest
!
            CALL MPI_IPROBE(0                                           &  !<-- The message is sent by moving child N's fcst task 0.
                           ,TIME_TAG                                    &  !<-- Tag associated with nest N's move timestep
                           ,COMM_TO_MY_PARENT                           &  !<-- MPI communicator between this nest and its parent
                           ,MOVE_TIME_IS_PRESENT                        &  !<-- Is the nest's new move time now available?
                           ,JSTAT                                       &
                           ,IERR)
!
!---------------------------------------------------------
!***  If NEXT_MOVE_TIMESTEP_PARENT from the parent is now
!***  present then use it to obtain the child timestep
!***  at which new data from the parent for the move is
!***  available.
!---------------------------------------------------------
!
            IF(MOVE_TIME_IS_PRESENT)THEN
              CALL MPI_RECV(NEXT_MOVE_TIMESTEP_PARENT                   &  !<-- Recv the message to clear the nest's ISEND
                           ,1                                           &  !<-- # of words in the message
                           ,MPI_INTEGER                                 &  !<-- The timestep is an integer
                           ,0                                           &  !<-- Local rank of the parent task sending the word
                           ,TIME_TAG                                    &  !<-- Arbitrary tag used for this data exchange
                           ,COMM_TO_MY_PARENT                           &  !<-- MPI communicator betweenthis nest and its parent
                           ,JSTAT                                       &
                           ,IERR )
!
              NEXT_MOVE_TIMESTEP=(NEXT_MOVE_TIMESTEP_PARENT+1)*         &  !<-- Child will move at beginning of this timestep
                                  TIME_RATIO_MY_PARENT
              RECVD_MOVE_TIMESTEP=.TRUE.
!
            ELSE
              RECVD_MOVE_TIMESTEP=.FALSE.
            ENDIF
!
          ENDIF nest_task_0_a
!
!---------------------------------------------------------
!***  Child task 0 needs to inform the other child tasks
!***  if the new move timestep has been received from
!***  the parent.  If it has then that new timestep is
!***  broadcast to the nest tasks.  That will be the 
!***  point in time at which BC data for the new nest
!***  location is received along with internal shift 
!***  data from the parent.
!---------------------------------------------------------
!
          CALL MPI_BCAST(RECVD_MOVE_TIMESTEP                            &  !<-- Has the new move timestep been received?
                        ,1                                              &  !<-- The timestep is one word
                        ,MPI_LOGICAL                                    &  !<-- The timestep is type Logical
                        ,0                                              &  !<-- Broadcast from nest forecast task 0
                        ,MPI_COMM_COMP                                  &  !<-- MPI communicator for this nest's forecast tasks
                        ,IRTN)
!
          IF(RECVD_MOVE_TIMESTEP)THEN
            CALL MPI_BCAST(NEXT_MOVE_TIMESTEP                           &  !<-- The next timestep to move 
                          ,1                                            &  !<-- The timestep is one word
                          ,MPI_INTEGER                                  &  !<-- The timestep is type Integer
                          ,0                                            &  !<-- Broadcast from nest forecast task 0
                          ,MPI_COMM_COMP                                &  !<-- MPI communicator for this nest's forecast tasks
                          ,IRTN)
!
            LAST_STEP_MOVED=NEXT_MOVE_TIMESTEP                             !<-- Reset to this newest move even if in the future
          ENDIF
!
        ENDIF new_move_time
!
!---------------------------------------------------------
!***  If this is now the point in time at which the
!***  parent learned the child was wanting to move,
!***  and is thus the time at which the parent prepared
!***  internal and boundary data for the child's new
!***  position, then the child moves now.  It will receive
!***  the data prepared for it by its parent for the 
!***  child's new position and will also shift all of its 
!***  internal data that remains within the lateral extent
!***  of the child's grid position prior to the move.
!---------------------------------------------------------
!
#ifdef ESMF_3
        MOVE_NOW=ESMF_FALSE
#else
        MOVE_NOW=.FALSE.
#endif
!
        the_child_moves: IF(NTIMESTEP==NEXT_MOVE_TIMESTEP)THEN
!
#ifdef ESMF_3
          MOVE_NOW=ESMF_TRUE                                               !<-- Yes, the child moves at beginning of this timestep
#else
          MOVE_NOW=.TRUE.                                                  !<-- Yes, the child moves at beginning of this timestep
#endif
!
          I_WANT_TO_MOVE=.FALSE.                                           !<-- Reset the 'move' flag
          MOVE_FLAG_SENT=.FALSE.                                           !<-- Reset the flag for ISending the move flag
          RECVD_MOVE_TIMESTEP=.FALSE.                                      !<-- Reset flag for receiving the move timestep
!
          I_SHIFT_CHILD=(I_SW_PARENT_NEW-I_SW_PARENT_CURRENT)           &  !<-- The shift in I in this domain's space
                        *SPACE_RATIO_MY_PARENT
          J_SHIFT_CHILD=(J_SW_PARENT_NEW-J_SW_PARENT_CURRENT)           &  !<-- The shift in J in this domain's space
                        *SPACE_RATIO_MY_PARENT
!
!-----------------------------------------------------------------------
!***  As stated at its outset the fundamental purpose of this routine
!***  is for a nested domain to recv update data from its parent.
!***  In the special case in which this child domain is also a parent
!***  then one parent-specific action needs to be taken here.  The
!***  domain is just now moving at this point in time and that will
!***  change its tasks' relationships with the boundary tasks of its
!***  (moving) children who remain in place with respect to the earth
!***  and atmosphere.  Therefore this domain as a parent must signal
!***  its children now that this is when their inter-task relationships
!***  change.  This will force the children to wait to recv the new
!***  parent-child task specifications before they execute their normal
!***  recvs of BC data updates from the future at the end of this
!***  routine.
!-----------------------------------------------------------------------
!
          IF(NUM_CHILDREN>0.AND.MYPE==0)THEN
!
            IF(.NOT.ALLOCATED(HANDLE_PARENT_SHIFT))THEN
              ALLOCATE(HANDLE_PARENT_SHIFT(1:NUM_CHILDREN))
              DO N=1,NUM_CHILDREN
                HANDLE_PARENT_SHIFT(N)=MPI_REQUEST_NULL
             ENDDO
            ENDIF
!
            PARENT_SHIFT(1)=I_SHIFT_CHILD                                    !<-- Parent's I shift in its space
            PARENT_SHIFT(2)=J_SHIFT_CHILD                                    !<-- Parent's J shift in its space
!
            DO N=1,NUM_CHILDREN
!
              CALL MPI_WAIT(HANDLE_PARENT_SHIFT(N)                        &  !<-- Handle for ISend of parent's shift
                           ,JSTAT                                         &  !<-- MPI status
                           ,IERR)
!
              NTAG0=PARENT_SHIFT_TAG+NTIMESTEP                               !<-- Unique timestep-dependent MPI tag
!
              CALL MPI_ISEND(PARENT_SHIFT                                 &  !<-- Send parent's shift to all its children
                            ,2                                            &  !<-- There are 2 words in the message
                            ,MPI_INTEGER                                  &  !<-- The shift increments are integers
                            ,0                                            &  !<-- Signal sent to all lead child tasks
                            ,NTAG0                                        &  !<-- Tag used for this data exchange
                            ,COMM_TO_MY_CHILDREN(N)                       &  !<-- MPI communicator between this parent and its children
                            ,HANDLE_PARENT_SHIFT(N)                       &  !<-- Handle for this ISend to children
                            ,IERR )
!
            ENDDO
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  That completes the action of this domain as a parent.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Deallocate this moving nest's working arrays/pointers whose 
!***  dimensions are functions of moving nests' positions.  They will 
!***  be reallocated with dimensions appropriate for the new positions.
!***  For static nests the same nest arrays/pointers never become
!***  invalid and thus are not deallocated/reallocated.
!-----------------------------------------------------------------------
!
          CALL DEALLOC_WORK_CHILDREN                                       !<-- Reset this child's working pointers for new location
!
!-----------------------------------------------------------------------
!***  Each child boundary task receives small information packets
!***  from the parent tasks that cover them.  That information
!***  can change with each move and includes the identities of 
!***  those parent tasks that will be sending boundary data
!***  updates along with the index limits on the child task of
!***  that boundary data from each parent task.
!-----------------------------------------------------------------------
!
          CALL CHILD_RECVS_CHILD_DATA_LIMITS(EXP_STATE)                    !<-- Recv specs of new parent/child task associations
!
!-----------------------------------------------------------------------
!***  Receive standard boundary data update from the parent valid for 
!***  the current timestep but now at the child's new location.  This
!***  data will be for time N in the boundary tendency computation:
!***  dX/dt = [ X(N+1) - X(N) ] / DT_parent
!-----------------------------------------------------------------------
!
          CALL RECV_NEST_BC_DATA('Current')                                !<-- Recv parent data for new nest boundary after move
!
!-----------------------------------------------------------------------
!***  Receive update data for all interior points on the nest domain
!***  that have moved outside of the nest's pre-move footprint.  Those
!***  points can only be updated by the parent.  The index limits of
!***  the parent update regions on the nest tasks are identical for
!***  H and V points therefore the nest needs to call its bookkeeping
!***  only once.
!-----------------------------------------------------------------------
!
          btim=timef()
          CALL MOVING_NEST_BOOKKEEPING(I_SHIFT_CHILD                    &
                                      ,J_SHIFT_CHILD                    &
                                      ,I_SW_PARENT_NEW                  &
                                      ,J_SW_PARENT_NEW                  &
                                      ,NUM_TASKS_PARENT                 &
                                      ,INPES_PARENT                     &
                                      ,ITS_PARENT                       &
                                      ,ITE_PARENT                       &
                                      ,JTS_PARENT                       &
                                      ,JTE_PARENT                       &
                                      ,SPACE_RATIO_MY_PARENT            &
                                      ,NROWS_P_UPD_W                    &
                                      ,NROWS_P_UPD_E                    &
                                      ,NROWS_P_UPD_S                    &
                                      ,NROWS_P_UPD_N                    &
                                      ,SEND_TASK                        &
                                      ,ITS,ITE,JTS,JTE                  &
                                      ,IMS,IME,JMS,JME                  &
                                      ,IDS,IDE,JDS,JDE                  &
                                        )
!
          btim=timef()
          CALL MOVING_NEST_RECV_DATA(COMM_TO_MY_PARENT                  &
                                    ,NTIMESTEP                          &
                                    ,NUM_FIELDS_MOVE_2D_H_I             &
                                    ,NUM_FIELDS_MOVE_2D_X_I             &
                                    ,NUM_FIELDS_MOVE_2D_H_R             &
                                    ,NUM_FIELDS_MOVE_2D_X_R             &
                                    ,NUM_LEVELS_MOVE_3D_H               &
                                    ,NUM_FIELDS_MOVE_2D_V               &
                                    ,NUM_LEVELS_MOVE_3D_V               &
                                    ,SEND_TASK                          &
                                    ,EXP_STATE                          &
                                        )
          moving_nest_update_tim=moving_nest_update_tim+(timef()-btim)
!
          I_SW_PARENT_CURRENT=I_SW_PARENT_NEW
          J_SW_PARENT_CURRENT=J_SW_PARENT_NEW
!
        ENDIF the_child_moves
!
!-----------------------------------------------------------------------
!***  Load the Attribute into the coupler export state indicating that 
!***  the child is or is not moving at this timestep.  This information 
!***  will be used in the transfer of the new data from the Parent-Child
!***  coupler export state to the DOMAIN import state to the Dynamics
!***  import state.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Child Inserts Move Flag into Cpl Export State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent-Child coupler export state
                              ,name ='MOVE_NOW'                         &  !<-- Name of the attribute to insert
                              ,value=MOVE_NOW                           &  !<-- Is the child moving right now?
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------------------------
!***  The moving nests load their motion in I and J
!***  into the coupler export state so it can be sent 
!***  to the DOMAIN component where intra-task and
!***  inter-task updates on the nest domains can be
!***  carried out.  Those updates have nothing to do
!***  with their parents.
!---------------------------------------------------------
!
#ifdef ESMF_3
        IF(MOVE_NOW==ESMF_TRUE)THEN
#else
        IF(MOVE_NOW)THEN
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Insert I_SHIFT/J_SHIFT in RUN_RECV"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeSet(state=EXP_STATE                        &  !<-- The Parent-Child coupler export state
                                ,name ='I_SHIFT'                        &  !<-- Insert Attribute with this name
                                ,value=I_SHIFT_CHILD                    &  !<-- Motion of nest in I on its grid
                                ,rc   =RC )
!
          CALL ESMF_AttributeSet(state=EXP_STATE                        &  !<-- The Parent-Child coupler export state
                                ,name ='J_SHIFT'                        &  !<-- Insert Attribute with this name
                                ,value=J_SHIFT_CHILD                    &  !<-- Motion of nest in J on its grid
                                ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDIF
!
!---------------------------------------------------------
!
      ENDIF moving_children_a
!
!-----------------------------------------------------------------------
!***  If this is a (moving) child of a moving parent then it needs
!***  to watch for the appearance of a signal from the parent
!***  indicating that it (the parent) has just moved.  If it has
!***  then the nest receives the parent's shift and the new specifics
!***  of the parent task and nest task associations.
!-----------------------------------------------------------------------
!
      moving_parent: IF(MY_PARENT_MOVES)THEN
!
!-----------------------------------------------------------------------
!
        IF(MYPE==0)THEN
          NTAG0=PARENT_SHIFT_TAG+NTIMESTEP/TIME_RATIO_MY_PARENT            !<-- Unique timestep-dependent MPI tag
!
          CALL MPI_IPROBE(0                                             &  !<-- The message is sent by moving parent's fcst task 0.
                         ,NTAG0                                         &  !<-- Tag associated with parent's shift
                         ,COMM_TO_MY_PARENT                             &  !<-- MPI communicator between this nest and its parent
                         ,PARENT_SHIFT_IS_PRESENT                       &  !<-- Is the parent's shift now available?
                         ,JSTAT                                         &
                         ,IERR)
!
          IF(PARENT_SHIFT_IS_PRESENT)THEN
            CALL MPI_RECV(PARENT_SHIFT                                  &  !<-- Recv the parent's shift
                         ,2                                             &  !<-- # of words in the message
                         ,MPI_INTEGER                                   &  !<-- The shifts in I and J are integers
                         ,0                                             &  !<-- Local rank of the parent task sending the word
                         ,NTAG0                                         &  !<-- Tag used for this data exchange
                         ,COMM_TO_MY_PARENT                             &  !<-- MPI communicator betweenthis nest and its parent
                         ,JSTAT                                         &
                         ,IERR )
          ENDIF
!
        ENDIF
!
        CALL MPI_BCAST(PARENT_SHIFT_IS_PRESENT                          &  !<-- Has the new move timestep been received?
                      ,1                                                &  !<-- The timestep is one word
                      ,MPI_LOGICAL                                      &  !<-- The timestep is type Logical
                      ,0                                                &  !<-- Broadcast from nest forecast task 0
                      ,MPI_COMM_COMP                                    &  !<-- MPI communicator for this nest's forecast tasks
                      ,IRTN)
!
        IF(PARENT_SHIFT_IS_PRESENT)THEN
!
          CALL MPI_BCAST(PARENT_SHIFT                                   &  !<-- Broadcast the parent shift
                        ,2                                              &  !<-- # of words in the message
                        ,MPI_INTEGER                                    &  !<-- The shifts in I and J are integers
                        ,0                                              &  !<-- Broadcast from nest forecast task 0
                        ,MPI_COMM_COMP                                  &  !<-- MPI communicator for this nest's forecast tasks
                        ,IRTN)
!
          I_SW_PARENT_CURRENT=I_SW_PARENT_CURRENT-PARENT_SHIFT(1)          !<-- Parent I,J of SW corner of this child's domain
          J_SW_PARENT_CURRENT=J_SW_PARENT_CURRENT-PARENT_SHIFT(2)          !<--
!
          I_SW_PARENT_NEW=I_SW_PARENT_NEW-PARENT_SHIFT(1)                  !<-- Shifted parent I,J of SW corner of this child's domain
          J_SW_PARENT_NEW=J_SW_PARENT_NEW-PARENT_SHIFT(2)                  !<--
!
          CALL DEALLOC_WORK_CHILDREN                                       !<-- Reset this child's working pointers for 'new' location
!
          CALL CHILD_RECVS_CHILD_DATA_LIMITS(EXP_STATE)                    !<-- Parent/child bndry task associations
!
        ENDIF
!
!---------------------------------------------------------
!
      ENDIF moving_parent
!
!---------------------------------------------------------
!***  All children receive boundary data from their
!***  parents from one parent timestep in the future
!***  which will be put into the Parent-Child coupler
!***  export state on its way to the Dynamics where it
!***  will be used to compute boundary tendencies through
!***  the next NN timesteps.  NN is the number of child
!***  timesteps per parent timestep.  This data will be
!***  for time N+1 in the boundary tendency computation:
!***  dX/dt = [ X(N+1) - X(N) ] / DT_parent
!
!***  Note that if this is a timestep at which a child
!***  has just moved, the child's location-dependent
!***  working pointers for the boundary data were already
!***  reset for the new location in the IF block for
!***  NTIMESTEP==NEXT_MOVE_TIMESTEP above.
!---------------------------------------------------------
!
      CALL RECV_NEST_BC_DATA('Future')
!
!---------------------------------------------------------
!***  Moving nests now decide if they want to move and
!***  then tell their parents if the answer is yes.
!***  The move will be executed in one of the next two
!***  parent timesteps in this routine immediately above.
!***  The move will be executed after the parent has 
!***  learned that the nest wants to move and has then
!***  prepared BC and internal shift data valid at the
!***  timestep that the parent received the message.
!***  If the nest is still waiting to shift following
!***  the most recent decision to do so then it does not
!***  enter the storm motion routine.
!---------------------------------------------------------
!
#ifdef ESMF_3
      moving_children_b: IF(MY_DOMAIN_MOVES==ESMF_TRUE)THEN                !<-- Select the moving nests
#else
      moving_children_b: IF(MY_DOMAIN_MOVES)THEN                           !<-- Select the moving nests
#endif
!
!---------------------------------------------------------
!***  First receive the flag from this moving nest's
!***  parent indicating that all of the parent's moving
!***  children have been sent their BC update data.
!---------------------------------------------------------
!
        IF(MYPE==0)THEN
!
          CALL MPI_RECV(BC_UPDATE_FLAG                                  &  !<-- All moving siblings have received BC updates
                       ,1                                               &  !<-- # of words sent
                       ,MPI_INTEGER                                     &  !<-- The flag is an integer.
                       ,0                                               &  !<-- Receiving from parent task 0
                       ,MOVING_BC_TAG                                   &  !<-- Tag associated with moving nests' recvs of BC data
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between this nest and its parent
                       ,JSTAT                                           &
                       ,IERR )
!
        ENDIF
!
!---------------------------------------------------------
!
        IF(.NOT.I_WANT_TO_MOVE)THEN
!
          CALL COMPUTE_STORM_MOTION(NTIMESTEP                           &
                                   ,LAST_STEP_MOVED                     &
                                   ,DT_DOMAIN(MY_DOMAIN_ID)             &
                                   ,NUM_PES_FCST                        &
                                   ,FIS                                 &
                                   ,PD                                  &
                                   ,PINT                                &
                                   ,T                                   &
                                   ,Q                                   &
                                   ,CW                                  &
                                   ,U                                   &
                                   ,V                                   &
                                   ,DSG2                                &
                                   ,PDSG1                               &
                                   ,DXH                                 &
                                   ,DYH                                 &
                                   ,SM                                  &
                                   ,I_SW_PARENT_CURRENT                 &
                                   ,J_SW_PARENT_CURRENT                 &
                                   ,I_WANT_TO_MOVE                      &
                                   ,I_SW_PARENT_NEW                     &
                                   ,J_SW_PARENT_NEW )
!
        ENDIF
!
!-----------------------------------------------------------------------
!
        nest_task_0_b: IF(MYPE==0)THEN                                     !<-- Lead forecast task on this moving nest
!
          IF(I_WANT_TO_MOVE.AND..NOT.MOVE_FLAG_SENT)THEN                   !<-- Nest wants to move; send flag if not already sent
!
            CALL MPI_WAIT(HANDLE_MOVE_FLAG                              &  !<-- Handle for ISend of child's move flag to parent
                         ,JSTAT                                         &  !<-- MPI status
                         ,IERR)
!
            MY_MOVE_FLAG=I_WANT_TO_MOVE                                    !<-- Safe to reload the buffer variable
!
            CALL MPI_ISEND(MY_MOVE_FLAG                                 &  !<-- Signal to parent:  Does this child want to move?
                          ,1                                            &  !<-- There is 1 word in the flag
                          ,MPI_LOGICAL                                  &  !<-- Signal is type Logical
                          ,0                                            &  !<-- Signal sent to parent task 0
                          ,MOVE_TAG                                     &  !<-- Arbitrary tag used for this data exchange
                          ,COMM_TO_MY_PARENT                            &  !<-- MPI communicator between this child and its parent
                          ,HANDLE_MOVE_FLAG                             &  !<-- Handle for ISend to parent
                          ,IERR )
!
            MOVE_FLAG_SENT=.TRUE.
!
!-----------------------------------------------------------------------
!***  If child wants to move, send the new location to parent.
!-----------------------------------------------------------------------
!
            CALL MPI_WAIT(HANDLE_SW_CORNER                              &  !<-- Handle for ISend of child's new position to parent
                         ,JSTAT                                         &  !<-- MPI status
                         ,IERR)
!
            IJ_SHIFT_CHILD(1)=I_SW_PARENT_NEW-I_SW_PARENT_CURRENT          !<-- Safe to reload shift buffers
            IJ_SHIFT_CHILD(2)=J_SW_PARENT_NEW-J_SW_PARENT_CURRENT          !<--
!
            CALL MPI_ISEND(IJ_SHIFT_CHILD                               &  !<-- Change in parent I,J of this child's SW corner after move
                          ,2                                            &  !<-- Sending 2 words
                          ,MPI_INTEGER                                  &  !<-- Indices are integers
                          ,0                                            &  !<-- Indices sent to parent task 0
                          ,9001                                         &  !<-- Arbitrary tag used for this data exchange
                          ,COMM_TO_MY_PARENT                            &  !<-- MPI communicator between this child and its parent
                          ,HANDLE_SW_CORNER                             &  !<-- Handle for this ISend to parent
                          ,IERR )
!
!-----------------------------------------------------------------------
!***  After the parent receives the nest's move flag the parent will
!***  send back the value of its own timestep at that moment.  That
!***  will be the point in time at which the nest will actually be
!***  able to execute its move since that is when the parent will
!***  have prepared BC data for the nest's new position as well as
!***  internal shift data for updating nest points outisde of the
!***  pre-move footprint.
!-----------------------------------------------------------------------
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDIF nest_task_0_b
!
!-----------------------------------------------------------------------
!
      ENDIF moving_children_b
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE RECV_NEST_BC_DATA(TIME_FLAG)
!
!-----------------------------------------------------------------------
!***  A nest receives boundary data from its parent so that
!***  it can compute boundary tendencies for its integration.
!***  This is an internal subroutine to PARENT_CHILD_CPL_RUN_RECV.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      CHARACTER(len=*),INTENT(IN) :: TIME_FLAG                             !<-- BC data valid for current or future timestep?
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IERR,N,NP_H,NP_V,NTAG
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Each child task that holds part of the domain boundary receives
!***  data from the parent tasks that cover its boundary points.  
!***  This occurs for two different situations:
!***   (1) All nests receive boundary data from their parents that was
!***       sent from one parent timestep in the future.  That allows
!***       each nest to compute boundary value tendencies that are
!***       applied through the next NN child timesteps of integration
!***       where NN is the number of child timesteps within each parent
!***       timestep.
!***   (2) Immediately after a moving nest moves, it needs new boundary
!***       values for that current time at the new location.  The
!***       structure of that boundary data is the same as in (1) so
!***       it is received from the parents in the same way.  However
!***       that data then needs to be stored as the values for current
!***       parent timestep N where the boundary tendency for variable X
!***       is [X(N+1)-X(N)]/DT(parent).  The values from the future 
!***       timestep N+1 will subsequently be received as usual.
!
!***  Thus we simply need to know whether the incoming data is for the
!***  future time (#1 above) or for the current time (#2 above).  That
!***  information is given by this subroutine's input argument.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!***  Now for each side of the nests' boundaries:
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
      cpl1_prelim_tim=cpl1_prelim_tim+(timef()-btim0)
!
!--------------------
!***  South H Points
!--------------------
!
      btim0=timef()
!
      NP_H=NUM_PARENT_TASKS_SENDING_H%SOUTH                                !<-- # of parent tasks sending south boundary H data
!
      IF(NP_H>0)THEN
!
        NTAG=NTIMESTEP+101                                                 !<-- Add 101 to obtain a unique south H tag
!
        DO N=1,NP_H                                                        !<-- Loop over each parent task sending Sboundary H data
!         call date_and_time(values=values)     
!         write(0,123)n,parent_task(n)%south_h%id_source,values(5),values(6),values(7),values(8)
! 123     format(' Ready to recv South_H from parent task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!     write(0,*)' ready to recv South_H from parent task #',n,' id=',parent_task(n)%south_h%id_source,' from ',time_flag
!     write(0,*)' # of words=',parent_task(n)%south_h%length,' ntag=',ntag,' comm=',comm_to_my_parent

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
!     write(0,*)' recvd South_H from parent task #',n,' id=',parent_task(n)%south_h%id_source,' for ',time_flag
!
!-----------------------------------------------------------------------
!
          btim=timef()
          CALL CHILD_DATA_FROM_STRING(length_data=PARENT_TASK(N)%SOUTH_H%LENGTH     &  !<-- Length of parent datastring
                                     ,datastring =PARENT_TASK(N)%SOUTH_H%STRING     &  !<-- Parent datastring of child task bndry segment 
                                     ,ilim_lo    =INDX_MIN_H%SOUTH                  &  !<-- Lower I limit of child's segment of boundary
                                     ,ilim_hi    =INDX_MAX_H%SOUTH                  &  !<-- Upper I limit of child's segment of boundary
                                     ,jlim_lo    =1                                 &  !<-- Lower J limit of child's segment of boundary
                                     ,jlim_hi    =N_BLEND_H                         &  !<-- Upper J limit of child's segment of boundary
                                     ,i_start    =PARENT_TASK(N)%SOUTH_H%INDX_START &  !<-- Child's segment Istart on each parent task
                                     ,i_end      =PARENT_TASK(N)%SOUTH_H%INDX_END   &  !<-- Child's segment Iend on each parent task
                                     ,j_start    =1                                 &  !<-- Child's segment Jstart on each parent task
                                     ,j_end      =N_BLEND_H                         &  !<-- Child's segment Jend on each parent task
                                     ,pdb        =PDB_S                             &  !<-- Child's 1-D segment of PD on this side of bndry
                                     ,tb         =TB_S                              &  !<-- Child's 1-D segment of T on this side of bndry
                                     ,qb         =QB_S                              &  !<-- Child's 1-D segment of Q on this side of bndry
                                     ,cwb        =CWB_S )                              !<-- Child's 1-D segment of CW on this side of bndry
!!!                                  ,len_2d     =LENGTH_BND_SEG_H%SOUTH )             !<-- Length of the 1-D segment in I and J
!
          cpl1_south_h_undo_tim=cpl1_south_h_undo_tim+(timef()-btim)
!         write(0,*)' south_h recv time=',((timef()-btim))*1.e-3
!
!         call date_and_time(values=values)     
!         write(0,224)n,parent_task(n)%south_h%id_source,values(5),values(6),values(7),values(8)
! 224     format(' After South_H Recv data_from_string task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
        ENDDO
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL EXPORT_CHILD_BOUNDARY(pdb         =PDB_S                   &  !<-- Child's 1-D segment of PD on this side of bndry
                                  ,tb          =TB_S                    &  !<-- Child's 1-D segment of T on this side of bndry
                                  ,qb          =QB_S                    &  !<-- Child's 1-D segment of Q on this side of bndry
                                  ,cwb         =CWB_S                   &  !<-- Child's 1-D segment of CW on this side of bndry
                                  ,ilim_lo     =INDX_MIN_H%SOUTH        &  !<-- Lower I limit of child's segment of boundary
                                  ,ilim_hi     =INDX_MAX_H%SOUTH        &  !<-- Upper I limit of child's segment of boundary
                                  ,jlim_lo     =1                       &  !<-- Lower J limit of child's segment of boundary
                                  ,jlim_hi     =N_BLEND_H               &  !<-- Upper J limit of child's segment of boundary
                                  ,data_name   ='SOUTH_H_'//TIME_FLAG   &  !<-- Name attached to the combined exported data
                                  ,data_exp    =BOUND_1D_SOUTH_H        &  !<-- Combined boundary segment H data for child task
                                  ,export_state=EXP_STATE )                !<-- The Parent-Child Coupler export state
!
          cpl1_south_h_exp_tim=cpl1_south_h_exp_tim+(timef()-btim)
!         call date_and_time(values=values)     
!         write(0,225)values(5),values(6),values(7),values(8)
! 225     format(' After South_H Recv export at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!       write(0,*)' Recvd South_H for ',time_flag
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
      NP_V=NUM_PARENT_TASKS_SENDING_V%SOUTH                                !<-- # of parent tasks sending south boundary V data
!     write(0,*)' PARENT_CHILD_CPL_RUN_RECV south np_v=',np_v
!
      IF(NP_V>0)THEN
        NTAG=NTIMESTEP+102                                                 !<-- Add 102 to obtain a unique south H tag
!
        DO N=1,NP_V                                                        !<-- Loop over each parent task sending Sboundary V data
!         call date_and_time(values=values)     
!         write(0,125)n,parent_task(n)%south_v%id_source,values(5),values(6),values(7),values(8)
! 125     format(' Ready to recv South_V from parent task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
!     write(0,*)' RECV_NEST_BC_DATA for V south ready to recv ',PARENT_TASK(N)%SOUTH_V%LENGTH,' words from parent task #',n &
!              ,' with ID=',PARENT_TASK(N)%SOUTH_V%ID_SOURCE
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
!-----------------------------------------------------------------------
!
          btim=timef()
!     write(0,*)' RECV_NEST_BC_DATA for V south recvd ',PARENT_TASK(N)%SOUTH_V%LENGTH,' words from parent task #',n &
!              ,' with ID=',PARENT_TASK(N)%SOUTH_V%ID_SOURCE
          CALL CHILD_DATA_FROM_STRING(length_data=PARENT_TASK(N)%SOUTH_V%LENGTH     &  !<-- Length of parent datastring
                                     ,datastring =PARENT_TASK(N)%SOUTH_V%STRING     &  !<-- Parent datastring of child task bndry segment
                                     ,ilim_lo    =INDX_MIN_V%SOUTH                  &  !<-- Lower I limit of child's segment of boundary
                                     ,ilim_hi    =INDX_MAX_V%SOUTH                  &  !<-- Upper I limit of child's segment of boundary
                                     ,jlim_lo    =1                                 &  !<-- Lower J limit of child's segment of boundary
                                     ,jlim_hi    =N_BLEND_V                         &  !<-- Upper J limit of child's segment of boundary
                                     ,i_start    =PARENT_TASK(N)%SOUTH_V%INDX_START &  !<-- Child's segment Istart on each parent task
                                     ,i_end      =PARENT_TASK(N)%SOUTH_V%INDX_END   &  !<-- Child's segment Iend on each parent task
                                     ,j_start    =1                                 &  !<-- Child's segment Jstart on each parent task
                                     ,j_end      =N_BLEND_V                         &  !<-- Child's segment Jend on each parent task
                                     ,ub         =UB_S                              &  !<-- Child's 1-D segment of U on this side of bndry
                                     ,vb         =VB_S )                               !<-- Child's 1-D segment of V on this side of bndry
!!!                                  ,len_2d     =LENGTH_BND_SEG_V%SOUTH )             !<-- Length of the 1-D segment in I and J
          cpl1_south_v_undo_tim=cpl1_south_v_undo_tim+(timef()-btim)

!         call date_and_time(values=values)     
!         write(0,324)n,parent_task(n)%south_v%id_source,values(5),values(6),values(7),values(8)
! 324     format(' After South_V Recv data_from_string task #',i1,' id=',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        btim=timef()
        CALL EXPORT_CHILD_BOUNDARY(ub          =UB_S                    &  !<-- Child's 1-D segment of U on this side of bndry
                                  ,vb          =VB_S                    &  !<-- Child's 1-D segment of V on this side of bndry
                                  ,ilim_lo     =INDX_MIN_V%SOUTH        &  !<-- Lower I limit of child's segment of boundary
                                  ,ilim_hi     =INDX_MAX_V%SOUTH        &  !<-- Upper I limit of child's segment of boundary
                                  ,jlim_lo     =1                       &  !<-- Lower J limit of child's segment of boundary
                                  ,jlim_hi     =N_BLEND_V               &  !<-- Upper J limit of child's segment of boundary
                                  ,data_name   ='SOUTH_V_'//TIME_FLAG   &  !<-- Name attached to the combined exported data
                                  ,data_exp    =BOUND_1D_SOUTH_V        &  !<-- Combined boundary segment V data for child task
                                  ,export_state=EXP_STATE )                !<-- The Parent-Child Coupler export state
!
        cpl1_south_v_exp_tim=cpl1_south_v_exp_tim+(timef()-btim)
!
!       call date_and_time(values=values)     
!       write(0,325)values(5),values(6),values(7),values(8)
! 325   format(' After South_V Recv export at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!       write(0,*)' Recvd South_V for ',time_flag
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
      NP_H=NUM_PARENT_TASKS_SENDING_H%NORTH                                !<-- # of parent tasks sending north boundary H data
!     write(0,*)' PARENT_CHILD_CPL_RUN_RECV north np_h=',np_h
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
!-----------------------------------------------------------------------
!
          CALL CHILD_DATA_FROM_STRING(length_data=PARENT_TASK(N)%NORTH_H%LENGTH     &
                                     ,datastring =PARENT_TASK(N)%NORTH_H%STRING     & 
                                     ,ilim_lo    =INDX_MIN_H%NORTH                  &
                                     ,ilim_hi    =INDX_MAX_H%NORTH                  &
                                     ,jlim_lo    =1                                 &
                                     ,jlim_hi    =N_BLEND_H                         &
                                     ,i_start    =PARENT_TASK(N)%NORTH_H%INDX_START &
                                     ,i_end      =PARENT_TASK(N)%NORTH_H%INDX_END   &
                                     ,j_start    =1                                 &
                                     ,j_end      =N_BLEND_H                         &
                                     ,pdb        =PDB_N                             &
                                     ,tb         =TB_N                              &
                                     ,qb         =QB_N                              &
                                     ,cwb        =CWB_N )
!!!                                  ,len_2d     =LENGTH_BND_SEG_H%NORTH )
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        CALL EXPORT_CHILD_BOUNDARY(pdb         =PDB_N                   &  !<-- Child's 1-D segment of PD on this side of bndry
                                  ,tb          =TB_N                    &  !<-- Child's 1-D segment of T on this side of bndry
                                  ,qb          =QB_N                    &  !<-- Child's 1-D segment of Q on this side of bndry
                                  ,cwb         =CWB_N                   &  !<-- Child's 1-D segment of CW on this side of bndry
                                  ,ilim_lo     =INDX_MIN_H%NORTH        &  !<-- Lower I limit of child's segment of boundary
                                  ,ilim_hi     =INDX_MAX_H%NORTH        &  !<-- Upper I limit of child's segment of boundary
                                  ,jlim_lo     =1                       &  !<-- Lower J limit of child's segment of boundary
                                  ,jlim_hi     =N_BLEND_H               &  !<-- Upper J limit of child's segment of boundary
                                  ,data_name   ='NORTH_H_'//TIME_FLAG   &  !<-- Name attached to the combined exported data
                                  ,data_exp    =BOUND_1D_NORTH_H        &  !<-- Combined boundary segment H data for child task
                                  ,export_state=EXP_STATE )                !<-- The Parent-Child Coupler export state
!
!       write(0,*)' Recvd North_H for ',time_flag
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
      NP_V=NUM_PARENT_TASKS_SENDING_V%NORTH                                !<-- # of parent tasks sending north boundary V data
!     write(0,*)' PARENT_CHILD_CPL_RUN_RECV north np_v=',np_v
!
      IF(NP_V>0)THEN
        NTAG=NTIMESTEP+104                                                 !<-- Add 104 to obtain a unique north H tag
!
        DO N=1,NP_V                                                        !<-- Loop over each parent task sending Nboundary V data
!     write(0,*)' RECV_NEST_BC_DATA for V north ready to recv ',PARENT_TASK(N)%NORTH_V%LENGTH,' words from parent task #',n &
!              ,' with ID=',PARENT_TASK(N)%NORTH_V%ID_SOURCE
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
!     write(0,*)' RECV_NEST_BC_DATA for V north recvd ',PARENT_TASK(N)%NORTH_V%LENGTH,' words from parent task #',n &
!              ,' with ID=',PARENT_TASK(N)%NORTH_V%ID_SOURCE
!
!-----------------------------------------------------------------------
!
          CALL CHILD_DATA_FROM_STRING(length_data=PARENT_TASK(N)%NORTH_V%LENGTH     &
                                     ,datastring =PARENT_TASK(N)%NORTH_V%STRING     & 
                                     ,ilim_lo    =INDX_MIN_V%NORTH                  &
                                     ,ilim_hi    =INDX_MAX_V%NORTH                  &
                                     ,jlim_lo    =1                                 &
                                     ,jlim_hi    =N_BLEND_V                         &
                                     ,i_start    =PARENT_TASK(N)%NORTH_V%INDX_START &
                                     ,i_end      =PARENT_TASK(N)%NORTH_V%INDX_END   &
                                     ,j_start    =1                                 &
                                     ,j_end      =N_BLEND_V                         &
                                     ,ub         =UB_N                              &
                                     ,vb         =VB_N )
!!!                                  ,len_2d     =LENGTH_BND_SEG_V%NORTH )
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        CALL EXPORT_CHILD_BOUNDARY(ub          =UB_N                    &  !<-- Child's 1-D segment of U on this side of bndry
                                  ,vb          =VB_N                    &  !<-- Child's 1-D segment of V on this side of bndry
                                  ,ilim_lo     =INDX_MIN_V%NORTH        &  !<-- Lower I limit of child's segment of boundary
                                  ,ilim_hi     =INDX_MAX_V%NORTH        &  !<-- Upper I limit of child's segment of boundary
                                  ,jlim_lo     =1                       &  !<-- Lower J limit of child's segment of boundary
                                  ,jlim_hi     =N_BLEND_V               &  !<-- Upper J limit of child's segment of boundary
                                  ,data_name   ='NORTH_V_'//TIME_FLAG   &  !<-- Name attached to the combined exported data
                                  ,data_exp    =BOUND_1D_NORTH_V        &  !<-- Combined boundary segment V data for child task
                                  ,export_state=EXP_STATE )                !<-- The Parent-Child Coupler export state
!
      ENDIF
!       write(0,*)' Recvd North_V for ',time_flag
!
      cpl1_north_v_tim=cpl1_north_v_tim+(timef()-btim0)
!
!-------------------
!***  West H Points
!-------------------
!
      btim0=timef()
!
      NP_H=NUM_PARENT_TASKS_SENDING_H%WEST                                 !<-- # of parent tasks sending west boundary H data
!     write(0,*)' PARENT_CHILD_CPL_RUN_RECV west np_h=',np_h
!
      IF(NP_H>0)THEN
        NTAG=NTIMESTEP+105                                                 !<-- Add 105 to obtain a unique west H tag
!
        DO N=1,NP_H                                                        !<-- Loop over each parent task sending Wboundary H data
!     write(0,*)' RECV_NEST_BC_DATA for H west ready to recv ',PARENT_TASK(N)%WEST_H%LENGTH,' words from parent task #',n &
!              ,' with ID=',PARENT_TASK(N)%WEST_H%ID_SOURCE,' ntag=',ntag
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
!-----------------------------------------------------------------------
!
          CALL CHILD_DATA_FROM_STRING(length_data=PARENT_TASK(N)%WEST_H%LENGTH      &
                                     ,datastring =PARENT_TASK(N)%WEST_H%STRING      & 
                                     ,ilim_lo    =1                                 &
                                     ,ilim_hi    =N_BLEND_H                         &
                                     ,jlim_lo    =INDX_MIN_H%WEST                   &
                                     ,jlim_hi    =INDX_MAX_H%WEST                   &
                                     ,i_start    =1                                 &
                                     ,i_end      =N_BLEND_H                         &
                                     ,j_start    =PARENT_TASK(N)%WEST_H%INDX_START  &
                                     ,j_end      =PARENT_TASK(N)%WEST_H%INDX_END    &
                                     ,pdb        =PDB_W                             &
                                     ,tb         =TB_W                              &
                                     ,qb         =QB_W                              &
                                     ,cwb        =CWB_W )
!!!                                  ,len_2d     =LENGTH_BND_SEG_H%WEST )
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        CALL EXPORT_CHILD_BOUNDARY(pdb         =PDB_W                   &  !<-- Child's 1-D segment of PD on this side of bndry
                                  ,tb          =TB_W                    &  !<-- Child's 1-D segment of T on this side of bndry
                                  ,qb          =QB_W                    &  !<-- Child's 1-D segment of Q on this side of bndry
                                  ,cwb         =CWB_W                   &  !<-- Child's 1-D segment of CW on this side of bndry
                                  ,ilim_lo     =1                       &  !<-- Lower I limit of child's segment of boundary
                                  ,ilim_hi     =N_BLEND_H               &  !<-- Upper I limit of child's segment of boundary
                                  ,jlim_lo     =INDX_MIN_H%WEST         &  !<-- Lower J limit of child's segment of boundary
                                  ,jlim_hi     =INDX_MAX_H%WEST         &  !<-- Upper J limit of child's segment of boundary
                                  ,data_name   ='WEST_H_'//TIME_FLAG    &  !<-- Name attached to the combined exported data
                                  ,data_exp    =BOUND_1D_WEST_H         &  !<-- Combined boundary segment H data for child task
                                  ,export_state=EXP_STATE )                !<-- The Parent-Child Coupler export state
!
!       write(0,*)' Recvd West_H for ',time_flag
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
      NP_V=NUM_PARENT_TASKS_SENDING_V%WEST                                 !<-- # of parent tasks sending west boundary V data
!     write(0,*)' PARENT_CHILD_CPL_RUN_RECV west np_v=',np_v
!
      IF(NP_V>0)THEN
        NTAG=NTIMESTEP+106                                                 !<-- Add 106 to obtain a unique west H tag
!
        DO N=1,NP_V                                                        !<-- Loop over each parent task sending Sboundary V data
!     write(0,*)' RECV_NEST_BC_DATA for V west ready to recv ',PARENT_TASK(N)%WEST_V%LENGTH,' words from parent task #',n &
!              ,' with ID=',PARENT_TASK(N)%WEST_V%ID_SOURCE
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
!     write(0,*)' RECV_NEST_BC_DATA for V west recvd ',PARENT_TASK(N)%WEST_V%LENGTH,' words from parent task #',n &
!              ,' with ID=',PARENT_TASK(N)%WEST_V%ID_SOURCE
!
!-----------------------------------------------------------------------
!
          CALL CHILD_DATA_FROM_STRING(length_data=PARENT_TASK(N)%WEST_V%LENGTH      &
                                     ,datastring =PARENT_TASK(N)%WEST_V%STRING      & 
                                     ,ilim_lo    =1                                 &
                                     ,ilim_hi    =N_BLEND_V                         &
                                     ,jlim_lo    =INDX_MIN_V%WEST                   &
                                     ,jlim_hi    =INDX_MAX_V%WEST                   &
                                     ,i_start    =1                                 &
                                     ,i_end      =N_BLEND_V                         &
                                     ,j_start    =PARENT_TASK(N)%WEST_V%INDX_START  &
                                     ,j_end      =PARENT_TASK(N)%WEST_V%INDX_END    &
                                     ,ub         =UB_W                              &
                                     ,vb         =VB_W )
!!!                                  ,len_2d     =LENGTH_BND_SEG_V%WEST )
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        CALL EXPORT_CHILD_BOUNDARY(ub          =UB_W                    &  !<-- Child's 1-D segment of U on this side of bndry
                                  ,vb          =VB_W                    &  !<-- Child's 1-D segment of V on this side of bndry
                                  ,ilim_lo     =1                       &  !<-- Lower I limit of child's segment of boundary
                                  ,ilim_hi     =N_BLEND_V               &  !<-- Upper I limit of child's segment of boundary
                                  ,jlim_lo     =INDX_MIN_V%WEST         &  !<-- Lower J limit of child's segment of boundary
                                  ,jlim_hi     =INDX_MAX_V%WEST         &  !<-- Upper J limit of child's segment of boundary
                                  ,data_name   ='WEST_V_'//TIME_FLAG    &  !<-- Name attached to the combined exported data
                                  ,data_exp    =BOUND_1D_WEST_V         &  !<-- Combined boundary segment V data for child task
                                  ,export_state=EXP_STATE )                !<-- The Parent-Child Coupler export state
!
      ENDIF
!       write(0,*)' Recvd West_V for ',time_flag
!
      cpl1_west_v_tim=cpl1_west_v_tim+(timef()-btim0)
!
!-------------------
!***  East H Points
!-------------------
!
      btim0=timef()
!
      NP_H=NUM_PARENT_TASKS_SENDING_H%EAST                                 !<-- # of parent tasks sending east boundary H data
!     write(0,*)' PARENT_CHILD_CPL_RUN_RECV east np_h=',np_h
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
!-----------------------------------------------------------------------
!
          CALL CHILD_DATA_FROM_STRING(length_data=PARENT_TASK(N)%EAST_H%LENGTH      &
                                     ,datastring =PARENT_TASK(N)%EAST_H%STRING      & 
                                     ,ilim_lo    =1                                 &
                                     ,ilim_hi    =N_BLEND_H                         &
                                     ,jlim_lo    =INDX_MIN_H%EAST                   &
                                     ,jlim_hi    =INDX_MAX_H%EAST                   &
                                     ,i_start    =1                                 &
                                     ,i_end      =N_BLEND_H                         &
                                     ,j_start    =PARENT_TASK(N)%EAST_H%INDX_START  &
                                     ,j_end      =PARENT_TASK(N)%EAST_H%INDX_END    &
                                     ,pdb        =PDB_E                             &
                                     ,tb         =TB_E                              &
                                     ,qb         =QB_E                              &
                                     ,cwb        =CWB_E )
!!!                                  ,len_2d     =LENGTH_BND_SEG_H%EAST )
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        CALL EXPORT_CHILD_BOUNDARY(pdb         =PDB_E                   &  !<-- Child's 1-D segment of PD on this side of bndry
                                  ,tb          =TB_E                    &  !<-- Child's 1-D segment of T on this side of bndry
                                  ,qb          =QB_E                    &  !<-- Child's 1-D segment of Q on this side of bndry
                                  ,cwb         =CWB_E                   &  !<-- Child's 1-D segment of CW on this side of bndry
                                  ,ilim_lo     =1                       &  !<-- Lower I limit of child's segment of boundary
                                  ,ilim_hi     =N_BLEND_H               &  !<-- Upper I limit of child's segment of boundary
                                  ,jlim_lo     =INDX_MIN_H%EAST         &  !<-- Lower J limit of child's segment of boundary
                                  ,jlim_hi     =INDX_MAX_H%EAST         &  !<-- Upper J limit of child's segment of boundary
                                  ,data_name   ='EAST_H_'//TIME_FLAG    &  !<-- Name attached to the combined exported data
                                  ,data_exp    =BOUND_1D_EAST_H         &  !<-- Combined boundary segment H data for child task
                                  ,export_state=EXP_STATE )                !<-- The Parent-Child Coupler export state
!
!       write(0,*)' Recvd East_H for ',time_flag
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
      NP_V=NUM_PARENT_TASKS_SENDING_V%EAST                                 !<-- # of parent tasks sending east boundary V data
!     write(0,*)' PARENT_CHILD_CPL_RUN_RECV east np_v=',np_v
!
      IF(NP_V>0)THEN
        NTAG=NTIMESTEP+108                                                 !<-- Add 108 to obtain a unique east H tag
!
        DO N=1,NP_V                                                        !<-- Loop over each parent task sending Eboundary V data
!     write(0,*)' RECV_NEST_BC_DATA for V east ready to recv ',PARENT_TASK(N)%EAST_V%LENGTH,' words from parent task #',n &
!              ,' with ID=',PARENT_TASK(N)%EAST_V%ID_SOURCE
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
!     write(0,*)' RECV_NEST_BC_DATA for V east recvd ',PARENT_TASK(N)%EAST_V%LENGTH,' words from parent task #',n &
!              ,' with ID=',PARENT_TASK(N)%EAST_V%ID_SOURCE
!
!-----------------------------------------------------------------------
!
          CALL CHILD_DATA_FROM_STRING(length_data=PARENT_TASK(N)%EAST_V%LENGTH      &
                                     ,datastring =PARENT_TASK(N)%EAST_V%STRING      & 
                                     ,ilim_lo    =1                                 &
                                     ,ilim_hi    =N_BLEND_V                         &
                                     ,jlim_lo    =INDX_MIN_V%EAST                   &
                                     ,jlim_hi    =INDX_MAX_V%EAST                   &
                                     ,i_start    =1                                 &
                                     ,i_end      =N_BLEND_V                         &
                                     ,j_start    =PARENT_TASK(N)%EAST_V%INDX_START  &
                                     ,j_end      =PARENT_TASK(N)%EAST_V%INDX_END    &
                                     ,ub         =UB_E                              &
                                     ,vb         =VB_E )
!!!                                  ,len_2d     =LENGTH_BND_SEG_V%EAST )
!
        ENDDO
!
!-----------------------------------------------------------------------
!
        CALL EXPORT_CHILD_BOUNDARY(ub          =UB_E                    &  !<-- Child's 1-D segment of U on this side of bndry
                                  ,vb          =VB_E                    &  !<-- Child's 1-D segment of V on this side of bndry
                                  ,ilim_lo     =1                       &  !<-- Lower I limit of child's segment of boundary
                                  ,ilim_hi     =N_BLEND_V               &  !<-- Upper I limit of child's segment of boundary
                                  ,jlim_lo     =INDX_MIN_V%EAST         &  !<-- Lower J limit of child's segment of boundary
                                  ,jlim_hi     =INDX_MAX_V%EAST         &  !<-- Upper J limit of child's segment of boundary
                                  ,data_name   ='EAST_V_'//TIME_FLAG    &  !<-- Name attached to the combined exported data
                                  ,data_exp    =BOUND_1D_EAST_V         &  !<-- Combined boundary segment V data for child task
                                  ,export_state=EXP_STATE )                !<-- The Parent-Child Coupler export state
!
      ENDIF
!       write(0,*)' Recvd East_V for ',time_flag
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
!     write(0,*)' exit RECV_NEST_BC_DATA for ',time_flag
!-----------------------------------------------------------------------
!
      END SUBROUTINE RECV_NEST_BC_DATA
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
!***  Run the coupler step where parents send data to their children.
!***  This is phase 2 of the coupler since it occurs at the end of
!***  the timesteps.  The children receive the data in phase 1 at
!***  the beginning of their appropriate timesteps.  Only parent
!***  tasks enter this routine.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_CplComp) :: CPL_COMP                                       !<-- The Parent-Child Coupler Component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Coupler's Import State
                         ,EXP_STATE                                        !<-- The Coupler's Export State
!
      TYPE(ESMF_Clock) :: CLOCK                                            !<-- The NMM Clock for this parent domain
!
      INTEGER,INTENT(OUT) :: RC_FINAL
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I_PARENT_SW_OLD,J_PARENT_SW_OLD             &
                           ,KOUNT_MOVING,N,N_MOVING                     &
                           ,N_UPDATE_CHILD_TASKS                        &
!!!                        ,N_UPDATE_CHILD_TASKS_V                      &
                           ,NR,NRES                                     &
                           ,NTIMESTEP                                   &
                           ,NTIMESTEP_CHILD                             &
                           ,NTIMESTEP_MOVE                              &
                           ,NUM_CHILD_TASKS,SPACE_RATIO
!
      INTEGER(kind=KINT) :: IERR,IRTN,RC,RC_CPL_RUN
!
      INTEGER(kind=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE,SAVE ::               &
                                                  PROCEED_AFTER_BC_RECV
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER :: CHILD_TASK_LIMITS
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      LOGICAL(kind=KLOG) :: EXCH_DONE,MOVE_FLAG_IS_PRESENT              &
                           ,PARENT_MOVED
!
      integer,dimension(8) :: values
      integer(kind=kint) :: kount_h_links
      type(child_update_link),pointer :: xxx
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
!***  Initialize the error signal variables.
!-----------------------------------------------------------------------
!
      RC        =ESMF_SUCCESS
      RC_FINAL  =ESMF_SUCCESS
      RC_CPL_RUN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      KOUNT_MOVING=0                                                       !<-- Keep track of nests who want to move.
!
      CALL ESMF_CLOCKGET(CLOCK       =CLOCK                             &  !<-- The ESMF Clock
                        ,advanceCount=ntimestep_esmf                    &  !<-- The parent's current timestep (ESMF)
                        ,rc          =rc)
!
      NTIMESTEP=NTIMESTEP_ESMF                                             !<-- The parent just finished the current timestep
!
!-----------------------------------------------------------------------
!***  The following block is for the very special setup in which
!***  the parent domain is a moving nest and it contains children
!***  (which must be moving and cannot be static).
!***  If the parent just moved at the beginning of the current timestep
!***  it must adjust its location of its children.
!***  Note that I_SHIFT_CHILD and J_SHIFT_CHILD here are the shift
!***  values of this parent domain on its own grid inherited from
!***  subroutine PARENT_CHILD_CPL_RUN_RECV where it was a child.
!-----------------------------------------------------------------------
!
      PARENT_MOVED=.FALSE.
!
#ifdef ESMF_3
      parent_moves: IF(MY_DOMAIN_MOVES==ESMF_TRUE)THEN
#else
      parent_moves: IF(MY_DOMAIN_MOVES)THEN                                !<-- Does this parent domain move?
#endif
!
        this_timestep: IF(NTIMESTEP==NEXT_MOVE_TIMESTEP)THEN
!
          PARENT_MOVED=.TRUE.                                              !<-- Parent moved at beginning of current timestep
!
          DO N=1,NUM_CHILDREN
            I_PARENT_SW(N)=I_PARENT_SW(N)-I_SHIFT_CHILD                    !<-- Child N's new SW corner I after parent moved
            J_PARENT_SW(N)=J_PARENT_SW(N)-J_SHIFT_CHILD                    !<-- Child N's new SW corner J after parent moved
          ENDDO
!
!-----------------------------------------------------------------------
!***  Now fill the parent's data objects that hold the nest-resolution
!***  topography at child H points and child V points.  Since the
!***  parent just moved then its MPI task subdomains need the data
!***  for their new locations.  Note that in this situation the total
!***  number of children must equal the number of moving children.
!-----------------------------------------------------------------------
!
          CALL PARENT_READS_MOVING_CHILD_TOPO(NUM_CHILDREN              &
                                             ,LINK_MRANK_RATIO          &
                                             ,LIST_OF_RATIOS            &
                                             ,M_NEST_RATIO              &
                                             ,KOUNT_RATIOS_MN           &
                                             ,IM_1,JM_1                 &
                                             ,TPH0_1,TLM0_1             &
                                             ,SB_1,WB_1                 &
                                             ,RECIP_DPH_1,RECIP_DLM_1   &
                                             ,GLAT,GLON                 &
                                             ,NEST_FIS_ON_PARENT_BNDS   &
                                             ,NEST_FIS_ON_PARENT        &
                                             ,NEST_FIS_V_ON_PARENT      &
                                             ,IDS,IDE,IMS,IME,ITS,ITE   &
                                             ,JDS,JDE,JMS,JME,JTS,JTE)
!
        ENDIF  this_timestep
!
      ENDIF  parent_moves
!
!-----------------------------------------------------------------------
!***  The parent generates new boundary data for all of its children
!***  given their domains' positions at the beginning of this parent
!***  timestep and sends it to the children so they can form time
!***  tendencies for their boundary variables as they integrate through
!***  this parent timestep.  This is relevant for all children, both
!***  static and moving.  However, if the parent learned a timestep
!***  ago that a child wants to move then it must now reset those
!***  working pointers/arrays that are used for the preparation of
!***  the standard child boundary updates that are sent back in time
!***  to all children every parent timestep so the child can generate
!***  its boundary tendencies.  The reset is needed because the child's
!***  boundary has different associations with the parent tasks after
!***  the move.  The same work is needed if the parent domain moved
!***  at the beginning of this timestep since that also changes the
!***  association of parent tasks and child boundary tasks.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_CHILDREN
!
        IF(STATIC_OR_MOVING(N)=='Moving')THEN                              !<-- Select the children who can move.
          KOUNT_MOVING=KOUNT_MOVING+1
!
          IF(MOVE_FLAG(KOUNT_MOVING).OR.PARENT_MOVED)THEN                  !<-- If true, this child just moved relative to parent
            NRES=LINK_MRANK_RATIO(KOUNT_MOVING)                            !<-- Rank of space ratio value among the moving children
            CALL RESET_WORK_PARENT(N,NRES,'Future',PARENT_MOVED)           !<-- Reset working arrays for this moving nest.
            MOVE_FLAG(KOUNT_MOVING)=.FALSE.                                !<-- Flag will be 'false' until the next move. 
          ENDIF
!
        ENDIF
!
        CALL COMPUTE_SEND_NEST_BC_DATA(N,'Future')                         !<-- Parent sends BC data to all children from their future.
!
      ENDDO
!
      IF(PARENT_MOVED)PARENT_MOVED=.FALSE.
!
!-----------------------------------------------------------------------
!***  We are at the end of a parent timestep.  If the parent has
!***  children who move then:
!***
!***  (1) The parent receives a message from each child who can move
!***      only when the child wants to move.  A 2nd message contains
!***      the new location of the nest domain's SW corner on the 
!***      parent grid.  The signal is received by the parent at the
!***      end of a parent timestep while the child will have sent
!***      these messages from the beginning of this parent timestep
!***      or from the beginning of the previous one depending on
!***      the relative speeds of parent and child.
!***  (2) When the child sends its move message it can't know the
!***      parent's current timestep which is the timestep at which
!***      the parent will generate new internal nest data as well as
!***      new boundary data to begin forming time tendencies at the
!***      child's new location.
!***      Let's assume a moving child reaches the beginning of parent
!***      timestep N and decides it wants to move.  It sends a move
!***      signal to the parent which will receive it at either the
!***      end of timestep N or N+1.  The parent then sends back to
!***      the child the value of the parent's current timestep so
!***      the child will know where in time to receive from its
!***      parent the new internal data and starting boundary data
!***      for its new location.  That timestep then must be the one 
!***      at which the child will actually be able to execute its move.
!***  (3) The parent computes and sends new information to the moving
!***      children regarding the association of parent and child tasks
!***      for the children's new locations after they move.  Then the
!***      parent computes and sends the new internal child data for
!***      those child gridpoints that have moved over a new region of
!***      the parent grid as well as the new starting boundary data
!***      for their grids' new locations.
!-----------------------------------------------------------------------
!
      moving_children: IF(NUM_MOVING_CHILDREN>0)THEN                       !<-- Select all of this parent's moving children
!
!-----------------------------------------------------------------------
!
        CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)                               !<-- Syncs Probe below with BC ISends above; required
!
        EXCH_DONE=.FALSE.                                                  !<-- Initialize flag for parent halo exchanges
!
!-----------------------------------------------------------------------
!***  We must keep some of a parent's moving children from greatly
!***  outpacing the others or else the parent-child exchanges of
!***  information become scrambled.  Therefore force the moving children
!***  to wait until all their moving siblings have also received BC
!***  updates before any of them can proceed.  Do this with another
!***  flag.
!-----------------------------------------------------------------------
!
        parent_task_0: IF(MYPE==0)THEN         
!
!-----------------------------------------------------------------------
!
          IF(.NOT.ALLOCATED(PROCEED_AFTER_BC_RECV))THEN
            ALLOCATE(PROCEED_AFTER_BC_RECV(1:NUM_MOVING_CHILDREN))
            DO N=1,NUM_MOVING_CHILDREN
              PROCEED_AFTER_BC_RECV(N)=1
            ENDDO
          ENDIF
!
          child_loop_0: DO N=1,NUM_MOVING_CHILDREN                         !<-- Loop through this parent's moving children
!
            N_MOVING=RANK_MOVING_CHILD(N)                                  !<-- In the list of this parent's children, these are moving.
!
            CALL MPI_WAIT(HANDLE_BC_UPDATE(N)                           &  !<-- Handle for ISend of moving child N's BC update flag
                         ,JSTAT                                         &  !<-- MPI status
                         ,IERR)
!
            CALL MPI_ISEND(PROCEED_AFTER_BC_RECV(N)                     &  !<-- Parent task 0 sends BC update flag to moving child N
                          ,1                                            &  !<-- # of words sent
                          ,MPI_INTEGER                                  &  !<-- The flag is an integer.
                          ,0                                            &  !<-- Sending to moving child N's task 0.
                          ,MOVING_BC_TAG                                &  !<-- Tag associated with moving nests' recvs of BC data
                          ,COMM_TO_MY_CHILDREN(N_MOVING)                &  !<-- MPI communicator between parent and moving child N
                          ,HANDLE_BC_UPDATE(N)                          &  !<-- Handle for this ISend of moving child N's BC update flag
                          ,IERR)
!
          ENDDO child_loop_0
!
!-----------------------------------------------------------------------
!
          child_loop_1: DO N=1,NUM_MOVING_CHILDREN                         !<-- Loop through this parent's moving children
!
!-----------------------------------------------------------------------
!
            N_MOVING=RANK_MOVING_CHILD(N)                                  !<-- In the list of this parent's children, these are moving.
!
!-----------------------------------------------------------------------
!
      btim=timef()
!
            MOVE_FLAG(N)=.FALSE.
!
            CALL MPI_IPROBE(0                                           &  !<-- Is move message present from moving child N's fcst task 0?
                           ,MOVE_TAG                                    &  !<-- Tag associated with nest N's move flag
                           ,COMM_TO_MY_CHILDREN(N_MOVING)               &  !<-- MPI communicator between parent and moving child N
                           ,MOVE_FLAG_IS_PRESENT                        &  !<-- Is the nest's move flag now available?
                           ,JSTAT                                       &
                           ,IERR)
!
!-----------------------------------------------------------------------
!
            IF(MOVE_FLAG_IS_PRESENT)THEN
!
              MOVE_FLAG(N)=.TRUE.                                          !<-- Moving child N is saying it wants to move
!
              CALL MPI_RECV(MOVE_FLAG(N)                                &  !<-- Recv the message to clear the nest's ISEND 
                           ,1                                           &  !<-- # of words in message
                           ,MPI_LOGICAL                                 &  !<-- The message is type Logical.
                           ,0                                           &  !<-- The message was sent by moving child N's fcst task 0.
                           ,MOVE_TAG                                    &  !<-- Arbitrary tag used for this data exchange
                           ,COMM_TO_MY_CHILDREN(N_MOVING)               &  !<-- MPI communicator between parent and moving child N
                           ,JSTAT                                       &
                           ,IERR)
!
              CALL MPI_WAIT(HANDLE_TIMESTEP(N)                          &  !<-- Handle for ISend of parent's timestep for nest move
                           ,JSTAT                                       &  !<-- MPI status
                           ,IERR)
!
              NTIMESTEP_MOVE=NTIMESTEP                                     !<-- Safe to reload the buffer variable
!
              CALL MPI_ISEND(NTIMESTEP_MOVE                             &  !<-- Parent task 0 sends its current timestep.
                            ,1                                          &  !<-- # of words sent
                            ,MPI_INTEGER                                &  !<-- The timestep is an integer.
                            ,0                                          &  !<-- Sending to moving child N's task 0.
                            ,TIME_TAG                                   &  !<-- Tag associated with nest N's move timestep
                            ,COMM_TO_MY_CHILDREN(N_MOVING)              &  !<-- MPI communicator between parent and moving child N
                            ,HANDLE_TIMESTEP(N)                         &  !<-- Handle for this ISend to moving child N's task 0
                            ,IERR)
!
            ENDIF
!
      t0_recv_move_tim=t0_recv_move_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!
          ENDDO child_loop_1
!
!-----------------------------------------------------------------------
!
        ENDIF parent_task_0
!
!-----------------------------------------------------------------------
!
        child_loop_2: DO N=1,NUM_MOVING_CHILDREN                           !<-- Loop through this parent's moving children
!
!-----------------------------------------------------------------------
!***  Parent task 0 informs the other parent tasks if moving child N 
!***  has signaled that it wants to move and if so what is the new 
!***  location of that child's southwest corner.
!-----------------------------------------------------------------------
!
          CALL MPI_BCAST(MOVE_FLAG(N)                                   &  !<-- Moving child N's signal:  Does it want to move?
                        ,1                                              &  !<-- The timestep is one word
                        ,MPI_LOGICAL                                    &  !<-- The signal is type Logical
                        ,0                                              &  !<-- Broadcast from parent forecast task 0
                        ,MPI_COMM_COMP                                  &  !<-- MPI communicator for this parent's forecast tasks
                        ,IRTN )
!
          N_MOVING=RANK_MOVING_CHILD(N)                                    !<-- In the list of this parent's children, these are moving.
          NTIMESTEP_CHILD=(NTIMESTEP+1)                                 &  !<-- The nest's timestep at which it will recv parent data.
                          *PARENT_CHILD_SPACE_RATIO(N_MOVING)
!
!-----------------------------------------------------------------------
!***  Now all parent tasks consider the children who want to move.
!-----------------------------------------------------------------------
!
          child_moves: IF(MOVE_FLAG(N))THEN                                !<-- If true, child N wants to move.
!
            I_PARENT_SW_OLD=I_PARENT_SW(N_MOVING)                          !<-- Save the previous location of the nest.
            J_PARENT_SW_OLD=J_PARENT_SW(N_MOVING)                          !<--
!
!-----------------------------------------------------------------------
!***  Parent task 0 needs to broadcast the moving child's new SW corner
!***  location to the other parent tasks.  Then the parent tasks who
!***  will contain child boundary tasks after the child moves must
!***  fill their working arrays with the child boundary topography.
!-----------------------------------------------------------------------
!
            IF(MYPE==0)THEN
              CALL MPI_RECV(IJ_SHIFT_CHILD                              &  !<-- Moving child N's change in SW corner I,J on parent grid
                           ,2                                           &  !<-- # of words in message
                           ,MPI_INTEGER                                 &  !<-- The message is type Integer
                           ,0                                           &  !<-- The message was sent by moving child N's fcst task 0.
                           ,9001                                        &  !<-- Arbitrary tag used for this data exchange
                           ,COMM_TO_MY_CHILDREN(N_MOVING)               &  !<-- MPI communicator between parent and moving child N
                           ,JSTAT                                       &
                           ,IERR)
            ENDIF
!
            CALL MPI_BCAST(IJ_SHIFT_CHILD                               &  !<-- Moving child N's change in SW corner
                          ,2                                            &  !<-- Two indices
                          ,MPI_INTEGER                                  &  !<-- The signal is type Integer
                          ,0                                            &  !<-- Broadcast from parent forecast task 0
                          ,MPI_COMM_COMP                                &  !<-- MPI communicator for this parent's forecast tasks
                          ,IRTN )
!
            I_PARENT_SW(N_MOVING)=I_PARENT_SW_OLD+IJ_SHIFT_CHILD(1)        !<-- Child N to move its SW corner to this parent I
            J_PARENT_SW(N_MOVING)=J_PARENT_SW_OLD+IJ_SHIFT_CHILD(2)        !<-- Child N to move its SW corner to this parent J
!
!-----------------------------------------------------------------------
!***  If this child wants to move then reset the working arrays/pointers
!***  used to generate values interpolated from the parent to child's
!***  boundary immediately after a move.  This set of working objects
!***  is separate from the standard ones used to interpolate boundary
!***  data for all nests since when a nest moves we need to have the
!***  objects in place for both the old location and the new until we
!***  know for certain those for the old location have been received by
!***  the moving child.
!
!***  Note that N_MOVING is the rank of the moving child among ALL of
!***  children and NRES is the rank of the moving child's space ratio
!***  in the list of different space ratios for all moving children.
!-----------------------------------------------------------------------
!
            NRES=LINK_MRANK_RATIO(N)                                       !<-- Rank of space ratio value among the moving children
!
            CALL RESET_WORK_PARENT(N_MOVING,NRES,'Current',PARENT_MOVED)
!
!-----------------------------------------------------------------------
!***  The parent generates and sends new boundary data for the child's 
!***  new position that it will move to when it reaches this point in
!***  time that the parent is at now.
!-----------------------------------------------------------------------
!
            CALL COMPUTE_SEND_NEST_BC_DATA(N_MOVING,'Current')
!
!-----------------------------------------------------------------------
!***  Parent tasks determine the index limits of the regions on
!***  moving nest N's task subdomains that they are responsible
!***  for updating after the nest moves.  Those index limits are
!***  identical for H and V points.
!-----------------------------------------------------------------------
!
            CHILD_TASK_LIMITS=>CTASK_LIMITS(N_MOVING)%LIMITS
            NUM_CHILD_TASKS=FTASKS_DOMAIN(MY_CHILDREN_ID(N_MOVING))
            SPACE_RATIO=PARENT_CHILD_SPACE_RATIO(N_MOVING)
!
            N_UPDATE_CHILD_TASKS=0
!
            btim=timef()
            CALL PARENT_BOOKKEEPING_MOVING(I_PARENT_SW(N_MOVING)        &  !<-- SW corner of nest is on this parent I after move
                                          ,J_PARENT_SW(N_MOVING)        &  !<-- SW corner of nest is on this parent J after move
                                          ,I_PARENT_SW_OLD              &  !<-- SW corner of nest is on this parent I before move
                                          ,J_PARENT_SW_OLD              &  !<-- SW corner of nest is on this parent J before move
                                          ,ITS,ITE,JTS,JTE              &  !<-- ITS,ITE,JTS,JTE for this parent task
                                          ,NUM_CHILD_TASKS              &  !<-- # of child forecast tasks
                                          ,CHILD_TASK_LIMITS            &  !<-- ITS,ITE,JTS,JTE for each child forecast task
                                          ,SPACE_RATIO                  &  !<-- # of child grid increments in one of parent's
                                          ,NHALO                        &  !<-- # of halo points
                                          ,NROWS_P_UPD_W                &  !<-- Moving nest footprint W bndry rows updated by parent
                                          ,NROWS_P_UPD_E                &  !<-- Moving nest footprint E bndry rows updated by parent
                                          ,NROWS_P_UPD_S                &  !<-- Moving nest footprint S bndry rows updated by parent
                                          ,NROWS_P_UPD_N                &  !<-- Moving nest footprint N bndry rows updated by parent
                                          ,N_UPDATE_CHILD_TASKS         &  !<-- # of moving nest tasks updated by this parent task
!!!                                       ,TASK_UPDATE_SPECS_H(N)       &  !<-- Linked list of nest task update region specs for H points
                                          ,TASK_UPDATE_SPECS(N)         &  !<-- Linked list of nest task update region specs
                                          ,HANDLE_MOVE_DATA(N)%NTASKS_TO_RECV &  !<-- MPI handles for update data ISent to nest N's tasks
                                          ,MOVING_CHILD_UPDATE(N)       &  !<-- Composite H/V update data for tasks on moving child N
                                            )
!
!-----------------------------------------------------------------------
!***  When a parent needs to update data on some of its moving nests
!***  and any of those nest points lie between ITE/JTE on one parent
!***  task and ITS/JTS on an adjacent parent task then values from 
!***  the parent tasks' halo regions must be used.  However some of
!***  the variables that need updating are not computed in the halo
!***  regions.  That means that prior to proceeding with moving nest
!***  updates the parent needs to do special halo exchanges for all
!***  those variables required for moving nest updates but for which
!***  halo exchanges were not performed during the normal integration.
!***  Of course these parent tasks' halo exchanges need to be done
!***  only once in a timestep in which any number of its nests move. 
!-----------------------------------------------------------------------
!
            IF(.NOT.EXCH_DONE)THEN
!
              CALL PARENT_UPDATES_HALOS('H'                               &
                                        ,MOVE_BUNDLE_H                    &
                                        ,NUM_FIELDS_MOVE_3D_H             &
                                        ,NUM_FIELDS_MOVE_2D_H_R           &
                                        ,nflds_2di=NUM_FIELDS_MOVE_2D_H_I)
!
              CALL PARENT_UPDATES_HALOS('V'                             &
                                        ,MOVE_BUNDLE_V                  &
                                        ,NUM_FIELDS_MOVE_3D_V           &
                                        ,NUM_FIELDS_MOVE_2D_V )
!
              EXCH_DONE=.TRUE.
!
            ENDIF
!
     parent_bookkeep_moving_tim=parent_bookkeep_moving_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  While the index limits of each parent update region of each 
!***  moving nest are identical for H and V points the routine that
!***  performs the updating will be called separately for H and V
!***  points.  That is because of the different physical locations
!***  of H versus V points which must be accounted for when finding 
!***  the parent's four surrounding points for bilinear interpolations.
!-----------------------------------------------------------------------
!
            btim=timef()
            IF(N_UPDATE_CHILD_TASKS>0)THEN
!
!-----------------------------------------------------------------------
!***  First do the H point updates for moving nest N.
!-----------------------------------------------------------------------
!
              NR=M_NEST_RATIO(N)
!
              CALL PARENT_UPDATES_MOVING('H'                                 &
                                        ,N_UPDATE_CHILD_TASKS                &
                                        ,SPACE_RATIO                         &
                                        ,TIME_RATIO_MY_CHILDREN(N_MOVING)    &
                                        ,NTIMESTEP_CHILD                     &
                                        ,I_PARENT_SW(N_MOVING)               &
                                        ,J_PARENT_SW(N_MOVING)               &
                                        ,PT,PDTOP,PSGML1,SGML2,SG1,SG2       &
                                        ,DSG2,PDSG1                          &
                                        ,FIS,PD                              &
                                        ,T                                   &
                                        ,TRACERS(:,:,:,INDX_Q)               &
                                        ,TRACERS(:,:,:,INDX_CW)              &
                                        ,NUM_CHILD_TASKS                     &  !<-- # of child forecast tasks
                                        ,CHILD_TASK_LIMITS                   &  !<-- ITS,ITE,JTS,JTE for each child forecast task
                                        ,IMS,IME,JMS,JME                     &  !<-- Subdomain memory limits for parent tasks
                                        ,IDS,IDE,JDS,JDE                     &  !<-- Full parent domain limits
                                        ,LM                                  &
                                        ,NEST_FIS_ON_PARENT_BNDS(NRES)%LBND1 &
                                        ,NEST_FIS_ON_PARENT_BNDS(NRES)%UBND1 &
                                        ,NEST_FIS_ON_PARENT_BNDS(NRES)%LBND2 &
                                        ,NEST_FIS_ON_PARENT_BNDS(NRES)%UBND2 &
                                        ,NEST_FIS_ON_PARENT(NRES)%DATA       &
                                        ,COMM_TO_MY_CHILDREN(N)              &
                                        ,HANDLE_MOVE_DATA(N)%NTASKS_TO_RECV  &
                                        ,MOVE_BUNDLE_H                       &
                                        ,NUM_FIELDS_MOVE_2D_H_I              &
                                        ,NUM_FIELDS_MOVE_2D_X_I              &
                                        ,NUM_FIELDS_MOVE_2D_H_R              &
                                        ,NUM_FIELDS_MOVE_2D_X_R              &
                                        ,NUM_FIELDS_MOVE_3D_H                &
                                        ,NUM_LEVELS_MOVE_3D_H                &
                                        ,NUM_FIELDS_MOVE_2D_V                &
                                        ,NUM_FIELDS_MOVE_3D_V                &
                                        ,NUM_LEVELS_MOVE_3D_V                &
!!!                                     ,TASK_UPDATE_SPECS_H(N)              &  !<-- Linked list of nest task update region specs for H 
!!!                                     ,TASK_UPDATE_SPECS_V(N)              &  !<-- Linked list of nest task update region specs for V 
                                        ,TASK_UPDATE_SPECS(N)                &  !<-- Linked list of nest task update region specs 
                                        ,MOVING_CHILD_UPDATE(N)              &  !<-- Composite H/V update data for nest task N
                                          )
!
!-----------------------------------------------------------------------
!***  Now the parent does the V point updates for moving nest N
!***  and then sends all H and V update data to that nest.
!-----------------------------------------------------------------------
!
              NR=M_NEST_RATIO(N)
!
              CALL PARENT_UPDATES_MOVING('V'                                 &
                                        ,N_UPDATE_CHILD_TASKS                &
                                        ,SPACE_RATIO                         &
                                        ,TIME_RATIO_MY_CHILDREN(N_MOVING)    &
                                        ,NTIMESTEP_CHILD                     &
                                        ,I_PARENT_SW(N_MOVING)               &
                                        ,J_PARENT_SW(N_MOVING)               &
                                        ,PT,PDTOP,PSGML1,SGML2,SG1,SG2       &
                                        ,DSG2,PDSG1                          &
                                        ,FIS,PD                              &
                                        ,T                                   &
                                        ,TRACERS(:,:,:,INDX_Q)               &
                                        ,TRACERS(:,:,:,INDX_CW)              &
                                        ,NUM_CHILD_TASKS                     &  !<-- # of child forecast tasks
                                        ,CHILD_TASK_LIMITS                   &  !<-- ITS,ITE,JTS,JTE for each child forecast task
                                        ,IMS,IME,JMS,JME                     &  !<-- Subdomain memory limits for parent tasks
                                        ,IDS,IDE,JDS,JDE                     &  !<-- Full parent domain limits
                                        ,LM                                  &
                                        ,NEST_FIS_ON_PARENT_BNDS(NRES)%LBND1 &
                                        ,NEST_FIS_ON_PARENT_BNDS(NRES)%UBND1 &
                                        ,NEST_FIS_ON_PARENT_BNDS(NRES)%LBND2 &
                                        ,NEST_FIS_ON_PARENT_BNDS(NRES)%UBND2 &
                                        ,NEST_FIS_V_ON_PARENT(NRES)%DATA     &
                                        ,COMM_TO_MY_CHILDREN(N)              &
                                        ,HANDLE_MOVE_DATA(N)%NTASKS_TO_RECV  &
                                        ,MOVE_BUNDLE_V                       &
                                        ,NUM_FIELDS_MOVE_2D_H_I              &
                                        ,NUM_FIELDS_MOVE_2D_X_I              &
                                        ,NUM_FIELDS_MOVE_2D_H_R              &
                                        ,NUM_FIELDS_MOVE_2D_X_R              &
                                        ,NUM_FIELDS_MOVE_3D_H                &
                                        ,NUM_LEVELS_MOVE_3D_H                &
                                        ,NUM_FIELDS_MOVE_2D_V                &
                                        ,NUM_FIELDS_MOVE_3D_V                &
                                        ,NUM_LEVELS_MOVE_3D_V                &
                                        ,TASK_UPDATE_SPECS(N)                &  !<-- Linked list of nest task update region specs
                                        ,MOVING_CHILD_UPDATE(N)              &  !<-- Composite H/V update data for nest task N
                                          )
            ENDIF
!
            parent_update_moving_tim=parent_update_moving_tim           &
                                    +(timef()-btim)
!
!-----------------------------------------------------------------------
!
          ENDIF child_moves
!
!-----------------------------------------------------------------------
!
        ENDDO child_loop_2
!
!-----------------------------------------------------------------------
!
      ENDIF moving_children
!
!-----------------------------------------------------------------------
!***  Insert clocktimes into the coupler's export state that are related
!***  to this phase.
!-----------------------------------------------------------------------
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
     CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The parent-child coupler export state
                            ,name ='parent_bookkeep_moving_tim'        &  !<-- Name of the attribute to insert
                            ,value=parent_bookkeep_moving_tim          &  !<-- moving nest bookeeping time
                            ,rc   =RC)

     CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The parent-child coupler export state
                            ,name ='parent_update_moving_tim'        &    !<-- Name of the attribute to insert
                            ,value=parent_update_moving_tim          &    !<-- moving nest update time
                            ,rc   =RC)

     CALL ESMF_AttributeSet(state=EXP_STATE                            &  !<-- The parent-child coupler export state
                            ,name ='t0_recv_move_tim'                  &  !<-- Name of the attribute to insert
                            ,value=t0_recv_move_tim                    &  !<-- task 0 time to process receive of move flag
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE RESET_WORK_PARENT(N_CHILD,N_RATIO                      &
                                  ,TIME_FLAG,PARENT_MOVED)
!
!-----------------------------------------------------------------------
!***  A parent resets its working pointers/arrays that depend on a
!***  moving child's location to get ready to generate values on
!***  that child's boundary.  This routine is not called for static
!***  nests since there is nothing to reset for them.
!***  This is an internal subroutine to PARENT_CHILD_CPL_RUN_SEND.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: N_CHILD                          &  !<-- Rank of nest in list of ALL children
                                      ,N_RATIO                             !<-- Rank of space ratio value among the moving children
!
      CHARACTER(len=*),INTENT(IN) :: TIME_FLAG                             !<-- Child to recv data from its present or future
!
      LOGICAL(kind=KLOG),INTENT(IN) :: PARENT_MOVED                        !<-- Did this parent just shift its own domain?
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IERR,INDX,N,NR,NT,NTAG,NUM_CHILD_TASKS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  If this child wants to move then deallocate parent working 
!***  arrays/pointers whose dimensions are functions of moving nests' 
!***  positions.  They will be reallocated with dimensions appropriate
!***  for the new positions.  (For static nests these parent arrays 
!***  are used over and over and are not deallocated/reallocated.)
!***  If this is the first time we have reached this point though
!***  then nothing has been allocated yet so skip the deallocation.
!***  Note however that if this parent can move then it must call
!***  this routine when it shifts since that will also mean that
!***  its children's positions have changed with respect to the
!***  parent's grid.
!-----------------------------------------------------------------------
!
      INDX=1
      IF(TIME_FLAG=='Current')INDX=2
!
      CALL DEALLOC_WORK_PARENTS(N_CHILD,TIME_FLAG)
!
!-----------------------------------------------------------------------
!***  We now compute various indices and weights needed by the parents
!***  to compute boundary data for their children.  It is here that
!***  location-dependent interpolation information is determined 
!***  regarding the parent and nests.  Parents need to call these
!***  routines only for children who have moved because this work was
!***  done once and for all for static nests in the coupler's Init
!***  step.
!-----------------------------------------------------------------------
!
      CALL PREPARE_NEST_INTERP_FACTORS(N_CHILD)
!
      CALL POINT_INTERP_DATA_TO_MEMORY(N_CHILD,TIME_FLAG)
!
!-----------------------------------------------------------------------
!***  The parent determines the new association between its tasks
!***  and those of its moving child's then sends the information
!***  to that child so the child will know exactly how to receive
!***  the new internal and boundary data from its parent when the
!***  child arrives at this point in time and executes its move.
!***  This only needs to be done when the nest has just moved, i.e.,
!***  when the time flag has switched to 'Current'.  When it goes
!***  back to 'Future' we do not need to send the information again
!***  since the nest has not moved again and thus the associations
!***  remain the same.
!-----------------------------------------------------------------------
!
      IF(TIME_FLAG=='Current'.OR.PARENT_MOVED)THEN
        CALL PARENT_SENDS_CHILD_DATA_LIMITS(N_CHILD)
      ENDIF
!
!-----------------------------------------------------------------------
!***  The parent determines the child's boundary topography at the
!***  new location after the child moves.  This is needed to maintain
!***  hydrostatic balance when parent data is interpolated to child
!***  boundaries where the terrain is different.
!-----------------------------------------------------------------------
!
      NR=N_RATIO                 
!
      CALL PARENT_COMPUTES_CHILD_TOPO(N_CHILD                           &
                                     ,I_PARENT_SW(N_CHILD)              &
                                     ,J_PARENT_SW(N_CHILD)              &
                                     ,IM_CHILD(N_CHILD)                 &
                                     ,JM_CHILD(N_CHILD)                 &
                                     ,N_BLEND_H_CHILD(N_CHILD)          &
                                     ,NEST_FIS_ON_PARENT_BNDS(NR)%LBND1 &
                                     ,NEST_FIS_ON_PARENT_BNDS(NR)%UBND1 &
                                     ,NEST_FIS_ON_PARENT_BNDS(NR)%LBND2 &
                                     ,NEST_FIS_ON_PARENT_BNDS(NR)%UBND2 &
                                     ,NEST_FIS_ON_PARENT(NR)%DATA       &
                                       )
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE RESET_WORK_PARENT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE COMPUTE_SEND_NEST_BC_DATA(N_CHILD,TIME_FLAG)
!
!-----------------------------------------------------------------------
!***  A parent generates and sends boundary data to a child.
!***  This is an internal subroutine to PARENT_CHILD_CPL_RUN_SEND.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: N_CHILD                             !<-- Compute/send this child's boundary conditions
!
      CHARACTER(len=*),INTENT(IN) :: TIME_FLAG                             !<-- Child to recv data from its present or future
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IERR,INDX2,N,NT,NTAG,NUM_CHILD_TASKS
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME) :: PD_V
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
! 
      N=N_CHILD
!
!-----------------------------------------------------------------------
!***  Select the appropriate part of the working array depending on
!***  whether we are now concerned with children's boundaries for
!***  their current time or from their future.
!-----------------------------------------------------------------------
!
      IF(TIME_FLAG=='Future')THEN
        INDX2=1
      ELSEIF(TIME_FLAG=='Current')THEN
        INDX2=2
      ENDIF
!
!-----------------------------------------------------------------------
!***  Before parents can generate new boundary data for their children
!***  we must check to be sure the previous set of ISend's from the
!***  parent tasks to the children's boundary tasks have completed.
!-----------------------------------------------------------------------
!
      btim=timef()
!
!-------------
!***  South H
!-------------
!
      IF(NUM_TASKS_SEND_H_S(N)>0)THEN                                      !<-- Parent task has Sbndry H data to send to child tasks?
        DO NT=1,NUM_TASKS_SEND_H_S(N)
          CALL MPI_WAIT(HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
!-------------
!***  South V
!-------------
!
      IF(NUM_TASKS_SEND_V_S(N)>0)THEN                                      !<-- Parent task has Sbndry V data to send to child tasks?
        DO NT=1,NUM_TASKS_SEND_V_S(N)
          CALL MPI_WAIT(HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
!-------------
!***  North H
!-------------
!
      IF(NUM_TASKS_SEND_H_N(N)>0)THEN                                      !<-- Parent task has Nbndry H data to send to child tasks?
        DO NT=1,NUM_TASKS_SEND_H_N(N)
          CALL MPI_WAIT(HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
!-------------
!***  North V
!-------------
!
      IF(NUM_TASKS_SEND_V_N(N)>0)THEN                                      !<-- Parent task has Nbndry V data to send to child tasks?
        DO NT=1,NUM_TASKS_SEND_V_N(N)
          CALL MPI_WAIT(HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
!------------
!***  West H
!------------
!
      IF(NUM_TASKS_SEND_H_W(N)>0)THEN                                      !<-- Parent task has Wbndry H data to send to child tasks?
        DO NT=1,NUM_TASKS_SEND_H_W(N)
          CALL MPI_WAIT(HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV(NT)       &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
!------------
!***  West V
!------------
!
      IF(NUM_TASKS_SEND_V_W(N)>0)THEN                                      !<-- Parent task has Wbndry V data to send to child tasks?
        DO NT=1,NUM_TASKS_SEND_V_W(N)
          CALL MPI_WAIT(HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV(NT)       &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
!------------
!***  East H
!------------
!
      IF(NUM_TASKS_SEND_H_E(N)>0)THEN                                      !<-- Parent task has Ebndry H data to send to child tasks?
        DO NT=1,NUM_TASKS_SEND_H_E(N)
          CALL MPI_WAIT(HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV(NT)       &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
!------------
!***  East V
!------------
!
      IF(NUM_TASKS_SEND_V_E(N)>0)THEN                                      !<-- Parent task has Ebndry V data to send to child tasks?
        DO NT=1,NUM_TASKS_SEND_V_E(N)
          CALL MPI_WAIT(HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV(NT)       &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
      cpl2_wait_tim=cpl2_wait_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  The parents can now compute the new surface pressure on the
!***  nests' boundary points (overwriting the previous values).
!-----------------------------------------------------------------------
!
      NUM_CHILD_TASKS=FTASKS_DOMAIN(MY_CHILDREN_ID(N))
!
!--------
!***  PD
!--------
!
      btim=timef()
      CALL PARENT_UPDATE_CHILD_PSFC(FIS,PD,T,TRACERS(:,:,:,INDX_Q)      &  !<-- Native parent values
                                   ,PT,PDTOP                            &  !<-- Domain PT and PDTOP
                                   ,SG1,SG2                             &  !<-- General vertical structure (shared by all domains)
                                   ,IMS,IME,JMS,JME                     &  !<-- Parent task subdomain lateral memory dimensions
                                   ,LM                                  &  !<-- # of model layers
!
                                   ,NUM_CHILD_TASKS                     &  !<-- # of fcst tasks on child N
                                   ,CTASK_LIMITS(N)%LIMITS              &  !<-- Integration limits on each task of child N
!
                                   ,FIS_CHILD_SOUTH(N)%TASKS            &  !<-- Sfc geopotential on Sbndry points on child tasks
                                   ,FIS_CHILD_NORTH(N)%TASKS            &  !<-- Sfc geopotential on Nbndry points on child tasks
                                   ,FIS_CHILD_WEST(N)%TASKS             &  !<-- Sfc geopotential on Wbndry points on child tasks
                                   ,FIS_CHILD_EAST(N)%TASKS             &  !<-- Sfc geopotential on Ebndry points on child tasks
! 
                                   ,CHILDTASK_BNDRY_H_RANKS(N)%SOUTH    &  !<-- Ranks of child N's fcst tasks on its Sbndry
                                   ,CHILDTASK_BNDRY_H_RANKS(N)%NORTH    &  !<-- Ranks of child N's fcst tasks on its Nbndry
                                   ,CHILDTASK_BNDRY_H_RANKS(N)%WEST     &  !<-- Ranks of child N's fcst tasks on its Wbndry
                                   ,CHILDTASK_BNDRY_H_RANKS(N)%EAST     &  !<-- Ranks of child N's fcst tasks on its Ebndry
!
                                   ,NUM_TASKS_SEND_H_S(N)               &  !<-- # of child tasks with south boundary segments
                                   ,NUM_TASKS_SEND_H_N(N)               &  !<-- # of child tasks with north boundary segments
                                   ,NUM_TASKS_SEND_H_W(N)               &  !<-- # of child tasks with west boundary segments
                                   ,NUM_TASKS_SEND_H_E(N)               &  !<-- # of child tasks with east boundary segments
!
                                   ,PARENT_4_INDICES_H(N)%I_INDX_SBND   &  !<-- Parent I's west and east of each child Sbndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_NBND   &  !<-- Parent I's west and east of each child Nbndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_WBND   &  !<-- Parent I's west and east of each child Wbndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_EBND   &  !<-- Parent I's west and east of each child Ebndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_SBND   &  !<-- Parent J's south and north of each child Sbndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_NBND   &  !<-- Parent J's south and north of each child Nbndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_WBND   &  !<-- Parent J's south and north of each child Wbndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_EBND   &  !<-- Parent J's south and north of each child Ebndry point
!
                             ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH            &  !<-- Starting I on each south boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH            &  !<-- Ending I on each south boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER   &  !<-- Ending I on each south boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_LO_NORTH            &  !<-- Starting I on each north boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_NORTH            &  !<-- Ending I on each north boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER   &  !<-- Ending I on each north boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_LO_WEST             &  !<-- Starting J on each west boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_WEST             &  !<-- Ending J on each west boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER    &  !<-- Ending J on each west boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_LO_EAST             &  !<-- Starting J on each east boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_EAST             &  !<-- Ending J on each east boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER    &  !<-- Ending J on each east boundary child task
!
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND   &  !<-- Bilinear interpolation wgts of the four parent
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND   &  !    points surrounding each child bndry point
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND   &  !    on each side of the child boundary.
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND   &  !
!
                                   ,N_BLEND_H_CHILD(N)                  &  !<-- Width of boundary blending region for mass points
                                   ,IM_CHILD(N)                         &  !<-- East-west points on child domain
                                   ,JM_CHILD(N)                         &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             Input
!                                                                          --------------
!                                                                             Output
!                                                                               |   
!                                                                               v   
                                   ,CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS  &  !<-- 1-D H-pt Sbndry datastring to be sent by parent to child
                                   ,CHILD_BOUND_H_NORTH(N,INDX2)%TASKS  &  !<-- 1-D H-pt Nbndry datastring to be sent by parent to child
                                   ,CHILD_BOUND_H_WEST(N,INDX2)%TASKS   &  !<-- 1-D H-pt Wbndry datastring to be sent by parent to child
                                   ,CHILD_BOUND_H_EAST(N,INDX2)%TASKS   &  !<-- 1-D H-pt Ebndry datastring to be sent by parent to child
!
                                   ,PD_B_SOUTH(N)%TASKS                 &  !<-- Updated sigma domain pressure (Pa) on nest bndry points
                                   ,PD_B_NORTH(N)%TASKS                 &  !    for all four sides of nest N's boundary.
                                   ,PD_B_WEST(N)%TASKS                  &  !
                                   ,PD_B_EAST(N)%TASKS )                   !<--
!
!-----------------------------------------------------------------------
!***  Now compute the new mass variable values in the columns above
!***  the nest boundary points.
!-----------------------------------------------------------------------
!
!-----------------
!***  Temperature
!-----------------
!
      CALL PARENT_UPDATE_CHILD_BNDRY(T                                  &  !<-- Parent sensible temperature (K)
                                    ,'T'                                &  !<-- Parent interpolating T
!
                                    ,PD,PT,PDTOP                        &  !<-- Parent PD; domain PT and PDTOP
                                    ,PSGML1,SGML2,SG1,SG2               &  !<-- General vertical structure (shared by all domains)
!
                                    ,PD_B_SOUTH(N)%TASKS                &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                    ,PD_B_NORTH(N)%TASKS                &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                    ,PD_B_WEST(N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                    ,PD_B_EAST(N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                    ,IMS,IME,JMS,JME                    &  !<-- Parent task subdomain lateral memory dimensions
                                    ,LM                                 &  !<-- # of model layers
                                    ,0                                  &  !<-- # of rows to ignore on north/east nest boundaries
!
                                    ,NUM_TASKS_SEND_H_S(N)              &  !<-- # of child tasks with south boundary segments
                                    ,NUM_TASKS_SEND_H_N(N)              &  !<-- # of child tasks with north boundary segments
                                    ,NUM_TASKS_SEND_H_W(N)              &  !<-- # of child tasks with west boundary segments
                                    ,NUM_TASKS_SEND_H_E(N)              &  !<-- # of child tasks with east boundary segments
!
                                   ,PARENT_4_INDICES_H(N)%I_INDX_SBND   &  !<-- Parent I's west and east of each child S bndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_NBND   &  !<-- Parent I's west and east of each child N bndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_WBND   &  !<-- Parent I's west and east of each child W bndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_EBND   &  !<-- Parent I's west and east of each child E bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_SBND   &  !<-- Parent J's south and north of each child S bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_NBND   &  !<-- Parent J's south and north of each child N bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_WBND   &  !<-- Parent J's south and north of each child W bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_EBND   &  !<-- Parent J's south and north of each child E bndry point
!
                             ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH            &  !<-- Starting I on each south boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH            &  !<-- Ending I on each south boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER   &  !<-- Ending I for transfer to child on each Sbndry child task
                             ,CHILDTASK_H_SAVE(N)%I_LO_NORTH            &  !<-- Starting I on each north boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_NORTH            &  !<-- Ending I on each north boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER   &  !<-- Ending I for transfer to child on each Nbndry child task
                             ,CHILDTASK_H_SAVE(N)%J_LO_WEST             &  !<-- Starting J on each west boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_WEST             &  !<-- Ending J on each west boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER    &  !<-- Ending J for transfer to child on each Wbndry child task
                             ,CHILDTASK_H_SAVE(N)%J_LO_EAST             &  !<-- Starting J on each east boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_EAST             &  !<-- Ending J on each east boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER    &  !<-- Ending J for transfer to child on each Ebndry child task
!
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND   &  !<-- Bilinear interpolation wgts of the four parent
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND   &  !    points surrounding each child bndry point.
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND   &  !
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND   &  !<--
!
                                    ,N_BLEND_H_CHILD(N)                 &  !<-- Width of boundary blending region
                                    ,IM_CHILD(N)                        &  !<-- East-west points on child domain
                                    ,JM_CHILD(N)                        &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             Input
!                                                                          --------------
!                                                                             Output
!                                                                               |   
!                                                                               v   
                                    ,T_B_SOUTH(N)%TASKS                 &  !<-- Updated sensible temperature (K) on nest bndry points.
                                    ,T_B_NORTH(N)%TASKS                 &  !
                                    ,T_B_WEST(N)%TASKS                  &  !
                                    ,T_B_EAST(N)%TASKS )                   !<-- 
!
!-----------------
!***  Specific humidity
!-----------------
!
      CALL PARENT_UPDATE_CHILD_BNDRY(TRACERS(:,:,:,INDX_Q)              &  !<-- Parent specific humidity (kg/kg) 
                                    ,'Q'                                &  !<-- Parent interpolating Q
!
                                    ,PD,PT,PDTOP                        &  !<-- Parent PD; domain PT and PDTOP
                                    ,PSGML1,SGML2,SG1,SG2               &  !<-- General vertical structure (shared by all domains)
!
                                    ,PD_B_SOUTH(N)%TASKS                &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                    ,PD_B_NORTH(N)%TASKS                &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                    ,PD_B_WEST(N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                    ,PD_B_EAST(N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                    ,IMS,IME,JMS,JME                    &  !<-- Parent task subdomain lateral memory dimensions
                                    ,LM                                 &  !<-- # of model layers
                                    ,0                                  &  !<-- # of rows to ignore on north/east nest boundaries
!
                                    ,NUM_TASKS_SEND_H_S(N)              &  !<-- # of child tasks with south boundary segments
                                    ,NUM_TASKS_SEND_H_N(N)              &  !<-- # of child tasks with north boundary segments
                                    ,NUM_TASKS_SEND_H_W(N)              &  !<-- # of child tasks with west boundary segments
                                    ,NUM_TASKS_SEND_H_E(N)              &  !<-- # of child tasks with east boundary segments
!
                                   ,PARENT_4_INDICES_H(N)%I_INDX_SBND   &  !<-- Parent I's west and east of each child S bndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_NBND   &  !<-- Parent I's west and east of each child N bndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_WBND   &  !<-- Parent I's west and east of each child W bndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_EBND   &  !<-- Parent I's west and east of each child E bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_SBND   &  !<-- Parent J's south and north of each child S bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_NBND   &  !<-- Parent J's south and north of each child N bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_WBND   &  !<-- Parent J's south and north of each child W bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_EBND   &  !<-- Parent J's south and north of each child E bndry point
!
                             ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH            &  !<-- Starting I on each south boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH            &  !<-- Ending I on each south boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER   &  !<-- Ending I for transfer to child on each Sbndry child task
                             ,CHILDTASK_H_SAVE(N)%I_LO_NORTH            &  !<-- Starting I on each north boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_NORTH            &  !<-- Ending I on each north boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER   &  !<-- Ending I for transfer to child on each Nbndry child task
                             ,CHILDTASK_H_SAVE(N)%J_LO_WEST             &  !<-- Starting J on each west boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_WEST             &  !<-- Ending J on each west boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER    &  !<-- Ending J for transfer to child on each Wbndry child task
                             ,CHILDTASK_H_SAVE(N)%J_LO_EAST             &  !<-- Starting J on each east boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_EAST             &  !<-- Ending J on each east boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER    &  !<-- Ending J for transfer to child on each Ebndry child task
!
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND   &  !<-- Bilinear interpolation wgts of the four parent
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND   &  !    points surrounding each child bndry point.
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND   &  !
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND   &  !<--
!
                                    ,N_BLEND_H_CHILD(N)                 &  !<-- Width of boundary blending region
                                    ,IM_CHILD(N)                        &  !<-- East-west points on child domain
                                    ,JM_CHILD(N)                        &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             Input
!                                                                          --------------
!                                                                             Output
!                                                                               |   
!                                                                               v   
                                    ,Q_B_SOUTH(N)%TASKS                 &  !<-- Updated specific humidity (kg/kg) on nest bndry points.
                                    ,Q_B_NORTH(N)%TASKS                 &  !
                                    ,Q_B_WEST(N)%TASKS                  &  !
                                    ,Q_B_EAST(N)%TASKS )                   !<-- 
!
!-----------------------------------------------------------------------
!
!----------------------
!***  Cloud Condensate 
!----------------------
!
      CALL PARENT_UPDATE_CHILD_BNDRY(TRACERS(:,:,:,INDX_CW)             &  !<-- Parent cloud condensate (kg/kg)
                                    ,'CW'                               &  !<-- Parent interpolating cloud condensate
!
                                    ,PD,PT,PDTOP                        &  !<-- Parent PD; domain PT and PDTOP
                                    ,PSGML1,SGML2,SG1,SG2               &  !<-- General vertical structure (shared by all domains)
!
                                    ,PD_B_SOUTH(N)%TASKS                &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                    ,PD_B_NORTH(N)%TASKS                &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                    ,PD_B_WEST(N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                    ,PD_B_EAST(N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                    ,IMS,IME,JMS,JME                    &  !<-- Parent task subdomain lateral memory dimensions
                                    ,LM                                 &  !<-- # of model layers
                                    ,0                                  &  !<-- # of rows to ignore on north/east nest boundaries
!
                                    ,NUM_TASKS_SEND_H_S(N)              &  !<-- # of child tasks with south boundary segments
                                    ,NUM_TASKS_SEND_H_N(N)              &  !<-- # of child tasks with north boundary segments
                                    ,NUM_TASKS_SEND_H_W(N)              &  !<-- # of child tasks with west boundary segments
                                    ,NUM_TASKS_SEND_H_E(N)              &  !<-- # of child tasks with east boundary segments
!
                                   ,PARENT_4_INDICES_H(N)%I_INDX_SBND   &  !<-- Parent I's west and east of each child S bndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_NBND   &  !<-- Parent I's west and east of each child N bndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_WBND   &  !<-- Parent I's west and east of each child W bndry point
                                   ,PARENT_4_INDICES_H(N)%I_INDX_EBND   &  !<-- Parent I's west and east of each child E bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_SBND   &  !<-- Parent J's south and north of each child S bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_NBND   &  !<-- Parent J's south and north of each child N bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_WBND   &  !<-- Parent J's south and north of each child W bndry point
                                   ,PARENT_4_INDICES_H(N)%J_INDX_EBND   &  !<-- Parent J's south and north of each child E bndry point
!
                             ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH            &  !<-- Starting I on each south boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH            &  !<-- Ending I on each south boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER   &  !<-- Ending I for transfer to child on each Sbndry child task
                             ,CHILDTASK_H_SAVE(N)%I_LO_NORTH            &  !<-- Starting I on each north boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_NORTH            &  !<-- Ending I on each north boundary child task
                             ,CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER   &  !<-- Ending I for transfer to child on each Nbndry child task
                             ,CHILDTASK_H_SAVE(N)%J_LO_WEST             &  !<-- Starting J on each west boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_WEST             &  !<-- Ending J on each west boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER    &  !<-- Ending J for transfer to child on each Wbndry child task
                             ,CHILDTASK_H_SAVE(N)%J_LO_EAST             &  !<-- Starting J on each east boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_EAST             &  !<-- Ending J on each east boundary child task
                             ,CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER    &  !<-- Ending J for transfer to child on each Ebndry child task
!
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_SBND   &  !<-- Bilinear interpolation wgts of the four parent
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_NBND   &  !    points surrounding each child bndry point.
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_WBND   &  !
                                  ,PARENT_4_WEIGHTS_H(N)%WEIGHTS_EBND   &  !<--
!
                                    ,N_BLEND_H_CHILD(N)                 &  !<-- Width of boundary blending region
                                    ,IM_CHILD(N)                        &  !<-- East-west points on child domain
                                    ,JM_CHILD(N)                        &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             Input
!                                                                          --------------
!                                                                             Output
!                                                                               |   
!                                                                               v   
                                    ,CW_B_SOUTH(N)%TASKS                &  !<-- Updated cloud condensate (kg/kg) on nest bndry points.
                                    ,CW_B_NORTH(N)%TASKS                &  !
                                    ,CW_B_WEST(N)%TASKS                 &  !
                                    ,CW_B_EAST(N)%TASKS )                  !<-- 
!
!-----------------------------------------------------------------------
!***  Before we can let the parent proceed to compute wind component
!***  updates on the nests' boundaries, we must generate the pressure
!***  on the parent's V points and at the nest boundary V points.
!***  This is begun by simple 4-pt averaging.  
!-----------------------------------------------------------------------
!
      CALL PRESSURE_ON_NEST_BNDRY_V(PD                                  &  !<-- Sigma domain pressure (Pa) on parent mass points
                                   ,IMS,IME,JMS,JME                     &  !<-- Memory dimensions of PD
!
                                   ,PD_B_SOUTH(N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Sbndry mass points
                                   ,PD_B_NORTH(N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Nbndry mass points
                                   ,PD_B_WEST (N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Wbndry mass points
                                   ,PD_B_EAST (N)%TASKS                 &  !<-- Sigma domain pressure (Pa) on nest Ebndry mass points
! 
                                   ,NUM_TASKS_SEND_V_S(N)               &  !<-- # of child tasks with south boundary segments on V
                                   ,NUM_TASKS_SEND_V_N(N)               &  !<-- # of child tasks with north boundary segments on V
                                   ,NUM_TASKS_SEND_V_W(N)               &  !<-- # of child tasks with west boundary segments on V
                                   ,NUM_TASKS_SEND_V_E(N)               &  !<-- # of child tasks with east boundary segments on V
!
                                   ,CHILDTASK_V_SAVE(N)%I_LO_SOUTH      &  !<-- Starting I on each south V boundary child task 
                                   ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH      &  !<-- Ending I on each south V boundary child task
                                   ,CHILDTASK_V_SAVE(N)%I_LO_NORTH      &  !<-- Starting I on each north V boundary child task
                                   ,CHILDTASK_V_SAVE(N)%I_HI_NORTH      &  !<-- Ending I on each north V boundary child task
                                   ,CHILDTASK_V_SAVE(N)%J_LO_WEST       &  !<-- Starting J on each west V boundary child task
                                   ,CHILDTASK_V_SAVE(N)%J_HI_WEST       &  !<-- Ending J on each west V boundary child task
                                   ,CHILDTASK_V_SAVE(N)%J_LO_EAST       &  !<-- Starting J on each east V boundary child task
                                   ,CHILDTASK_V_SAVE(N)%J_HI_EAST       &  !<-- Ending J on each east V boundary child task
!
                                   ,CHILDTASK_H_SAVE(N)%I_LO_SOUTH      &  !<-- Starting I on each south H boundary child task
                                   ,CHILDTASK_H_SAVE(N)%I_HI_SOUTH      &  !<-- Ending I on each south H boundary child task
                                   ,CHILDTASK_H_SAVE(N)%I_LO_NORTH      &  !<-- Starting I on each north H boundary child task
                                   ,CHILDTASK_H_SAVE(N)%I_HI_NORTH      &  !<-- Ending I on each north H boundary child task
                                   ,CHILDTASK_H_SAVE(N)%J_LO_WEST       &  !<-- Starting J on each west H boundary child task
                                   ,CHILDTASK_H_SAVE(N)%J_HI_WEST       &  !<-- Ending J on each west H boundary child task
                                   ,CHILDTASK_H_SAVE(N)%J_LO_EAST       &  !<-- Starting J on each east H boundary child task
                                   ,CHILDTASK_H_SAVE(N)%J_HI_EAST       &  !<-- Ending J on each east H boundary child task
!
                                   ,N_BLEND_H_CHILD(N)                  &  !<-- H rows in nests' boundary regions 
                                   ,N_BLEND_V_CHILD(N)                  &  !<-- V rows in nests' boundary regions
                                   ,IM_CHILD(N)                         &  !<-- East-west points on child domain
                                   ,JM_CHILD(N)                         &  !<-- North-south points on child domain
!
                                   ,INC_FIX(N)                          &  !<-- Increment used to select nest tasks for averaging H to V
!                                                                               ^
!                                                                               |
!                                                                             Input
!                                                                          --------------
!                                                                             Output
!                                                                               |   
!                                                                               v   
                                   ,PD_V                                &  !<-- Sigma domain pressure (Pa) on parent V points
                                   ,PD_B_SOUTH_V(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Sbndry V points
                                   ,PD_B_NORTH_V(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Nbndry V points
                                   ,PD_B_WEST_V(N)%TASKS                &  !<-- Sigma domain pressure (Pa) on nest Wbndry V points
                                   ,PD_B_EAST_V(N)%TASKS )                 !<-- Sigma domain pressure (Pa) on nest Ebndry V points
!
!-----------------------------------------------------------------------
!
!----------------------
!***  U Wind Component 
!----------------------
!
      CALL PARENT_UPDATE_CHILD_BNDRY(U                                  &  !<-- Parent U wind component (m/s)
                                    ,'U'                                &  !<-- Parent interpolating U wind component
!
                                    ,PD_V,PT,PDTOP                      &  !<-- Parent PD; domain PT and PDTOP
                                    ,PSGML1,SGML2,SG1,SG2               &  !<-- General vertical structure (shared by all domains)
!
                                    ,PD_B_SOUTH_V(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                    ,PD_B_NORTH_V(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                    ,PD_B_WEST_V(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                    ,PD_B_EAST_V(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                    ,IMS,IME,JMS,JME                    &  !<-- Parent task subdomain lateral memory dimensions
                                    ,LM                                 &  !<-- # of model layers
                                    ,1                                  &  !<-- # of rows to ignore on north/east nest boundaries
!
                                    ,NUM_TASKS_SEND_V_S(N)              &  !<-- # of child tasks with south boundary segments
                                    ,NUM_TASKS_SEND_V_N(N)              &  !<-- # of child tasks with north boundary segments
                                    ,NUM_TASKS_SEND_V_W(N)              &  !<-- # of child tasks with west boundary segments
                                    ,NUM_TASKS_SEND_V_E(N)              &  !<-- # of child tasks with east boundary segments
!
                                   ,PARENT_4_INDICES_V(N)%I_INDX_SBND   &  !<-- Parent I's west and east of each child S bndry point
                                   ,PARENT_4_INDICES_V(N)%I_INDX_NBND   &  !<-- Parent I's west and east of each child N bndry point
                                   ,PARENT_4_INDICES_V(N)%I_INDX_WBND   &  !<-- Parent I's west and east of each child W bndry point
                                   ,PARENT_4_INDICES_V(N)%I_INDX_EBND   &  !<-- Parent I's west and east of each child E bndry point
                                   ,PARENT_4_INDICES_V(N)%J_INDX_SBND   &  !<-- Parent J's south and north of each child S bndry point
                                   ,PARENT_4_INDICES_V(N)%J_INDX_NBND   &  !<-- Parent J's south and north of each child N bndry point
                                   ,PARENT_4_INDICES_V(N)%J_INDX_WBND   &  !<-- Parent J's south and north of each child W bndry point
                                   ,PARENT_4_INDICES_V(N)%J_INDX_EBND   &  !<-- Parent J's south and north of each child E bndry point
!
                             ,CHILDTASK_V_SAVE(N)%I_LO_SOUTH            &  !<-- Starting I on each south boundary child task
                             ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH            &  !<-- Ending I on each south boundary child task
                             ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH_TRANSFER   &  !<-- Not relevant for V points
                             ,CHILDTASK_V_SAVE(N)%I_LO_NORTH            &  !<-- Starting I on each north boundary child task
                             ,CHILDTASK_V_SAVE(N)%I_HI_NORTH            &  !<-- Ending I on each north boundary child task
                             ,CHILDTASK_V_SAVE(N)%I_HI_NORTH_TRANSFER   &  !<-- Not relevant for V points
                             ,CHILDTASK_V_SAVE(N)%J_LO_WEST             &  !<-- Starting J on each west boundary child task
                             ,CHILDTASK_V_SAVE(N)%J_HI_WEST             &  !<-- Ending J on each west boundary child task
                             ,CHILDTASK_V_SAVE(N)%J_HI_WEST_TRANSFER    &  !<-- Not relevant for V points
                             ,CHILDTASK_V_SAVE(N)%J_LO_EAST             &  !<-- Starting J on each east boundary child task
                             ,CHILDTASK_V_SAVE(N)%J_HI_EAST             &  !<-- Ending J on each east boundary child task
                             ,CHILDTASK_V_SAVE(N)%J_HI_EAST_TRANSFER    &  !<-- Not relevant for V points
!
                                  ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_SBND   &  !<-- Bilinear interpolation wgts of the four parent
                                  ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_NBND   &  !    points surrounding each child bndry point.
                                  ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_WBND   &  !
                                  ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_EBND   &  !<--
!
                                    ,N_BLEND_V_CHILD(N)                 &  !<-- Width of boundary blending region
                                    ,IM_CHILD(N)                        &  !<-- East-west points on child domain
                                    ,JM_CHILD(N)                        &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             Input
!                                                                          --------------
!                                                                             Output
!                                                                               |   
!                                                                               v   
                                    ,U_B_SOUTH(N)%TASKS                 &  !<-- Updated U wind component (m/s) on nest bndry points.
                                    ,U_B_NORTH(N)%TASKS                 &  !
                                    ,U_B_WEST(N)%TASKS                  &  !
                                    ,U_B_EAST(N)%TASKS )                   !<-- 
!
!------------
!***  V Wind
!------------
!
      CALL PARENT_UPDATE_CHILD_BNDRY(V                                  &  !<-- Parent V wind component (m/s)
                                    ,'V'                                &  !<-- Parent interpolating V wind component
!
                                    ,PD_V,PT,PDTOP                      &  !<-- Parent PD; domain PT and PDTOP
                                    ,PSGML1,SGML2,SG1,SG2               &  !<-- General vertical structure (shared by all domains)
!
                                    ,PD_B_SOUTH_V(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Sbndry points
                                    ,PD_B_NORTH_V(N)%TASKS              &  !<-- Sigma domain pressure (Pa) on nest Nbndry points
                                    ,PD_B_WEST_V(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Wbndry points
                                    ,PD_B_EAST_V(N)%TASKS               &  !<-- Sigma domain pressure (Pa) on nest Ebndry points
!
                                    ,IMS,IME,JMS,JME                    &  !<-- Parent task subdomain lateral memory dimensions
                                    ,LM                                 &  !<-- # of model layers
                                    ,1                                  &  !<-- # of rows to ignore on north/east nest boundaries
!
                                    ,NUM_TASKS_SEND_V_S(N)              &  !<-- # of child tasks with south boundary segments
                                    ,NUM_TASKS_SEND_V_N(N)              &  !<-- # of child tasks with north boundary segments
                                    ,NUM_TASKS_SEND_V_W(N)              &  !<-- # of child tasks with west boundary segments
                                    ,NUM_TASKS_SEND_V_E(N)              &  !<-- # of child tasks with east boundary segments
!
                                   ,PARENT_4_INDICES_V(N)%I_INDX_SBND   &  !<-- Parent I's west and east of each child S bndry point
                                   ,PARENT_4_INDICES_V(N)%I_INDX_NBND   &  !<-- Parent I's west and east of each child N bndry point
                                   ,PARENT_4_INDICES_V(N)%I_INDX_WBND   &  !<-- Parent I's west and east of each child W bndry point
                                   ,PARENT_4_INDICES_V(N)%I_INDX_EBND   &  !<-- Parent I's west and east of each child E bndry point
                                   ,PARENT_4_INDICES_V(N)%J_INDX_SBND   &  !<-- Parent J's south and north of each child S bndry point
                                   ,PARENT_4_INDICES_V(N)%J_INDX_NBND   &  !<-- Parent J's south and north of each child N bndry point
                                   ,PARENT_4_INDICES_V(N)%J_INDX_WBND   &  !<-- Parent J's south and north of each child W bndry point
                                   ,PARENT_4_INDICES_V(N)%J_INDX_EBND   &  !<-- Parent J's south and north of each child E bndry point
!
                             ,CHILDTASK_V_SAVE(N)%I_LO_SOUTH            &  !<-- Starting I on each south boundary child task
                             ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH            &  !<-- Ending I on each south boundary child task
                             ,CHILDTASK_V_SAVE(N)%I_HI_SOUTH_TRANSFER   &  !<-- Not relevant for V points
                             ,CHILDTASK_V_SAVE(N)%I_LO_NORTH            &  !<-- Starting I on each north boundary child task
                             ,CHILDTASK_V_SAVE(N)%I_HI_NORTH            &  !<-- Ending I on each north boundary child task
                             ,CHILDTASK_V_SAVE(N)%I_HI_NORTH_TRANSFER   &  !<-- Not relevant for V points
                             ,CHILDTASK_V_SAVE(N)%J_LO_WEST             &  !<-- Starting J on each west boundary child task
                             ,CHILDTASK_V_SAVE(N)%J_HI_WEST             &  !<-- Ending J on each west boundary child task
                             ,CHILDTASK_V_SAVE(N)%J_HI_WEST_TRANSFER    &  !<-- Not relevant for V points
                             ,CHILDTASK_V_SAVE(N)%J_LO_EAST             &  !<-- Starting J on each east boundary child task
                             ,CHILDTASK_V_SAVE(N)%J_HI_EAST             &  !<-- Ending J on each east boundary child task
                             ,CHILDTASK_V_SAVE(N)%J_HI_EAST_TRANSFER    &  !<-- Not relevant for V points
!
                                  ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_SBND   &  !<-- Bilinear interpolation wgts of the four parent
                                  ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_NBND   &  !    points surrounding each child bndry point.
                                  ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_WBND   &  !
                                  ,PARENT_4_WEIGHTS_V(N)%WEIGHTS_EBND   &  !<--
!
                                    ,N_BLEND_V_CHILD(N)                 &  !<-- Width of boundary blending region
                                    ,IM_CHILD(N)                        &  !<-- East-west points on child domain
                                    ,JM_CHILD(N)                        &  !<-- North-south points on child domain
!                                                                               ^
!                                                                               |
!                                                                             Input
!                                                                          --------------
!                                                                             Output
!                                                                               |   
!                                                                               v   
                                    ,V_B_SOUTH(N)%TASKS                 &  !<-- Updated V wind component (m/s) on nest bndry points.
                                    ,V_B_NORTH(N)%TASKS                 &  !
                                    ,V_B_WEST(N)%TASKS                  &  !
                                    ,V_B_EAST(N)%TASKS )                   !<-- 
!
      cpl2_comp_tim=cpl2_comp_tim+(timef()-btim)
!
!-----------------------------------------------------------------------
!***  Parent tasks send data directly to child tasks whose boundary
!***  points the parent tasks contain. 
!-----------------------------------------------------------------------
!
      IF(TIME_FLAG=='Current')THEN
        NSTEP_CHILD_RECV(N)=(NTIMESTEP+1)*TIME_RATIO_MY_CHILDREN(N)        !<-- Child "N" is waiting at this timestep to recv its data
      ELSEIF(TIME_FLAG=='Future')THEN
        NSTEP_CHILD_RECV(N)=NTIMESTEP*TIME_RATIO_MY_CHILDREN(N)            !<-- Child "N" is waiting at this timestep to recv its data
      ENDIF
!
!-------------
!***  South H
!-------------
!
      NTAG=NSTEP_CHILD_RECV(N)+101                                         !<-- Add 101 to obtain a unique South H tag
!
      IF(NUM_TASKS_SEND_H_S(N)>0)THEN                                      !<-- Parent task has Sbndry H data to send to child tasks?     
!
        DO NT=1,NUM_TASKS_SEND_H_S(N)                                      !<-- Send to all appropriate child tasks
!
!         call date_and_time(values=values)
!         write(0,221)n,nt,childtask_bndry_h_ranks(n)%south(nt),values(5),values(6),values(7),values(8)
! 221     format(' Ready to send South_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
!     write(0,*)' ready to send South_H to child #',n,' task #',nt,' id ',childtask_bndry_h_ranks(n)%south(nt)
          btim=timef()
          CALL MPI_ISEND(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA    &  !<-- Child south boundary H data on child task NT
                        ,WORDS_BOUND_H_SOUTH(N)%TASKS(NT)               &  !<-- # of words in the data string
                        ,MPI_REAL                                       &  !<-- Datatype
                        ,CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NT)           &  !<-- Local rank of child to recv data
                        ,NTAG                                           &  !<-- MPI tag
                        ,COMM_TO_MY_CHILDREN(N)                         &  !<-- MPI communicator
                        ,HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV(NT)     &  !<-- Handle for ISend to child N's task NT
                        ,IERR )
          cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!     write(0,*)' isent South_H to child #',n,' task #',nt,' id ',childtask_bndry_h_ranks(n)%south(nt)
!     write(0,*)' # of words=',words_bound_h_south(n)%tasks(nt),' ntag=',ntag,' comm=',comm_to_my_children(n)
!
!         call date_and_time(values=values)
!         write(0,124)n,nt,childtask_bndry_h_ranks(n)%south(nt),values(5),values(6),values(7),values(8)
! 124     format(' Sent South_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        ENDDO 
!
      ENDIF
!
!-------------
!***  South V
!-------------
!
      NTAG=NSTEP_CHILD_RECV(N)+102                                         !<-- Add 102 to obtain a unique South V tag
!
      IF(NUM_TASKS_SEND_V_S(N)>0)THEN                                      !<-- Parent task has Sbndry V data to send to child tasks?     
        DO NT=1,NUM_TASKS_SEND_V_S(N)                                      !<-- Send to all appropriate child tasks
!
!         call date_and_time(values=values)
!         write(0,125)n,nt,childtask_bndry_v_ranks(n)%south(nt),values(5),values(6),values(7),values(8)
! 125     format(' Ready to send South_V to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          btim=timef()
          CALL MPI_ISEND(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(NT)%DATA    &  !<-- Child south boundary V data on child task NT
                        ,WORDS_BOUND_V_SOUTH(N)%TASKS(NT)               &  !<-- # of words in the data string
                        ,MPI_REAL                                       &  !<-- Datatype
                        ,CHILDTASK_BNDRY_V_RANKS(N)%SOUTH(NT)           &  !<-- Local rank of child to recv data
                        ,NTAG                                           &  !<-- MPI tag
                        ,COMM_TO_MY_CHILDREN(N)                         &  !<-- MPI communicator
                        ,HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV(NT)     &  !<-- Handle for ISend to child N's task NT
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
!-------------
!***  North H
!-------------
!
      NTAG=NSTEP_CHILD_RECV(N)+103                                         !<-- Add 103 to obtain a unique North H tag
!
      IF(NUM_TASKS_SEND_H_N(N)>0)THEN                                      !<-- Parent task has Nbndry H data to send to child tasks?     
        DO NT=1,NUM_TASKS_SEND_H_N(N)                                      !<-- Send to all appropriate child tasks
!
!           call date_and_time(values=values)
!           write(0,127)n,nt,childtask_bndry_h_ranks(n)%north(nt),values(5),values(6),values(7),values(8)
! 127       format(' Ready to send North_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          btim=timef()
          CALL MPI_ISEND(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA    &  !<-- Child north boundary H data on child task NT
                        ,WORDS_BOUND_H_NORTH(N)%TASKS(NT)               &  !<-- # of words in the data string
                        ,MPI_REAL                                       &  !<-- Datatype
                        ,CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NT)           &  !<-- Local rank of child to recv data
                        ,NTAG                                           &  !<-- MPI tag
                        ,COMM_TO_MY_CHILDREN(N)                         &  !<-- MPI communicator
                        ,HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV(NT)     &  !<-- Handle for ISend to child N's task NT
                        ,IERR )
          cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
!           call date_and_time(values=values)
!           write(0,128)n,nt,childtask_bndry_h_ranks(n)%north(nt),values(5),values(6),values(7),values(8)
! 128       format(' Sent North_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        ENDDO 
      ENDIF
!
!-------------
!***  North V
!-------------
!
      NTAG=NSTEP_CHILD_RECV(N)+104                                         !<-- Add 104 to obtain a unique North V tag
!
      IF(NUM_TASKS_SEND_V_N(N)>0)THEN                                      !<-- Parent task has Nbndry V data to send to child tasks?     
        DO NT=1,NUM_TASKS_SEND_V_N(N)                                      !<-- Send to all appropriate child tasks
!
!           call date_and_time(values=values)
!           write(0,129)n,nt,childtask_bndry_v_ranks(n)%north(nt),values(5),values(6),values(7),values(8)
! 129       format(' Ready to send North_V to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          btim=timef()
          CALL MPI_ISEND(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA    &  !<-- Child north boundary V data on child task NT
                        ,WORDS_BOUND_V_NORTH(N)%TASKS(NT)               &  !<-- # of words in the data string
                        ,MPI_REAL                                       &  !<-- Datatype
                        ,CHILDTASK_BNDRY_V_RANKS(N)%NORTH(NT)           &  !<-- Local rank of child to recv data
                        ,NTAG                                           &  !<-- MPI tag
                        ,COMM_TO_MY_CHILDREN(N)                         &  !<-- MPI communicator
                        ,HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV(NT)     &  !<-- Handle for ISend to child N's task NT
                        ,IERR )
          cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
!           call date_and_time(values=values)
!           write(0,130)n,nt,childtask_bndry_v_ranks(n)%north(nt),values(5),values(6),values(7),values(8)
! 130       format(' Sent North_V to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        ENDDO 
      ENDIF
!
!------------
!***  West H
!------------
!
      NTAG=NSTEP_CHILD_RECV(N)+105                                         !<-- Add 105 to obtain a unique West H tag
!
      IF(NUM_TASKS_SEND_H_W(N)>0)THEN                                      !<-- Parent task has Wbndry H data to send to child tasks?     
        DO NT=1,NUM_TASKS_SEND_H_W(N)                                      !<-- Send to all appropriate child tasks
!
!           call date_and_time(values=values)
!           write(0,131)n,nt,childtask_bndry_h_ranks(n)%west(nt),values(5),values(6),values(7),values(8)
! 131       format(' Ready to send West_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
!    write(0,*)' COMPUTE_SEND_NEST_BC_DATA for H west ready to send ',WORDS_BOUND_H_WEST(N)%TASKS(NT),' words to nest #',n &
!             ,' task #',nt,' task id ',CHILDTASK_BNDRY_H_RANKS(N)%WEST(NT),' ntag=',ntag
          btim=timef()
          CALL MPI_ISEND(CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA     &  !<-- Child west boundary H data on child task NT
                        ,WORDS_BOUND_H_WEST(N)%TASKS(NT)                &  !<-- # of words in the data string
                        ,MPI_REAL                                       &  !<-- Datatype
                        ,CHILDTASK_BNDRY_H_RANKS(N)%WEST(NT)            &  !<-- Local rank of child to recv data
                        ,NTAG                                           &  !<-- MPI tag
                        ,COMM_TO_MY_CHILDREN(N)                         &  !<-- MPI communicator
                        ,HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend to child N's task NT
                        ,IERR )
          cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
!           call date_and_time(values=values)
!           write(0,132)n,nt,childtask_bndry_h_ranks(n)%west(nt),values(5),values(6),values(7),values(8)
! 132       format(' Sent West_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        ENDDO 
      ENDIF
!
!------------
!***  West V
!------------
!
      NTAG=NSTEP_CHILD_RECV(N)+106                                         !<-- Add 106 to obtain a unique West V tag
!
      IF(NUM_TASKS_SEND_V_W(N)>0)THEN                                      !<-- Parent task has Wbndry V data to send to child tasks?     
        DO NT=1,NUM_TASKS_SEND_V_W(N)                                      !<-- Send to all appropriate child tasks
!
!           call date_and_time(values=values)
!           write(0,133)n,nt,childtask_bndry_v_ranks(n)%west(nt),values(5),values(6),values(7),values(8)
! 133       format(' Ready to send West_V to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          btim=timef()
          CALL MPI_ISEND(CHILD_BOUND_V_WEST(N,INDX2)%TASKS(NT)%DATA     &  !<-- Child west boundary V data on child task NT
                        ,WORDS_BOUND_V_WEST(N)%TASKS(NT)                &  !<-- # of words in the data string
                        ,MPI_REAL                                       &  !<-- Datatype
                        ,CHILDTASK_BNDRY_V_RANKS(N)%WEST(NT)            &  !<-- Local rank of child to recv data
                        ,NTAG                                           &  !<-- MPI tag
                        ,COMM_TO_MY_CHILDREN(N)                         &  !<-- MPI communicator
                        ,HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend to child N's task NT
                        ,IERR )
          cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
!           call date_and_time(values=values)
!           write(0,134)n,nt,childtask_bndry_v_ranks(n)%west(nt),values(5),values(6),values(7),values(8)
! 134       format(' Sent West_V to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        ENDDO 
      ENDIF
!
!------------
!***  East H
!------------
!
      NTAG=NSTEP_CHILD_RECV(N)+107                                         !<-- Add 107 to obtain a unique East H tag
!
      IF(NUM_TASKS_SEND_H_E(N)>0)THEN                                      !<-- Parent task has Ebndry H data to send to child tasks?     
        DO NT=1,NUM_TASKS_SEND_H_E(N)                                      !<-- Send to all appropriate child tasks
!
!           call date_and_time(values=values)
!           write(0,135)n,nt,childtask_bndry_h_ranks(n)%east(nt),values(5),values(6),values(7),values(8)
! 135       format(' Ready to send East_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          btim=timef()
          CALL MPI_ISEND(CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA     &  !<-- Child east boundary H data on child task NT
                        ,WORDS_BOUND_H_EAST(N)%TASKS(NT)                &  !<-- # of words in the data string
                        ,MPI_REAL                                       &  !<-- Datatype
                        ,CHILDTASK_BNDRY_H_RANKS(N)%EAST(NT)            &  !<-- Local rank of child to recv data
                        ,NTAG                                           &  !<-- MPI tag
                        ,COMM_TO_MY_CHILDREN(N)                         &  !<-- MPI communicator
                        ,HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend to child N's task NT
                        ,IERR )
          cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
!           call date_and_time(values=values)
!           write(0,136)n,nt,childtask_bndry_h_ranks(n)%east(nt),values(5),values(6),values(7),values(8)
! 136       format(' Sent East_H to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        ENDDO 
      ENDIF
!
!------------
!***  East V
!------------
!
      NTAG=NSTEP_CHILD_RECV(N)+108                                         !<-- Add 108 to obtain a unique East V tag
!
      IF(NUM_TASKS_SEND_V_E(N)>0)THEN                                      !<-- Parent task has Ebndry V data to send to child tasks?     
        DO NT=1,NUM_TASKS_SEND_V_E(N)                                      !<-- Send to all appropriate child tasks
!
!           call date_and_time(values=values)
!           write(0,137)n,nt,childtask_bndry_v_ranks(n)%east(nt),values(5),values(6),values(7),values(8)
! 137       format(' Ready to send East_V to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
          btim=timef()
          CALL MPI_ISEND(CHILD_BOUND_V_EAST(N,INDX2)%TASKS(NT)%DATA     &  !<-- Child east boundary V data on child task NT
                        ,WORDS_BOUND_V_EAST(N)%TASKS(NT)                &  !<-- # of words in the data string
                        ,MPI_REAL                                       &  !<-- Datatype
                        ,CHILDTASK_BNDRY_V_RANKS(N)%EAST(NT)            &  !<-- Local rank of child to recv data
                        ,NTAG                                           &  !<-- MPI tag
                        ,COMM_TO_MY_CHILDREN(N)                         &  !<-- MPI communicator
                        ,HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend to child N's task NT
                        ,IERR )
          cpl2_send_tim=cpl2_send_tim+(timef()-btim)
!
!           call date_and_time(values=values)
!           write(0,138)n,nt,childtask_bndry_v_ranks(n)%east(nt),values(5),values(6),values(7),values(8)
! 138       format(' Sent East_V to child #',i1,' task #',i1,' id #',i3.3,' at ',i2.2,':',i2.2,':',i2.2,'.',i3.3)
!
        ENDDO 
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE COMPUTE_SEND_NEST_BC_DATA
!
!-----------------------------------------------------------------------
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
!------------------------
!***  Argument Variables
!------------------------
!
      TYPE(ESMF_CplComp) :: CPL_COMP                                       !<-- The Parent-Child Coupler Component
!
      TYPE(ESMF_State) :: IMP_STATE                                     &  !<-- The Coupler's Import State
                         ,EXP_STATE                                        !<-- The Coupler's Export State
!
      TYPE(ESMF_Clock) :: CLOCK                                            !<-- The ESMF Clock
!
      INTEGER,INTENT(OUT) :: RC_FINAL
!
!---------------------
!***  Local variables
!---------------------
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
      WRITE(0,*)' '
#ifdef ESMF_3
      IF(I_AM_A_PARENT==ESMF_TRUE)THEN
#else
      IF(I_AM_A_PARENT)THEN
#endif
        WRITE(0,*)'   Cpl2 Parent Bookkeeping for Moving Nest='         &
                  ,parent_bookkeep_moving_tim*1.e-3
        WRITE(0,*)'   Cpl2 Parent Update for Moving Nest='              &
                  ,parent_update_moving_tim*1.e-3
      ENDIF
      WRITE(0,*)' '
#ifdef ESMF_3
      IF(MY_DOMAIN_MOVES==ESMF_TRUE)THEN
#else
      IF(MY_DOMAIN_MOVES)THEN
#endif
        WRITE(0,*)'   Cpl1 Moving Nest Bookkeeping='                    &
                  ,parent_bookkeep_moving_tim*1.e-3
        WRITE(0,*)'   Cpl1 Moving Nest Update='                         &
                  ,parent_update_moving_tim*1.e-3
      ENDIF
! 
!-----------------------------------------------------------------------
!
      IF(RC_CPL_FINAL==ESMF_SUCCESS)THEN
!       WRITE(0,*)"CPL FINALIZE STEP SUCCEEDED"
      ELSE
        WRITE(0,*)"CPL FINALIZE STEP FAILED"
      ENDIF
!
      RC_FINAL=RC_CPL_FINAL
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
                                           ,DOMAIN_GRID_COMP            &  !     |
                                           ,EXP_STATE_DOMAIN            &  !     |
                                           ,FTASKS_DOMAIN               &  !     |  
                                           ,ID_PARENTS_IN               &  !     |   
                                           ,DOMAIN_ID_TO_RANK           &  !     |   
                                           ,MAX_DOMAINS                 &  !   Input 
!                                                                           -----------
                                           ,IMP_STATE_CPL_NEST          &  !   Output
                                           ,EXP_STATE_CPL_NEST          &  !     |
                                           ,PARENT_CHILD_COUPLER_COMP )    !     v
!
!-----------------------------------------------------------------------
!***  Create the Parent-Child coupler through which they will
!***  communicate.  This coupler is called by the NMM component.
!***  Move data from the DOMAIN export state into the Parent-Child
!***  coupler import state that the coupler will need in order for
!***  parents to generate data for their children and for moving
!***  nests to determine when to move.
!-----------------------------------------------------------------------
!
      USE module_DOMAIN_INTERNAL_STATE,ONLY: DOMAIN_INTERNAL_STATE      &
                                            ,WRAP_DOMAIN_INTERNAL_STATE 
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: COMM_MY_DOMAIN                   &  !<-- MPI communicator for each individual domain
                                      ,COMM_TO_MY_PARENT                &  !<-- Current domain's MPI communicator to its parent
                                      ,MAX_DOMAINS                      &  !<-- Maximum # of domains 
                                      ,MY_DOMAIN_ID                     &  !<-- ID of current domain
                                      ,NUM_CHILDREN                     &  !<-- Current domain's number of children
                                      ,NUM_DOMAINS                         !<-- Total number of domains
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(IN) :: CHILD_ID             &  !<-- Domain IDs of current domain's children
                                                           ,COMM_TO_MY_CHILDREN  &  !<-- Current domain's MPI communicators to its children
                                                           ,FTASKS_DOMAIN        &  !<-- # of forecast tasks on each domain
                                                           ,ID_PARENTS_IN           !<-- IDs of parents of nested domains
!
      INTEGER(kind=KINT),DIMENSION(MAX_DOMAINS),INTENT(IN) :: DOMAIN_ID_TO_RANK  !<-- Configure file associated with each domain ID
!
      REAL(kind=KFPT),DIMENSION(1:NUM_DOMAINS),INTENT(IN) :: DT            !<-- Timesteps for all domains (DOMAIN Components)
!
      TYPE(ESMF_GridComp) :: DOMAIN_GRID_COMP                              !<-- The current DOMAIN component
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_DOMAIN                   !<-- Export state of the current DOMAIN Component
!
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_CPL_NEST              &  !<-- Parent-Child Coupler import state
                                       ,EXP_STATE_CPL_NEST                 !<-- Parent-Child Coupler export state
!
      TYPE(ESMF_CplComp),INTENT(OUT) :: PARENT_CHILD_COUPLER_COMP          !<-- Parent-Child Coupler Component
!
!---------------------
!***  Local Variables
!---------------------
!
!!    TYPE CHILD_TASKS
!!      INTEGER,DIMENSION(:,:),POINTER :: LIMITS
!!    END TYPE CHILD_TASKS
!
!!    TYPE HANDLES_MULTI     
!!      INTEGER,DIMENSION(:),POINTER :: RECV
!!    END TYPE HANDLES_MULTI
!
      INTEGER(kind=KINT) :: ITS,ITE,JTS,JTE,LM                          &
                           ,IDS,IDE,JDS,JDE                             &
                           ,ID_DOM,IERR                                 &
                           ,IHANDLE_RECV,IHANDLE_SEND                   &
                           ,INDX_CW,INDX_Q                              &
                           ,KOUNT,LMP1,MYPE                             &
                           ,N,N_BLEND_H,N_BLEND_V                       &
                           ,NKOUNT,NN,NX
!
      INTEGER(kind=KINT) :: RC,RC_NESTSET
!
      INTEGER(kind=KINT),DIMENSION(4) :: LIMITS
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: ISTAT,JSTAT
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: PARENT_CHILD_RATIO
!
!!!   INTEGER(kind=KINT),DIMENSION(:,:),ALLOCATABLE :: ISTAT_ALL
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: ARRAY_1D
!
!!!!  TYPE(HANDLES_RECV),DIMENSION(:),POINTER :: HANDLES
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      TYPE(WRAP_DOMAIN_INTERNAL_STATE) :: WRAP_DOMAIN
!
      TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE
!
      integer :: itemcount,nsize
      character(len=14),dimension(1:200) :: itemnamelist
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC        =ESMF_SUCCESS
      RC_NESTSET=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Create import/export states for the parent-child coupler.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Nesting"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_520r
      IMP_STATE_CPL_NEST=ESMF_StateCreate(Name     ="Nesting Coupler Import" &  !<-- The Nesting Coupler import state name
                                         ,statetype=ESMF_STATE_IMPORT        &
                                         ,rc       =RC)
!
      EXP_STATE_CPL_NEST=ESMF_StateCreate(Name     ="Nesting Coupler Export" &  !<-- The Nesting Coupler export state name
                                         ,statetype=ESMF_STATE_EXPORT        &
                                         ,rc       =RC)
#else
      IMP_STATE_CPL_NEST=ESMF_StateCreate(stateName="Nesting Coupler Import" &  !<-- The Nesting Coupler import state name
                                         ,statetype=ESMF_STATE_IMPORT        &
                                         ,rc       =RC)
!
      EXP_STATE_CPL_NEST=ESMF_StateCreate(stateName="Nesting Coupler Export" &  !<-- The Parent-Child Cpl export state name
                                         ,statetype=ESMF_STATE_EXPORT        &
                                         ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create the Parent-Child Coupler.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Parent-Child Coupler Component"
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
!***  Register the coupler's Init, Run, and Finalize steps.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the Parent-Child Coupler's Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_CplCompSetServices(comp          =PARENT_CHILD_COUPLER_COMP &  ! <-- The Parent-Child Coupler component
                                  ,subroutineName=PARENT_CHILD_CPL_REGISTER &  ! <-- The user's subroutineName
                                  ,rc            =RC)
#else
      CALL ESMF_CplCompSetServices(cplcomp       =PARENT_CHILD_COUPLER_COMP &  ! <-- The Nesting coupler component
                                  ,userRoutine   =PARENT_CHILD_CPL_REGISTER &  ! <-- The user's subroutineName
                                  ,rc            =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  From the DOMAIN export state find out if this task is a forecast
!***  task and if it is on a parent domain.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Parent-Child CPL Setup: Extract Fcst-or-Write Flag from DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='Fcst-or-Write Flag'                 &  !<-- Name of the attribute to extract
                            ,value=I_AM_A_FCST_TASK                     &  !<-- Am I a forecast task?
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Parent/Not-a-Parent Flag from DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='I-Am-A-Parent Flag'                 &  !<-- Name of the attribute to extract
                            ,value=I_AM_A_PARENT                        &  !<-- Am I on a nested domain?
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  We want the write tasks to see the Parent-Child Coupler but they
!***  should never actually use it therefore they may return now.
!-----------------------------------------------------------------------
!
#ifdef ESMF_3
!!!   IF(I_AM_A_FCST_TASK==ESMF_FALSE.OR.I_AM_A_PARENT==ESMF_FALSE)RETURN
      IF(I_AM_A_FCST_TASK==ESMF_FALSE)RETURN
#else
!!!   IF(.NOT.I_AM_A_FCST_TASK.OR..NOT.I_AM_A_PARENT)RETURN
      IF(.NOT.I_AM_A_FCST_TASK)RETURN
#endif
!
!-----------------------------------------------------------------------
!***  Load key variables into the coupler's import state.
!-----------------------------------------------------------------------
!
!-------------------------------
!***  Maximum number of domains
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Max # of Domains to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='MAX_DOMAINS'                        &  !<-- Maximum # of domains
                            ,value=MAX_DOMAINS                          &  !<-- Insert this into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------
!***  Current Domain's ID 
!-------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Current Domain ID to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='MY_DOMAIN_ID'                       &  !<-- Current domain's ID
                            ,value=MY_DOMAIN_ID                         &  !<-- Insert this into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------
!***  Total Number of Domains
!-----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Number of Domains to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='NUM_DOMAINS'                        &  !<-- Total number of domains
                            ,value=NUM_DOMAINS                          &  !<-- Insert this into the import state
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------------------------
!***  The association of domains and their configure files
!----------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Domain/ConfigFile Association to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='DOMAIN_ID_TO_RANK'              &  !<-- The association of domains and their config files
                            ,count    =MAX_DOMAINS                      &  !<-- Maximum # of domains
                            ,valueList=DOMAIN_ID_TO_RANK                &  !<-- Insert this into the import state
                            ,rc       =RC)
#else
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='DOMAIN_ID_TO_RANK'              &  !<-- The association of domains and their config files
                            ,itemCount=MAX_DOMAINS                      &  !<-- Maximum # of domains
                            ,valueList=DOMAIN_ID_TO_RANK                &  !<-- Insert this into the import state
                            ,rc       =RC)
#endif
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
      MESSAGE_CHECK="Add Number of Fcst Tasks Per Domain to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='FTASKS_DOMAIN'                  &  !<-- Number of forecast tasks on each domain
                            ,count    =NUM_DOMAINS                      &  !<-- Number of domains
                            ,valueList=FTASKS_DOMAIN                    &  !<-- Insert this into the import state
                            ,rc       =RC)
#else
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='FTASKS_DOMAIN'                  &  !<-- Number of forecast tasks on each domain
                            ,itemCount=NUM_DOMAINS                      &  !<-- Number of domains
                            ,valueList=FTASKS_DOMAIN                    &  !<-- Insert this into the import state
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------------
!***  Fundamental Timestep on Domains
!-------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Number of Fcst Tasks Per Domain to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='DOMAIN_DTs'                     &  !<-- Number of forecast tasks on each domain
                            ,count    =NUM_DOMAINS                      &  !<-- Number of domains
                            ,valueList=DT                               &  !<-- Insert this into the import state
                            ,rc       =RC)
#else
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='DOMAIN_DTs'                     &  !<-- Number of forecast tasks on each domain
                            ,itemCount=NUM_DOMAINS                      &  !<-- Number of domains
                            ,valueList=DT                               &  !<-- Insert this into the import state
                            ,rc       =RC)
#endif
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
      MESSAGE_CHECK="Add Domain IDs of Parents to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='ID_PARENTS'                     &  !<-- IDs of parent domain
                            ,count    =NUM_DOMAINS                      &  !<-- Number of domains
                            ,valueList=ID_PARENTS_IN                    &  !<-- Insert this into the import state
                            ,rc       =RC)
#else
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='ID_PARENTS'                     &  !<-- IDs of parent domain
                            ,itemCount=NUM_DOMAINS                      &  !<-- Number of domains
                            ,valueList=ID_PARENTS_IN                    &  !<-- Insert this into the import state
                            ,rc       =RC)
#endif
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
      MESSAGE_CHECK="Add Number of Children to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='NUM_CHILDREN'                       &  !<-- This DOMAIN component's # of children
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
      IF(NUM_CHILDREN>0)THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Add Parent-to-Child Communicators to the Parent-Child CPL Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
        CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST             &  !<-- The Nesting Coupler's import state
                              ,name     ='Parent-to-Child Comms'        &  !<-- Name of Attribute
                              ,count    =NUM_CHILDREN                   &  !<-- Length of inserted array
                              ,valueList=COMM_TO_MY_CHILDREN            &  !<-- Communicators to my children
                              ,rc       =RC)
#else
        CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST             &  !<-- The Nesting Coupler's import state
                              ,name     ='Parent-to-Child Comms'        &  !<-- Name of Attribute
                              ,itemCount=NUM_CHILDREN                   &  !<-- Length of inserted array
                              ,valueList=COMM_TO_MY_CHILDREN            &  !<-- Communicators to my children
                              ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!----------------------------
!***  Communicator to Parent
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Child-to-Parent Communicator to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
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
      MESSAGE_CHECK="Add Domain Intercommunicator to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='Single Domain Comm'                 &  !<-- Name of Attribute
                            ,value=COMM_MY_DOMAIN                       &  !<-- The communicator to my parent
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------
!***  Number of fcst tasks on this domain
!-----------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add NUM_PES_FCST to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='NUM_PES_FCST'                       &  !<-- The name of the Attribute 
                            ,value=NUM_PES_FCST                         &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='NUM_PES_FCST'                       &  !<-- Name of Attribute
                            ,value=NUM_PES_FCST                         &  !<-- The # of fcst tasks on this domain
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------------
!***  Subdomain integration limits on my task
!---------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Integration Limits from DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='ITS'                                &  !<-- The name of the Attribute 
                            ,value=ITS                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='ITE'                                &  !<-- The name of the Attribute
                            ,value=ITE                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JTS'                                &  !<-- The name of the Attribute
                            ,value=JTS                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JTE'                                &  !<-- The name of the Attribute
                            ,value=JTE                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='LM'                                 &  !<-- The name of the Attribute
                            ,value=LM                                   &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='NHALO'                              &  !<-- The name of the Attribute
                            ,value=NHALO                                &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Subdomain Integration Limits to the P-C Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='ITS'                                &  !<-- The name of the Attribute
                            ,value=ITS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='ITE'                                &  !<-- The name of the Attribute
                            ,value=ITE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='JTS'                                &  !<-- The name of the Attribute
                            ,value=JTS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='JTE'                                &  !<-- The name of the Attribute
                            ,value=JTE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='LM'                                 &  !<-- The name of the Attribute
                            ,value=LM                                   &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='NHALO'                              &  !<-- The name of the Attribute
                            ,value=NHALO                                &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------------
!***  Subdomain integration limits for all fcst tasks on my domain 
!------------------------------------------------------------------
!
!
      CALL ESMF_GridCompGetInternalState(DOMAIN_GRID_COMP               &  !<-- The DOMAIN component
                                        ,WRAP_DOMAIN                    &  !<-- Extract the pointer to my DOMAIN internal state
                                        ,RC )
!
      DOMAIN_INT_STATE=>wrap_domain%DOMAIN_INT_STATE                       !<-- Point at my DOMAIN internal state
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add All Fcst Task Integration Limits to the P-C Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='LOCAL ISTART'                   &  !<-- Name of Attribute
                            ,count    =FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- Length of inserted array (# of fcst tasks on domain)
                            ,valueList=domain_int_state%LOCAL_ISTART    &  !<-- Starting I's on my domain's fcst tasks
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='LOCAL IEND'                     &  !<-- Name of Attribute
                            ,count    =FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- Length of inserted array (# of fcst tasks on domain)
                            ,valueList=domain_int_state%LOCAL_IEND      &  !<-- Ending I's on my domain's fcst tasks
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='LOCAL JSTART'                   &  !<-- Name of Attribute
                            ,count    =FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- Length of inserted array (# of fcst tasks on domain)
                            ,valueList=domain_int_state%LOCAL_JSTART    &  !<-- Starting J's on my domain's fcst tasks
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='LOCAL JEND'                     &  !<-- Name of Attribute
                            ,count    =FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- Length of inserted array (# of fcst tasks on domain)
                            ,valueList=domain_int_state%LOCAL_JEND      &  !<-- Ending J's on my domain's fcst tasks
                            ,rc       =RC)
#else
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='LOCAL ISTART'                   &  !<-- Name of Attribute
                            ,itemCount=FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- Length of inserted array (# of fcst tasks on domain)
                            ,valueList=domain_int_state%LOCAL_ISTART    &  !<-- Starting I's on my domain's fcst tasks
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='LOCAL IEND'                     &  !<-- Name of Attribute
                            ,itemCount=FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- Length of inserted array (# of fcst tasks on domain)
                            ,valueList=domain_int_state%LOCAL_IEND      &  !<-- Ending I's on my domain's fcst tasks
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='LOCAL JSTART'                   &  !<-- Name of Attribute
                            ,itemCount=FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- Length of inserted array (# of fcst tasks on domain)
                            ,valueList=domain_int_state%LOCAL_JSTART    &  !<-- Starting J's on my domain's fcst tasks
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='LOCAL JEND'                     &  !<-- Name of Attribute
                            ,itemCount=FTASKS_DOMAIN(MY_DOMAIN_ID)      &  !<-- Length of inserted array (# of fcst tasks on domain)
                            ,valueList=domain_int_state%LOCAL_JEND      &  !<-- Ending J's on my domain's fcst tasks
                            ,rc       =RC)
#endif
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
      MESSAGE_CHECK="Extract Full Domain Dimensions from DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='IDS'                                &  !<-- The name of the Attribute 
                            ,value=IDS                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='IDE'                                &  !<-- The name of the Attribute 
                            ,value=IDE                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JDS'                                &  !<-- The name of the Attribute 
                            ,value=JDS                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JDE'                                &  !<-- The name of the Attribute 
                            ,value=JDE                                  &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Full Domain Dimensions to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='IDS'                                &  !<-- The name of the Attribute
                            ,value=IDS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='IDE'                                &  !<-- The name of the Attribute
                            ,value=IDE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='JDS'                                &  !<-- The name of the Attribute
                            ,value=JDS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
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
      MESSAGE_CHECK="Extract Boundary Blending Region Widths from DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='LNSH'                               &  !<-- The name of the Attribute 
                            ,value=N_BLEND_H                            &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='LNSV'                               &  !<-- The name of the Attribute
                            ,value=N_BLEND_V                            &  !<-- The Attribute to be retrieved
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Add Boundary Blending Region Widths to the Parent-Child Cpl Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='N_BLEND_H'                          &  !<-- The name of the Attribute
                            ,value=N_BLEND_H                            &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
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
!     call esmf_stateget(state   =exp_state_domain                      &  !<-- The Parent's DOMAIN export state
!                       ,itemcount=itemcount                            &  !<-- Extract FIS Field
!                       ,itemNamelist=itemnamelist                         &  !<-- Extract FIS Field
!                       ,rc      =RC)
!     if(mype==0)then
!       write(0,*)' X itemcount=',itemcount,' rc=',rc
!       write(0,*)' itemnamelist=',itemnamelist(1:itemcount)
!     endif
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract FIS from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='FIS'                                 &  !<-- Extract FIS Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
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
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------
!***  Transfer geographic latitude
!----------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract GLAT from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='GLAT'                                &  !<-- Extract GLAT Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert GLAT into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------
!***  Transfer geographic longitude
!-----------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract GLON from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='GLON'                                &  !<-- Extract GLON Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert GLON into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------------------------------
!***  Transfer PT,PDTOP,PSGML1,SG1,SG2,SGML2,DSG2,PDSG1
!-------------------------------------------------------
!
      LMP1=LM+1
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PT from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The Parent's DOMAIN export state
                            ,name ='PT'                                 &  !<-- Extract PT
                            ,value=PT                                   &  !<-- Put the extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='PT'                                 &  !<-- Insert PT
                            ,value=PT                                   &  !<-- Insert this value
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The Parent's DOMAIN export state
                            ,name ='PDTOP'                              &  !<-- Extract PDTOP
                            ,value=PDTOP                                &  !<-- Put the extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='PDTOP'                              &  !<-- Insert PDTOP
                            ,value=PDTOP                                &  !<-- Insert this value
                            ,rc   =RC)
!
      ALLOCATE(ARRAY_1D(1:LM))
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='PSGML1'                         &  !<-- Extract PGMSL1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='PSGML1'                         &  !<-- Insert PGMSL1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='SGML2'                          &  !<-- Extract SGML2
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='SGML2'                          &  !<-- Insert SGML2
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='DSG2'                           &  !<-- Extract DSG2
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='DSG2'                           &  !<-- Insert DSG2
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='PDSG1'                          &  !<-- Extract PDSG1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='PDSG1'                          &  !<-- Insert PDSG1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      DEALLOCATE(ARRAY_1D)
!
      ALLOCATE(ARRAY_1D(1:LMP1))
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='SG1'                            &  !<-- Extract SG1
                            ,count    =LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='SG1'                            &  !<-- Insert SG1
                            ,count    =LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='SG2'                            &  !<-- Extract SG2
                            ,count    =LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='SG2'                            &  !<-- Insert SG2
                            ,count    =LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
#else
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='PSGML1'                         &  !<-- Extract PGMSL1
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here 
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='PSGML1'                         &  !<-- Insert PGMSL1
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's ATM export state
                            ,name     ='SGML2'                          &  !<-- Extract SGML2
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='SGML2'                          &  !<-- Insert SGML2
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='DSG2'                           &  !<-- Extract DSG2
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='DSG2'                           &  !<-- Insert DSG2
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='PDSG1'                          &  !<-- Extract PDSG1
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='PDSG1'                          &  !<-- Insert PDSG1
                            ,itemCount=LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      DEALLOCATE(ARRAY_1D)
!
      ALLOCATE(ARRAY_1D(1:LMP1))
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's ATM export state
                            ,name     ='SG1'                            &  !<-- Extract SG1
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='SG1'                            &  !<-- Insert SG1
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's ATM export state
                            ,name     ='SG2'                            &  !<-- Extract SG2
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='SG2'                            &  !<-- Insert SG2
                            ,itemCount=LMP1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DEALLOCATE(ARRAY_1D)
!
!------------------------
!***  Transfer DY and DX
!------------------------
!
      NKOUNT=JDE-JDS+1
      ALLOCATE(ARRAY_1D(1:NKOUNT))
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract DY from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The Parent's DOMAIN export state
                            ,name ='DYH'                                &  !<-- Extract DYH
                            ,value=DYH                                  &  !<-- Put the extracted value here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='DYH'                                &  !<-- Insert DYH
                            ,value=DYH                                  &  !<-- Insert this value
                            ,rc   =RC)
!
#ifdef ESMF_3
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='DXH'                            &  !<-- Extract DXH
                            ,count    =NKOUNT                           &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Parent-Child Coupler's import state
                            ,name     ='DXH'                            &  !<-- Insert DXH
                            ,count    =NKOUNT                           &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
!
#else
      CALL ESMF_AttributeGet(state    =EXP_STATE_DOMAIN                 &  !<-- The Parent's DOMAIN export state
                            ,name     ='DXH'                            &  !<-- Extract DXH
                            ,itemCount=NKOUNT                           &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values here
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST               &  !<-- The Nesting Coupler's import state
                            ,name     ='DXH'                            &  !<-- Insert DXH
                            ,itemCount=NKOUNT                           &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Insert these values
                            ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DEALLOCATE(ARRAY_1D)
!
!-----------------------------------------------------------------------
!
#ifdef ESMF_3
      parents_only: IF(I_AM_A_PARENT==ESMF_TRUE)THEN             
#else
      parents_only: IF(I_AM_A_PARENT)THEN             
#endif
!
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
        MESSAGE_CHECK="Add Parent-Child DT Ratios to the Parent-Child Cpl Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
        CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST             &  !<-- The Parent-Child Coupler's import state
                              ,name     ='Parent-Child Time Ratio'      &  !<-- Name of Attribute
                              ,count    =NUM_CHILDREN                   &  !<-- Length of inserted array
                              ,valueList=PARENT_CHILD_RATIO             &  !<-- The communicator to my parent
                              ,rc       =RC)
#else
        CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST             &  !<-- The Nesting Coupler's import state
                              ,name     ='Parent-Child Time Ratio'      &  !<-- Name of Attribute
                              ,itemCount=NUM_CHILDREN                   &  !<-- Length of inserted array
                              ,valueList=PARENT_CHILD_RATIO             &  !<-- The communicator to my parent
                              ,rc       =RC)
#endif
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
        MESSAGE_CHECK="Add Domain IDs of Children to Parent-Child Cpl Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
        CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST             &  !<-- The Parent-Child Coupler's import state
                              ,name     ='CHILD_IDs'                    &  !<-- Name of Attribute
                              ,count    =NUM_CHILDREN                   &  !<-- Length of inserted array
                              ,valueList=CHILD_ID                       &  !<-- The children's IDs of this DOMAIN Component
                              ,rc       =RC)
#else
        CALL ESMF_AttributeSet(state    =IMP_STATE_CPL_NEST             &  !<-- The Nesting Coupler's import state
                              ,name     ='CHILD_IDs'                    &  !<-- Name of Attribute
                              ,itemCount=NUM_CHILDREN                   &  !<-- Length of inserted array
                              ,valueList=CHILD_ID                       &  !<-- The children's IDs of this ATM Component
                              ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF parents_only
!
!-----------------------------------------------------------------------
!***  Now transfer the parent's prognostic arrays from the DOMAIN export
!***  state to the Parent-Child coupler import state that will be
!***  required for the children's boundary data.
!-----------------------------------------------------------------------
!
!-----------------
!***  Transfer PD
!-----------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PD Field from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='PD'                                  &  !<-- Extract PD Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert PD into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------
!***  Transfer Layer Interface Pressures
!----------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PINT Field from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='PINT'                                &  !<-- Extract PINT Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert PINT into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
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
      MESSAGE_CHECK="Extract T Field from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='T'                                   &  !<-- Extract T Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert T into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
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
      MESSAGE_CHECK="Extract U Wind Field from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='U'                                   &  !<-- Extract U Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert U into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
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
      MESSAGE_CHECK="Extract V Wind Field from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='V'                                   &  !<-- Extract V Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert V into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------
!***  Transfer TRACERS Field
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Tracers Field from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='TRACERS'                             &  !<-- Extract TRACERS Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Tracers into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------
!***  Transfer Sea Mask Field
!-----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract SM Field from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DOMAIN                      &  !<-- The Parent's DOMAIN export state
                        ,itemName='SM'                                  &  !<-- Extract Seas Mask Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert SM into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=IMP_STATE_CPL_NEST                       &  !<-- The Parent-Child Coupler's import state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
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
      MESSAGE_CHECK="Extract INDX_Q from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The Parent's DOMAIN export state
                            ,name ='INDX_Q'                             &  !<-- Name of Attribute to extract
                            ,value=INDX_Q                               &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert INDX_Q into Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='INDX_Q'                             &  !<-- The name of the Attribute to insert
                            ,value=INDX_Q                               &  !<-- The Attribute to be inserted
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
      MESSAGE_CHECK="Extract INDX_CW from Parent DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DOMAIN                     &  !<-- The Parent's DOMAIN export state
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
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=IMP_STATE_CPL_NEST                   &  !<-- The Parent-Child Coupler's import state
                            ,name ='INDX_CW'                            &  !<-- The name of the Attribute to insert
                            ,value=INDX_CW                              &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  All parent tasks need to know the local subdomain limits of each
!***  task on their children.  
!-----------------------------------------------------------------------
!
      LIMITS(1)=ITS
      LIMITS(2)=ITE
      LIMITS(3)=JTS
      LIMITS(4)=JTE
!
!-----------------------------------------------------------------------
!***  Child tasks send their subdomain limits to parent task 0.
!-----------------------------------------------------------------------
!
      IF(COMM_TO_MY_PARENT/=-999)THEN                                      !<-- Select child tasks
        CALL MPI_COMM_RANK(COMM_TO_MY_PARENT,MYPE,IERR)                    !<-- Obtain my rank
        CALL MPI_SEND(LIMITS                                            &  !<-- Child task sends its subdomain limits
                     ,4                                                 &  !<-- # of indices sent
                     ,MPI_INTEGER                                       &  !<-- Indices are integers
                     ,0                                                 &  !<-- Indices sent to parent fcst task 0
                     ,MYPE                                              &  !<-- This task's rank
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator between parent and child
                     ,IERR)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Rank 0 parent tasks recv their children's tasks' subdomain limits
!***  then send the children the integration limits of every fcst task
!***  on the parent domain.
!-----------------------------------------------------------------------
!
!!!   IF(NUM_CHILDREN>0)THEN
#ifdef ESMF_3
      IF(I_AM_A_PARENT==ESMF_TRUE)THEN
#else
      IF(I_AM_A_PARENT)THEN
#endif
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
!
            DO NN=1,FTASKS_DOMAIN(ID_DOM)
              CALL MPI_RECV(CTASK_LIMITS(N)%LIMITS(1:4,NN)              &  !<-- Subdomain limits of fcst task NN on child N
                           ,4                                           &  !<-- # of index limits received
                           ,MPI_INTEGER                                 &  !<-- Indices are integers
                           ,NN-1                                        &  !<-- Local index of child fcst task NN (the sender)
                           ,NN-1                                        &  !<-- MPI tag
                           ,COMM_TO_MY_CHILDREN(N)                      &  !<-- MPI communicator between parent and child N
                           ,JSTAT                                       &
                           ,IERR)
            ENDDO
!
          ENDIF
!
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  For each child, parent task 0 sends all of the other parent tasks
!***  the local subdomain limits of each child task so the parent tasks
!***  will know how to divvy up the child domain data.
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
!***  If there are moving nests then their parents need the Bundles
!***  of 2-D and 3-D variables in the Dynamics and Physics internal
!***  states that need to be updated after the nests move.  The
!***  Bundles are unloaded from the DOMAIN export state and loaded
!***  into the Parent-Child coupler import state.  If there are no
!***  moving nests then the Bundles are empty.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Bundles for Moving Nests from Domain Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state      =EXP_STATE_DOMAIN                   &  !<-- The Parent's DOMAIN export state
                        ,itemName   ='Move_Bundle H'                    &  !<-- Name of Bundle of internal state H arrays
                        ,fieldbundle=MOVE_BUNDLE_H                      &  !<-- Put the extracted Bundle here
                        ,rc         =RC)
!
      CALL ESMF_StateGet(state      =EXP_STATE_DOMAIN                   &  !<-- The Parent's DOMAIN export state
                        ,itemName   ='Move_Bundle V'                    &  !<-- Name of Bundle of internal state V arrays
                        ,fieldbundle=MOVE_BUNDLE_V                      &  !<-- Put the extracted Bundle here
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Bundle for Moving Nests into P-C Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state      =IMP_STATE_CPL_NEST                 &  !<-- The Parent-Child Coupler's import state
                        ,fieldbundle=MOVE_BUNDLE_H                      &  !<-- The Bundle of internal state H arrays to update
                        ,rc         =RC)
!
      CALL ESMF_StateAdd(state      =IMP_STATE_CPL_NEST                 &  !<-- The Parent-Child Coupler's import state
                        ,fieldbundle=MOVE_BUNDLE_V                      &  !<-- The Bundle of internal state V arrays to update
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_NESTSET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
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
                                          ,N_BLEND_H_CHILD              &
                                          ,N_BLEND_V_CHILD              &
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
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: ITS,ITE,JTS,JTE                  &
                                      ,IDS,IDE,JDS,JDE                  &
                                      ,NUM_CHILDREN
!
      INTEGER(kind=KINT),DIMENSION(1:NUM_CHILDREN),INTENT(IN) :: IM_CHILD         &
                                                                ,JM_CHILD         &
                                                                ,N_BLEND_H_CHILD  &
                                                                ,N_BLEND_V_CHILD 
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(IN) :: FTASKS_DOMAIN  &
                                                           ,MY_CHILDREN_ID         
!
      TYPE(CHILD_TASKS),DIMENSION(:),POINTER,INTENT(IN) :: CTASK_LIMITS
!
      TYPE(ESMF_Config),DIMENSION(1:NUM_CHILDREN),INTENT(INOUT) :: CF
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: N,N_CHILD_TASKS,NUM_CHILD_TASKS             &
                           ,THIS_CHILD_ID
!
      INTEGER(kind=KINT) :: EAST_LIMIT1 ,EAST_LIMIT2                    &
                           ,WEST_LIMIT1 ,WEST_LIMIT2                    &
                           ,NORTH_LIMIT1,NORTH_LIMIT2                   &
                           ,SOUTH_LIMIT1,SOUTH_LIMIT2
!
      INTEGER(kind=KINT) :: RC,RC_SET 
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_SET=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Allocate the pointers that hold the four H and V parent points
!***  that surround each child point in the child's boundary region.
!-----------------------------------------------------------------------
!
      ALLOCATE(PARENT_4_INDICES_H(1:NUM_CHILDREN))
      ALLOCATE(PARENT_4_INDICES_V(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  Allocate the pointers that hold the weights of the four H and V 
!***  parent points that surround each child point in the child's 
!***  boundary region.
!-----------------------------------------------------------------------
!
      ALLOCATE(PARENT_4_WEIGHTS_H(1:NUM_CHILDREN))
      ALLOCATE(PARENT_4_WEIGHTS_V(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  Allocate the arrays that hold the number of child tasks 
!***  on each side of the child boundaries that will be sent
!***  data from the parent tasks.
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
!***  Allocate the pointers that will hold the ranks of all child tasks 
!***  on each side of the child boundaries that will be sent data
!***  from the parent tasks.
!-----------------------------------------------------------------------
!
      ALLOCATE(CHILDTASK_BNDRY_H_RANKS(1:NUM_CHILDREN))
      ALLOCATE(CHILDTASK_BNDRY_V_RANKS(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  Allocate the pointers for starting/ending I's and J's on each
!***  parent task for each side of the boundary.
!-----------------------------------------------------------------------
!
      ALLOCATE(CHILDTASK_H_SAVE(1:NUM_CHILDREN))
      ALLOCATE(CHILDTASK_V_SAVE(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  Allocate the arrays that hold the parent's grid location
!***  of each child's southwest corner.
!-----------------------------------------------------------------------
!
      ALLOCATE(I_PARENT_SW(1:NUM_CHILDREN)                              &
              ,J_PARENT_SW(1:NUM_CHILDREN))
!
!-----------------------------------------------------------------------
!***  Extract relevant information from the children's configure files.
!-----------------------------------------------------------------------
!
      child_loop_0: DO N=1,NUM_CHILDREN
!
        THIS_CHILD_ID=MY_CHILDREN_ID(N)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent-Child_Interp_Setup: Extract I of SW Mass Point on Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(N)                       &  !<-- The child's config object
                                    ,value =I_PARENT_SW(N)              &  !<-- The variable filled (parent I of child's SW corner)
                                    ,label ='i_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)  
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent-Child_Interp_Setup: Extract J of SW Mass Point on Parent"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(N)                       &  !<-- The child's config object
                                    ,value =J_PARENT_SW(N)              &  !<-- The variable filled (parent J of child's SW corner)
                                    ,label ='j_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Invert the Parent-to-Child space ratio for computation.
!-----------------------------------------------------------------------
!
        CHILD_PARENT_SPACE_RATIO(N)=1./REAL(PARENT_CHILD_SPACE_RATIO(N))
!
!-----------------------------------------------------------------------
!***  Allocate the individual pointers holding the four H points of
!***  the parent that surround this child's boundary region H points
!***  and the bilinear interpolation weights of the four parent points
!***  surrounding those same child points.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!   ***************************  NOTE  *****************************
!-----------------------------------------------------------------------
!     Although the H points in the nests' boundary region cover only
!     N_BLEND rows, we actually need to have the nests' PD values
!     one row further.  That is because we also need PD values at the 
!     V points in the nests' boundary region to perform the proper
!     hydrostatic updating of the winds by the parents there.  To
!     do the 4-point average needed to obtain PD on V points, we
!     necessarily must have them on mass points one row beyond where 
!     they are needed for the mass points alone.
!-----------------------------------------------------------------------
!   ***************************  NOTE  *****************************
!-----------------------------------------------------------------------
!
        SOUTH_LIMIT1=1
        SOUTH_LIMIT2=N_BLEND_H_CHILD(N)+1                                  !<-- Extend the region by 1 row for 4-point averaging of PD
!
        NORTH_LIMIT1=JM_CHILD(N)-N_BLEND_H_CHILD(N)                        !<-- Extend the region by 1 row for 4-point averaging of PD
        NORTH_LIMIT2=JM_CHILD(N)
!
        WEST_LIMIT1=1
        WEST_LIMIT2=N_BLEND_H_CHILD(N)+1                                   !<-- Extend the region by 1 row for 4-point averaging of PD
!
        EAST_LIMIT1=IM_CHILD(N)-N_BLEND_H_CHILD(N)                         !<-- Extend the region by 1 row for 4-point averaging of PD
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
!***  Allocate the individual pointers holding the four V points of
!***  the parent that surround this child's boundary region V points.
!-----------------------------------------------------------------------
!
        SOUTH_LIMIT1=1
        SOUTH_LIMIT2=N_BLEND_V_CHILD(N)
!
        NORTH_LIMIT1=JM_CHILD(N)-1-N_BLEND_V_CHILD(N)+1
        NORTH_LIMIT2=JM_CHILD(N)-1
!
        WEST_LIMIT1=1
        WEST_LIMIT2=N_BLEND_V_CHILD(N)
!
        EAST_LIMIT1=IM_CHILD(N)-1-N_BLEND_V_CHILD(N)+1
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
!***  What is the number of forecast tasks on the child domain?
!-----------------------------------------------------------------------
!
        NUM_CHILD_TASKS=FTASKS_DOMAIN(THIS_CHILD_ID)
!
!-----------------------------------------------------------------------
!***  Allocate the pointers for starting/ending I's and J's on each
!***  parent task for each side of the boundary.
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
!***  Allocate the pointers for the child task ID's that contain
!***  segments of the child boundary within a parent task for
!***  each side of the boundary.
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
      SUBROUTINE PREPARE_NEST_INTERP_FACTORS(NUM_CHILDREN)
!     SUBROUTINE PREPARE_NEST_INTERP_FACTORS(NUM_CHILDREN               &
!                                           ,NUM_DOMAINS                &
!                                           ,MY_CHILDREN_ID             &
!                                           ,FTASKS_DOMAIN              &
!                                           ,I_PARENT_SW                &
!                                           ,J_PARENT_SW                &
!                                           ,ITS,ITE,JTS,JTE            &
!                                           ,IDS,IDE,JDS,JDE )
!
!-----------------------------------------------------------------------
!***  Call the routine that computes the interpolation factors
!***  each parent needs in order to interpolate its data to 
!***  its nests' boundaries.
!
!***  Only parent tasks execute this routine.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: NUM_CHILDREN                        !<-- Number of nests on the parent
!     INTEGER(kind=KINT),INTENT(IN) :: NUM_CHILDREN                     &  !<-- Number of nests on the parent
!                                     ,NUM_DOMAINS                      &  !<-- Total # of domains
!                                     ,ITS,ITE,JTS,JTE                  &  !<-- Horizontal integration limits of parent task
!                                     ,IDS,IDE,JDS,JDE                     !<-- Horizontal dimensions of full parent domain
!
!     INTEGER(kind=KINT),DIMENSION(1:NUM_CHILDREN),INTENT(IN) ::        &
!                                                       I_PARENT_SW        !<-- Parent I of nests' SW corner H point
!                                                       J_PARENT_SW        !<-- Parent J of nests' SW corner H point
!                                                       MY_CHILDREN_ID     !<-- Domain IDs of each nest on this parent
!
!     INTEGER(kind=KINT),DIMENSION(1:NUM_DOMAINS),INTENT(IN) ::         &
!                                                       FTASKS_DOMAIN      !<-- # of fcst tasks on each domain
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: N,N_CHILD_TASKS,NUM_CHILD_TASKS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The parent sets up quantities to be used for general
!***  bilinear interpolation from the parent to its children's
!***  boundary regions.  These quantities are:
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
      child_loop: DO N=1,NUM_CHILDREN
!
        NUM_CHILD_TASKS=FTASKS_DOMAIN(MY_CHILDREN_ID(N))                   !<-- # of fcst tasks on this parent's Nth child
!
!------------------------------------------------------------
!***  Compute interpolation indices and weights for H points
!------------------------------------------------------------
!
        CALL PARENT_TO_CHILD_INTERP_FACTORS('H_POINTS'                  &
                                           ,I_PARENT_SW(N)              &
                                           ,J_PARENT_SW(N)              &
                                           ,N_BLEND_H_CHILD(N)          &
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
                                           ,IDS,IDE,JDS,JDE             & !     Input
!                                                                          -----------------
                                           ,NUM_TASKS_SEND_H_S(N)       & !     Output
                                           ,NUM_TASKS_SEND_H_N(N)       & !       |
                                           ,NUM_TASKS_SEND_H_W(N)       & !       |
                                           ,NUM_TASKS_SEND_H_E(N)       & !       v
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
!------------------------------------------------------------
!***  Compute interpolation indices and weights for V points
!------------------------------------------------------------
!
        CALL PARENT_TO_CHILD_INTERP_FACTORS('V_POINTS'                  &
                                           ,I_PARENT_SW(N)              &
                                           ,J_PARENT_SW(N)              &
                                           ,N_BLEND_V_CHILD(N)          &
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
                                           ,IDS,IDE,JDS,JDE             & !     Input
!                                                                        ----------------
                                           ,NUM_TASKS_SEND_V_S(N)       & !     Output
                                           ,NUM_TASKS_SEND_V_N(N)       & !       |
                                           ,NUM_TASKS_SEND_V_W(N)       & !       |
                                           ,NUM_TASKS_SEND_V_E(N)       & !       v
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
!***  For V point variables, the number of points to be transferred
!***  from parents to their children's boundaries is the same as
!***  the number of computation points (no extensions as is needed
!***  for PDB).
!-----------------------------------------------------------------------
!
        NUM_CHILD_TASKS=FTASKS_DOMAIN(MY_CHILDREN_ID(N))
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
      ENDDO child_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PREPARE_NEST_INTERP_FACTORS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_TO_CHILD_INTERP_FACTORS(FLAG_H_OR_V             &
                                               ,I_PARENT_SW             &
                                               ,J_PARENT_SW             &
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
                                               ,IDS,IDE,JDS,JDE         & !    Input
!                                                                          --------------
                                               ,NUM_TASKS_SEND_S        & !    Output
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
!***  Parent components compute various indices, weights, etc.
!***  needed to generate boundary point data for the given
!***  child throughout the upcoming forecast.
!***  Only parent tasks execute this routine.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_PARENT_SW,J_PARENT_SW          &  !<-- SW corner of nest lies on this I,J of parent
                                      ,IM_CHILD,JM_CHILD                &  !<-- Horizontal dimensions of nest domain
                                      ,N_BLEND                          &  !<-- Width (in rows) of boundary's blending region
                                      ,NUM_CHILD_TASKS                     !<-- # of fcst tasks on the child's domain
! 
      INTEGER(kind=KINT),INTENT(IN) :: ITE,ITS,JTE,JTS                  &  !<-- Index limits on parent task subdomain
                                      ,IDE,IDS,JDE,JDS                     !<-- Full dimensions of parent domain
!
      INTEGER(kind=KINT),DIMENSION(1:4,1:NUM_CHILD_TASKS),INTENT(IN) :: &
                                                                 LIMITS    !<-- ITS,ITE,JTS,JTE on each task of the child
!
      REAL(kind=KFPT),INTENT(IN) :: CHILD_PARENT_SPACE_RATIO               !<-- Ratio of nest grid increment to parent's increment
!
      CHARACTER(*),INTENT(IN) :: FLAG_H_OR_V                               !<-- Are we dealing with H or V child boundary points?
!
      INTEGER(kind=KINT),INTENT(OUT) :: NUM_TASKS_SEND_S                &  !<-- # of child tasks with S bndry segments on this parent task
                                       ,NUM_TASKS_SEND_N                &  !<-- # of child tasks with N bndry segments on this parent task
                                       ,NUM_TASKS_SEND_W                &  !<-- # of child tasks with W bndry segments on this parent task
                                       ,NUM_TASKS_SEND_E                   !<-- # of child tasks with E bndry segments on this parent task
!
      INTEGER(kind=KINT),DIMENSION(NUM_CHILD_TASKS),INTENT(OUT) ::      &
                                                     I_SAVE_LO_SOUTH    &  !<-- Child tasks' westernmost Sbndry I's on this parent task
                                                    ,I_SAVE_HI_SOUTH    &  !<-- Child tasks' easternmost Sbndry I's on this parent task
                                                    ,I_SAVE_LO_NORTH    &  !<-- Child tasks' westernmost Nbndry I's on this parent task
                                                    ,I_SAVE_HI_NORTH    &  !<-- Child tasks' easternmost Nbndry I's on this parent task
                                                    ,J_SAVE_LO_WEST     &  !<-- Child tasks' southernmost Wbndry J's on this parent task
                                                    ,J_SAVE_HI_WEST     &  !<-- Child tasks' northernmost Wbndry J's on this parent task
                                                    ,J_SAVE_LO_EAST     &  !<-- Child tasks' southernmost Ebndry J's on this parent task
                                                    ,J_SAVE_HI_EAST     &  !<-- Child tasks' northernmost Ebndry J's on this parent task
!
                                                   ,LOCAL_TASK_RANK_S   &  !<-- Child task ranks with S bndry on this parent task
                                                   ,LOCAL_TASK_RANK_N   &  !<-- Child task ranks with N bndry on this parent task
                                                   ,LOCAL_TASK_RANK_W   &  !<-- Child task ranks with W bndry on this parent task
                                                   ,LOCAL_TASK_RANK_E      !<-- Child task ranks with E bndry on this parent task
!
      INTEGER(kind=KINT),DIMENSION(:,:,:),POINTER,INTENT(OUT) ::        & 
                                                          I_INDX_SBND   &  !<-- Parent I west/east of child south boundary point 
                                                         ,I_INDX_NBND   &  !<-- Parent I west/east of child north boundary point
                                                         ,I_INDX_WBND   &  !<-- Parent I west/east of child west boundary point
                                                         ,I_INDX_EBND   &  !<-- Parent I west/east of child east boundary point
                                                         ,J_INDX_SBND   &  !<-- Parent J south/north of child south boundary point
                                                         ,J_INDX_NBND   &  !<-- Parent J south/north of child north boundary point
                                                         ,J_INDX_WBND   &  !<-- Parent J south/north of child west boundary point
                                                         ,J_INDX_EBND      !<-- Parent J south/north of child east boundary point
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER,INTENT(OUT) ::           &
                                                         WEIGHTS_SBND   &  !<-- Sbndry bilinear interp wghts for 4 surrounding parent points
                                                        ,WEIGHTS_NBND   &  !<-- Nbndry bilinear interp wghts for 4 surrounding parent points
                                                        ,WEIGHTS_WBND   &  !<-- Wbndry bilinear interp wghts for 4 surrounding parent points
                                                        ,WEIGHTS_EBND      !<-- Ebndry bilinear interp wghts for 4 surrounding parent points
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I_CHILD,IM_END                              &
                           ,J_CHILD,JM_END                              &
                           ,KOUNT_I,KOUNT_J                             &
                           ,N,N_ADD,NC                                  &
                           ,NC_LAST_S,NC_LAST_N,NC_LAST_W,NC_LAST_E     &
                           ,RATIO_P_C
!
      INTEGER(kind=KINT),DIMENSION(1:NUM_CHILD_TASKS) :: I_LIMIT_LO     &
                                                        ,I_LIMIT_HI     &
                                                        ,ITS_CHILD      &
                                                        ,ITE_CHILD      &
                                                        ,J_LIMIT_LO     &
                                                        ,J_LIMIT_HI     &
                                                        ,JTS_CHILD      &
                                                        ,JTE_CHILD      &
                                                        ,NC_HOLD_S      &
                                                        ,NC_HOLD_N      &
                                                        ,NC_HOLD_W      &
                                                        ,NC_HOLD_E
!
      REAL(kind=KFPT) :: ADD_INC,ARG1,ARG2                              &
                        ,R_ITS,R_ITE,R_IEND,R_JTS,R_JTE,R_JEND       
!
      REAL(kind=KFPT) :: PARENT_I_CHILD_EBND,PARENT_I_CHILD_WBND        &
                        ,PARENT_J_CHILD_NBND,PARENT_J_CHILD_SBND        &
                        ,PARENT_S_TASK_LIM_ON_NEST                      &
                        ,PARENT_W_TASK_LIM_ON_NEST                      &
                        ,RATIO_C_P                                      &
                        ,REAL_I_PARENT,REAL_I_START                     &
                        ,REAL_J_PARENT,REAL_J_START                     &
                        ,RECIP_SUM
!
      REAL(kind=KFPT) :: WEIGHT_NE,WEIGHT_NW,WEIGHT_SE,WEIGHT_SW
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RATIO_C_P=CHILD_PARENT_SPACE_RATIO                                   !<-- Child-to-Parent gridspace ratio
      RATIO_P_C=NINT(1./RATIO_C_P)                                         !<-- Parent-to-Child gridspace ratio
!
!-----------------------------------------------------------------------
!***  Create the Real index limits on the parent grid across which
!***  the children's boundary point values will be computed.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***                      !!!!! NOTE !!!!!
!
!***  For the purpose of handling child boundaries, parent tasks will
!***  "BEGIN" directly on their southernmost/westernmost H or V points.
!***  Each parent task covers the gap between itself and the next task
!***  on the parent grid in each direction.
!***  This means that if a child south boundary point lies exactly on 
!***  a parent task point that itself is on the westernmost side of
!***  that parent task's integration subdomain then that child point
!***  will be considered to lie on both that parent task and the
!***  parent task to the west simply because it is the intersection
!***  of the regions managed by both of those parent tasks.
!***  This same notion applies for all other directions and sides.
!-----------------------------------------------------------------------
!
      R_ITE =REAL(ITE)                                                     !<-- REAL Iend of parent task's subdomain
!
      R_JTE =REAL(JTE)                                                     !<-- REAL Jend of parent task's subdomain
!
!-----------------------------------------------------------------------
!***  Because each parent gridpoint covers the gap to the next parent
!***  gridpoint as explained above, increase the search limit for
!***  child boundary points.  That increase would be 1 for both H and V
!***  but due to the nature of the B-Grid layout and the fact that
!***  the I index of child V points on the west boundary and the
!***  J index of the child V points on the south boundary have smaller
!***  grid index values in terms of the parent indices, we must search
!***  for child H points 1/2+0.5*(space_ratio) grid increments further
!***  than for child V points in order to reach the same actual position.
!-----------------------------------------------------------------------
!
!***  In this diagram the H's and V's are points on the parent task's
!***  subdomain while the h's and v's are points on a nest.  It shows
!***  how each parent point must look eastward.  The same goes for
!***  looking northwward.  A parent/nest ratio of 3:1 is used in this
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
! must cover
! area to the
! next H.
! But V with
! the same I 
! as the next H
! is 1.5 farther
! than this H.
! If nest v at
! 1.5 past this
! H is on the east
! bndry of the nest
! then the east h
! on the bndry is
! 1+1/2+1/6.
! That is how
! far we must
! scan from
! this H.
!-----------------------------------------------------------------------
!
      IF(FLAG_H_OR_V=='H_POINTS')THEN
!       ADD_INC=1.5
        ADD_INC=1.5+0.5*RATIO_C_P+EPS

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
!***  What are the child I and J index limits of any sections of its
!***  (the child's) boundary that lie within a parent task's subdomain?
!
!***  What are the indices of the four parent gridpoints surrounding
!***  each child boundary point?
!
!***  What are the bilinear weights associated with each of the four
!***  surrounding parent points to obtain the child boundary point?
!
!***  The parent will use these pieces of information to interpolate
!***  from its grid to its children's boundary points.
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
        PARENT_J_CHILD_SBND=REAL(J_PARENT_SW)                              !<-- J index of parent H for child's south H boundary
        PARENT_J_CHILD_NBND=PARENT_J_CHILD_SBND+(JM_CHILD-1)*RATIO_C_P     !<-- J index of parent H for child's north H boundary
        PARENT_I_CHILD_WBND=REAL(I_PARENT_SW)                              !<-- I index of parent H for child's west H boundary
        PARENT_I_CHILD_EBND=PARENT_I_CHILD_WBND+(IM_CHILD-1)*RATIO_C_P     !<-- I index of parent H for child's east H boundary
        IM_END=IM_CHILD
        JM_END=JM_CHILD
        N_ADD=1                                                            !<-- Blending region along child's boundary
                                                                           !    increased by 1 row to allow 4-pt averaging of PD.
!
      ELSEIF(FLAG_H_OR_V=='V_POINTS')THEN
        PARENT_J_CHILD_SBND=REAL(J_PARENT_SW)-0.5+RATIO_C_P*0.5            !<-- J index of parent V for child's south V boundary
        PARENT_J_CHILD_NBND=PARENT_J_CHILD_SBND+(JM_CHILD-2)*RATIO_C_P     !<-- J index of parent V for child's north V boundary
        PARENT_I_CHILD_WBND=REAL(I_PARENT_SW)-0.5+RATIO_C_P*0.5            !<-- I index of parent V for child's west V boundary
        PARENT_I_CHILD_EBND=PARENT_I_CHILD_WBND+(IM_CHILD-2)*RATIO_C_P     !<-- I index of parent V for child's east V boundary
        IM_END=IM_CHILD-1
        JM_END=JM_CHILD-1
        N_ADD=0                                                            !<-- Blending region along child's boundary
                                                                           !    increased only for mass points (for PD averaging)
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
!
!-----------------------------------------------------------------------
!***  If the northernmost/easternmost extra row of H bndry points
!***  on a nest task coincides with the southern/western boundary 
!***  of a parent task then that parent task will not be associated
!***  with the nest since no bndry V points would be seen by the
!***  parent task.
!-----------------------------------------------------------------------
!
        IF(FLAG_H_OR_V=='H_POINTS')THEN
!
          PARENT_S_TASK_LIM_ON_NEST=REAL(JTS-J_PARENT_SW)*RATIO_P_C+1      !<-- South limit of parent task w/r to nest J
          IF(J_LIMIT_HI(N)-PARENT_S_TASK_LIM_ON_NEST<=0.5)THEN
            J_LIMIT_HI(N)=J_LIMIT_HI(N)-1
          ENDIF
!
          PARENT_W_TASK_LIM_ON_NEST=REAL(ITS-I_PARENT_SW)*RATIO_P_C+1      !<-- West limit of parent task w/r to nest I
          IF(I_LIMIT_HI(N)-PARENT_W_TASK_LIM_ON_NEST<=0.5)THEN
            I_LIMIT_HI(N)=I_LIMIT_HI(N)-1
          ENDIF
        ENDIF
!
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
      NC_LAST_S=-1
      NC_LAST_N=-1
!
      NUM_TASKS_SEND_S=0
      NUM_TASKS_SEND_N=0
!
      IF(FLAG_H_OR_V=='H_POINTS')THEN
        R_ITS=REAL(ITS)-EPS                                                !<-- REAL Istart of parent task's subdomain for H on B grid
      ELSEIF(FLAG_H_OR_V=='V_POINTS')THEN
        R_ITS=REAL(ITS-0.5)-EPS                                            !<-- REAL Istart of parent task's subdomain for V on B grid
      ENDIF
!
      ARG1=REAL(ITE)+ADD_INC
      ARG2=REAL(IDE)
      R_IEND=MIN(ARG1,ARG2)-EPS                                            !<-- REAL Iend of parent task's region for child N/S boundaries
!
!-----------------------------------------------------
!
      REAL_I_START=PARENT_I_CHILD_WBND
!
!-----------------------------------------------------------------------
      i_loop: DO I_CHILD=1,IM_END                                          !<-- Loop over child I's across its South/North boundaries
!-----------------------------------------------------------------------
!
        REAL_I_PARENT=REAL_I_START+(I_CHILD-1)*RATIO_C_P                   !<-- Parent I index coinciding with child domain point
!
!       i_block: IF(REAL_I_PARENT>=R_ITS.AND.REAL_I_PARENT<=R_IEND)THEN    !<-- Column (I) of child's S/N bndry point lies on parent task?
        i_block: IF(REAL_I_PARENT>=R_ITS.AND.REAL_I_PARENT< R_IEND)THEN    !<-- Column (I) of child's S/N bndry point lies on parent task?
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
          IF(FLAG_H_OR_V=='H_POINTS')THEN
            R_JTS=REAL(JTS)-EPS                                            !<-- REAL Jstart of parent task's subdomain for H on B grid
            R_JEND=REAL(MIN(JTE+1,JDE))-EPS                                !<-- Allow search for child H boundary points to go into
                                                                           !    the parent's halo.
          ELSEIF(FLAG_H_OR_V=='V_POINTS')THEN
            R_JTS =REAL(JTS-0.5)-EPS                                       !<-- REAL Jstart of parent task's subdomain for V on B grid
                                                                           !    (-0.5 yields same location on grid as R_JTS for H).
            R_JEND=REAL(MIN(REAL(JTE+0.5),REAL(JDE)))-EPS                  !<-- Use JTE+0.5 to stop V search at the row of the
                                                                           !    northernmost H that is searched; this ensures that
                                                                           !    a parent will send both H and V boundary points.
          ENDIF
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
!***  Also remember that the NMM-B boundary update routines go
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
                  I_SAVE_LO_SOUTH(NC_HOLD_S(NC))=I_CHILD                   !<-- Save westernmost Sbndry I of child task NC
                                                                           !    that is on this parent task.
                ENDIF
                I_SAVE_HI_SOUTH(NC_HOLD_S(NC))=I_CHILD                     !<-- Save easternmost Sbndry I of child task NC
                                                                           !    that is on this parent task.
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
              REAL_J_PARENT=REAL_J_START+(KOUNT_J-1)*RATIO_C_P             !<-- REAL parent J for this child's J
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
          REAL_J_START=PARENT_J_CHILD_NBND
          KOUNT_J=0
!
          IF(FLAG_H_OR_V=='H_POINTS')THEN
            R_JTS=REAL(JTS)+EPS                                            !<-- REAL Jstart of parent task's subdomain for H on B grid
            R_JEND=REAL(MIN(JTE+1,JDE))+EPS                                !<-- Allow search for child H boundary points to go into
                                                                           !    the parent's halo.
          ELSEIF(FLAG_H_OR_V=='V_POINTS')THEN
            R_JTS =REAL(JTS-0.5)+EPS                                       !<-- REAL Jstart of parent task's subdomain for V on B grid
                                                                           !    (-0.5 yields same location on grid as R_JTS for H).
            R_JEND=REAL(MIN(REAL(JTE+0.5),REAL(JDE)))+EPS                  !<-- Use JTE+0.5 to stop V search at the row of the
                                                                           !    northernmost H that is searched; this ensures that
                                                                           !    a parent will send both H and V boundary points.
          ENDIF
!
!-----------------------------------------------------------------------
!
          J_CHILD=JM_END                                                   
!
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
              IF(I_CHILD>=I_LIMIT_LO(NC).AND.                           &  !<-- Does current child boundary point on this
                 I_CHILD<=I_LIMIT_HI(NC)                                &  !    parent task lie on child task "NC"?
                 .AND.                                                  &  !
                 J_CHILD>=JTS_CHILD(NC).AND.                            &  !
                 J_CHILD<=JTE_CHILD(NC))THEN
!
                IF(NC>NC_LAST_N)THEN                                       !<-- Have we encountered a new child task holding this N bndry?
                  NUM_TASKS_SEND_N=NUM_TASKS_SEND_N+1                      !<-- Then increment the N bndry counter of the child tasks
                  LOCAL_TASK_RANK_N(NUM_TASKS_SEND_N)=NC-1                 !<-- Save this child task's local rank
                  NC_LAST_N=NC
                  NC_HOLD_N(NC)=NUM_TASKS_SEND_N
                ENDIF
!
                IF(I_SAVE_LO_NORTH(NC_HOLD_N(NC))<0)THEN
                  I_SAVE_LO_NORTH(NC_HOLD_N(NC))=I_CHILD                   !<-- Save westernmost Nbndry I of child task NC
                                                                           !    that is on this parent task.
                ENDIF
                I_SAVE_HI_NORTH(NC_HOLD_N(NC))=I_CHILD                     !<-- Save easternmost Nbndry I of child task NC
                                                                           !    that is on this parent task.
! 
              ENDIF
! 
            ENDDO
! 
!-----------------------------------------------------------------------
!
            j_north: DO J_CHILD=JM_END,JM_END-N_BLEND+1-N_ADD,-1           !<-- Blending region of child's northern boundary
!
!-----------------------------------------------------------------------
!
              KOUNT_J=KOUNT_J+1
              REAL_J_PARENT=REAL_J_START-(KOUNT_J-1)*RATIO_C_P             !<-- REAL parent J for this child's J
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
      NC_LAST_W=-1
      NC_LAST_E=-1
!
      NUM_TASKS_SEND_W=0                                                   !<-- Parent task sends to this many child tasks on W bndry
      NUM_TASKS_SEND_E=0                                                   !<-- Parent task sends to this many child tasks on E bndry
!
      IF(FLAG_H_OR_V=='H_POINTS')THEN
        R_JTS=REAL(JTS)-EPS                                                !<-- REAL Jstart of parent task's subdomain for H on B grid
      ELSEIF(FLAG_H_OR_V=='V_POINTS')THEN
        R_JTS=REAL(JTS-0.5)-EPS                                            !<-- REAL Jstart of parent task's subdomain for V on B grid
      ENDIF
!
      ARG1=REAL(JTE)+ADD_INC
      ARG2=REAL(JDE)
      R_JEND=MIN(ARG1,ARG2)-EPS                                            !<-- REAL Jend of parent task's region for child W/E boundaries
!
      REAL_J_START=PARENT_J_CHILD_SBND
!
!-----------------------------------------------------------------------
      j_loop: DO J_CHILD=1,JM_END                                          !<-- Loop through child J's across its W/E boundaries
!-----------------------------------------------------------------------

        REAL_J_PARENT=REAL_J_START+(J_CHILD-1)*RATIO_C_P                   !<-- Parent J index coinciding with child domain point
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
          IF(FLAG_H_OR_V=='H_POINTS')THEN
            R_ITS=REAL(ITS)-EPS                                            !<-- REAL Istart of parent task's subdomain for H on B grid
            R_IEND=REAL(MIN(ITE+1,IDE))-EPS                                !<-- Allow search for child H boundary points to go into
                                                                           !    the parent's halo.
          ELSEIF(FLAG_H_OR_V=='V_POINTS')THEN
            R_ITS =REAL(ITS-0.5)-EPS                                       !<-- REAL Istart of parent task's subdomain for V on B grid
                                                                           !    (-0.5 yields same location on grid as R_JTS for H).
            R_IEND=REAL(MIN(REAL(ITE+0.5),REAL(IDE)))-EPS                  !<-- Use ITE+0.5 to stop V search at the row of the
                                                                           !    northernmost H that is searched; this ensures that
                                                                           !    a parent will send both H and V boundary points.
          ENDIF
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
                  J_SAVE_LO_WEST(NC_HOLD_W(NC))=J_CHILD                    !<-- Save southernmost Wbndry I of child task NC
                                                                           !    that is on this parent task.
                ENDIF
                J_SAVE_HI_WEST(NC_HOLD_W(NC))=J_CHILD                      !<-- Save northernmost Wbndry I of child task NC
                                                                           !    that is on this parent task.
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
              REAL_I_PARENT=REAL_I_START+(KOUNT_I-1)*RATIO_C_P             !<-- REAL parent I for this child's I
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
          REAL_I_START=PARENT_I_CHILD_EBND
          KOUNT_I=0
!
          IF(FLAG_H_OR_V=='H_POINTS')THEN
            R_ITS=REAL(ITS)+EPS                                            !<-- REAL Istart of parent task's subdomain for H on B grid
            R_IEND=REAL(MIN(ITE+1,IDE))+EPS                                !<-- Allow search for child H boundary points to go into
                                                                           !    the parent's halo.
          ELSEIF(FLAG_H_OR_V=='V_POINTS')THEN
            R_ITS =REAL(ITS-0.5)+EPS                                       !<-- REAL Istart of parent task's subdomain for V on B grid
                                                                           !    (-0.5 yields same location on grid as R_JTS for H).
            R_IEND=REAL(MIN(REAL(ITE+0.5),REAL(IDE)))+EPS                  !<-- Use ITE+0.5 to stop V search at the row of the
                                                                           !    northernmost H that is searched; this ensures that
                                                                           !    a parent will send both H and V boundary points.
          ENDIF
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
          I_CHILD=IM_END
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
                  J_SAVE_LO_EAST(NC_HOLD_E(NC))=J_CHILD                    !<-- Save southernmost Ebndry I of child task NC
                                                                           !    that is on this parent task.
                ENDIF
                J_SAVE_HI_EAST(NC_HOLD_E(NC))=J_CHILD                      !<-- Save northernmost Ebndry I of child task NC
                                                                           !    that is on this parent task.
!
              ENDIF
!
            ENDDO
!
!-----------------------------------------------------------------------
!
            i_east: DO I_CHILD=IM_END,IM_END-N_BLEND+1-N_ADD,-1            !<-- Blending region of child's eastern boundary
!
              KOUNT_I=KOUNT_I+1
              REAL_I_PARENT=REAL_I_START-(KOUNT_I-1)*RATIO_C_P             !<-- REAL parent I for this child's I
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
      SUBROUTINE POINT_INTERP_DATA_TO_MEMORY(N_CHILD,TIME_FLAG)
!
!-----------------------------------------------------------------------
!***  Create unallocated working pointers for nest boundary variables
!***  and point them into the allocated composite pointer that holds
!***  all of a parent task's data it will send to each child boundary
!***  task it covers.  Nest boundary pressure though must be allocated
!***  because it contains more data than is transferred since we need
!***  extra points in order to do the 4-pt averaging to the nest
!***  boundary V points for hydrostatic balancing of the boundary
!***  data.
!
!***  Only parents execute this routine.
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument Variables
!-----------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: N_CHILD                             !<-- This child is being handled.
!
      CHARACTER(*),INTENT(IN) :: TIME_FLAG                                 !<-- Current or future boundary data for the child?
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I_END_TRANSFER,ITE_CHILD_X,INDX2            &
                           ,J_END_TRANSFER,JTE_CHILD_X                  &
                           ,N,N_TASK                                    &
                           ,NBASE,NBASE_3D,NBASE_EXP                    &
                           ,NCHILD_TASKS                                &
                           ,NLOC_1,NLOC_2,NLOC_2_EXP                    &
                           ,NN,NT,NWORDS
!
      INTEGER(kind=KINT) :: ISTAT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      N=N_CHILD
!
!-----------------------------------------------------------------------
!***  Select the appropriate part of the working array depending on
!***  whether we are now concerned with children's boundaries for
!***  their current time or for their future.
!-----------------------------------------------------------------------
!
      IF(TIME_FLAG=='Future')THEN
        INDX2=1
      ELSEIF(TIME_FLAG=='Current')THEN
        INDX2=2
      ENDIF
!
!-----------------------------------------------------------------------
!***  For each child domain on this parent, create the working pointers
!***  for the nest boundary variables for each child boundary task
!***  and point them into the allocated composite data pointer that
!***  holds all the data for transfer.  Nest boundary pressure is 
!***  treated differently by allocating it and eventually copying it
!***  directly into the composite data pointer.
!
!***  Set logical flags so parent tasks know if they must send any
!***  data at all to any nest boundary tasks.
!
!***  Allocate/nullify new MPI handles for the most recent association
!***  between parent tasks and nest boundary tasks for the ISends of
!***  data to the nest boundaries.
!-----------------------------------------------------------------------
!
!-----------
!***  South
!-----------
!
      south_h: IF(NUM_TASKS_SEND_H_S(N)>0)THEN                             !<-- Parent task has child south boundary H points?
!
        NCHILD_TASKS=NUM_TASKS_SEND_H_S(N)                                 !<-- # of Sbndry tasks on child N to recv H point data
        ALLOCATE(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(1:NCHILD_TASKS))       !<-- 1-D bndry data string for child tasks with Sbndry H points
        ALLOCATE(WORDS_BOUND_H_SOUTH(N)%TASKS(1:NCHILD_TASKS))             !<-- # of words in Sbndry H point 1-D data string
!
        ALLOCATE(PD_B_SOUTH(N)%TASKS(1:NCHILD_TASKS))                      !<-- PD_B_SOUTH for each child task
        ALLOCATE(T_B_SOUTH (N)%TASKS(1:NCHILD_TASKS))                      !<-- T_B_SOUTH for each child task
        ALLOCATE(Q_B_SOUTH (N)%TASKS(1:NCHILD_TASKS))                      !<-- Q_B_SOUTH for each child task
        ALLOCATE(CW_B_SOUTH(N)%TASKS(1:NCHILD_TASKS))                      !<-- CW_B_SOUTH for each child task
!
        nt_south_h: DO NT=1,NCHILD_TASKS
!
          N_TASK=CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NT)+1
          ITE_CHILD_X=CTASK_LIMITS(N)%LIMITS(2,N_TASK)
          I_END_TRANSFER=MIN(ITE_CHILD_X+2,IM_CHILD(N))
!
!------------------------------------------------------------------------
!***  The child I extent of words to be transferred from this parent task
!***  to child task NT is one less than the limit used for saving values
!***  of PDB on the child boundary.  We needed to save PDB at one point
!***  further east than the easternmost V in the segment to be able 
!***  to do 4-pt averaging of PDB onto the V points in order to do
!***  hydrostatic updating of V by the parent.  Now indicate that
!***  reduction in the number of points to be transferred.
!------------------------------------------------------------------------
!
          CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NT)=                   &   !<-- Sbndry I limit for transfer to child
            CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NT)-1
!
          IF(CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NT)==ITE_CHILD_X)            &   !<-- We do not reduce the area for H data   
            CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NT)=                 &   !    transfer if the bndry segment reaches
            CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NT)+1                    !    the physical limit of that bndry
!
          NBASE=CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NT)              &
               -CHILDTASK_H_SAVE(N)%I_LO_SOUTH(NT)+1          
!
          NBASE_3D=LM*NBASE*N_BLEND_H_CHILD(N)
          NWORDS  =(3*LM+1)*NBASE*N_BLEND_H_CHILD(N)                         !<-- # of Sbndry words to transfer from parent to child
          WORDS_BOUND_H_SOUTH(N)%TASKS(NT)=NWORDS                            !<-- Save total number of words
!
          CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA=>NULL()
!
          CALL CHECK_REAL(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA &
                        ,'CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA')
          ALLOCATE(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA(1:NWORDS))    !<-- 1-D bndry data string for child tasks
!                                                                                 with South boundary H points
          NLOC_1=1
          NLOC_2=NLOC_1+NBASE*N_BLEND_H_CHILD(N)-1
!
          NBASE_EXP=CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NT)                   &
                   -CHILDTASK_H_SAVE(N)%I_LO_SOUTH(NT)+1
!
          NLOC_2_EXP=NLOC_1+NBASE_EXP*(N_BLEND_H_CHILD(N)+1)-1               !<-- Extend PD_B_* to allow 4-pt averaging to V pts
!
          PD_B_SOUTH(N)%TASKS(NT)%DATA=>NULL()
!
          CALL CHECK_REAL(PD_B_SOUTH(N)%TASKS(NT)%DATA                   &
                        ,'PD_B_SOUTH(N)%TASKS(NT)%DATA')
          ALLOCATE(PD_B_SOUTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2_EXP),stat=ISTAT)
!
          NLOC_1=NLOC_2+1                                                    !<-- Start at NLOC_2, NOT NLOC_2_EXPAND
          NLOC_2=NLOC_1+NBASE_3D-1
          T_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- T_B_SOUTH storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          Q_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- Q_B_SOUTH storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          CW_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- CW_B_SOUTH storage location
!
        ENDDO  nt_south_h
!
      ELSE  south_h                                                          !<-- Dummy nonzero length 
!
        ALLOCATE(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(1:1))       
        CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(1)%DATA        &
                      ,'CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(1)%DATA')
        ALLOCATE(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(1)%DATA(1:1))
        ALLOCATE(PD_B_SOUTH(N)%TASKS(1:1))              
        PD_B_SOUTH(N)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(PD_B_SOUTH(N)%TASKS(1)%DATA                       &
                      ,'PD_B_SOUTH(N)%TASKS(1)%DATA')
        ALLOCATE(PD_B_SOUTH(N)%TASKS(1)%DATA(1:1),stat=ISTAT)              
        ALLOCATE(T_B_SOUTH (N)%TASKS(1:1))             
        ALLOCATE(Q_B_SOUTH (N)%TASKS(1:1))            
        ALLOCATE(CW_B_SOUTH(N)%TASKS(1:1))           
        ALLOCATE(WORDS_BOUND_H_SOUTH(N)%TASKS(1:1))           
!
      ENDIF south_h
!
!------------------------------------------------------------------------
!
      south_v: IF(NUM_TASKS_SEND_V_S(N)>0)THEN                              !<-- Parent task has child south boundary V points?
!
        NCHILD_TASKS=NUM_TASKS_SEND_V_S(N)
        ALLOCATE(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(1:NCHILD_TASKS))        !<-- 1-D bndry data string for child tasks with Sbndry V points
        ALLOCATE(WORDS_BOUND_V_SOUTH(N)%TASKS(1:NCHILD_TASKS))              !<-- # of words in Sbndry V point 1-D data string
!
        ALLOCATE(PD_B_SOUTH_V(N)%TASKS(1:NCHILD_TASKS))                     !<-- PD_B_SOUTH_V for each child task
        ALLOCATE(U_B_SOUTH(N)%TASKS(1:NCHILD_TASKS))                        !<-- U_B_SOUTH for each child task
        ALLOCATE(V_B_SOUTH(N)%TASKS(1:NCHILD_TASKS))                        !<-- V_B_SOUTH for each child task
!
        DO NT=1,NCHILD_TASKS
          NBASE   =CHILDTASK_V_SAVE(N)%I_HI_SOUTH(NT)                    &
                  -CHILDTASK_V_SAVE(N)%I_LO_SOUTH(NT)+1
          NBASE_3D=LM*NBASE*N_BLEND_V_CHILD(N)
          NWORDS  =2*NBASE_3D                                               !<-- Total number of V boundary words for child task's segment
          WORDS_BOUND_V_SOUTH(N)%TASKS(NT)=NWORDS                           !<-- Save total number of words
!
          CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(NT)%DATA    &
                        ,'CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(NT)%DATA')
          ALLOCATE(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(NT)%DATA(1:NWORDS))   !<-- 1-D bndry data string for child tasks with Sbndry V points
!
          PD_B_SOUTH_V(N)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(PD_B_SOUTH_V(N)%TASKS(NT)%DATA                 &
                        ,'PD_B_SOUTH_V(N)%TASKS(NT)%DATA')
          ALLOCATE(PD_B_SOUTH_V(N)%TASKS(NT)%DATA(1:NBASE*N_BLEND_V_CHILD(N)),stat=ISTAT)
!
          NLOC_1=1
          NLOC_2=NLOC_1+NBASE_3D-1
          U_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- U_B_SOUTH storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          V_B_SOUTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- V_B_SOUTH storage location
!
        ENDDO
!
      ELSE  south_v                                                         !<-- Dummy nonzero length
!
        ALLOCATE(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(1:1))       
        CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(1)%DATA       &
                      ,'CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(1)%DATA')
        ALLOCATE(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(1)%DATA(1:1))
        ALLOCATE(PD_B_SOUTH_V(N)%TASKS(1:1))
        PD_B_SOUTH_V(N)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(PD_B_SOUTH_V(N)%TASKS(1)%DATA                    &
                      ,'PD_B_SOUTH_V(N)%TASKS(1)%DATA')
        ALLOCATE(PD_B_SOUTH_V(N)%TASKS(1)%DATA(1:1),stat=ISTAT)
        ALLOCATE(U_B_SOUTH(N)%TASKS(1:1))
        ALLOCATE(V_B_SOUTH(N)%TASKS(1:1))
        ALLOCATE(WORDS_BOUND_V_SOUTH(N)%TASKS(1:1))           
!
      ENDIF south_v
!
!------------------------------------------------------------------------
!
!-----------
!***  North
!-----------
!
      north_h: IF(NUM_TASKS_SEND_H_N(N)>0)THEN                              !<-- Parent task has child north boundary H points?
!
        NCHILD_TASKS=NUM_TASKS_SEND_H_N(N)
        ALLOCATE(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(1:NCHILD_TASKS)      &  !<-- 1-D bndry data string for child tasks with Nbndry H points
                                                            ,stat=ISTAT)
        ALLOCATE(WORDS_BOUND_H_NORTH(N)%TASKS(1:NCHILD_TASKS))              !<-- # of words in Nbndry H point 1-D data string
!
        ALLOCATE(PD_B_NORTH(N)%TASKS(1:NCHILD_TASKS))                       !<-- PD_B_NORTH for each child task
        ALLOCATE(T_B_NORTH (N)%TASKS(1:NCHILD_TASKS))                       !<-- T_B_NORTH for each child task
        ALLOCATE(Q_B_NORTH (N)%TASKS(1:NCHILD_TASKS))                       !<-- Q_B_NORTH for each child task
        ALLOCATE(CW_B_NORTH(N)%TASKS(1:NCHILD_TASKS))                       !<-- CW_B_NORTH for each child task
!
        nt_north_h: DO NT=1,NCHILD_TASKS
!
          N_TASK=CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NT)+1
          ITE_CHILD_X=CTASK_LIMITS(N)%LIMITS(2,N_TASK)
          I_END_TRANSFER=MIN(ITE_CHILD_X+2,IM_CHILD(N))
!
!------------------------------------------------------------------------
!***  The child I extent of words to be transferred from this parent task
!***  to child task NT is one less than the limit used for saving values
!***  of PDB on the child boundary.  We needed to save PDB at one point
!***  further east than the easternmost V in the segment to be able 
!***  to do 4-pt averaging of PDB onto the V points in order to do
!***  hydrostatic updating of V by the parent.  Now indicate that
!***  reduction in the number of points to be transferred.
!------------------------------------------------------------------------
!
          CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NT)=                   &  !<-- Nbndry I limit for transfer to child
            CHILDTASK_H_SAVE(N)%I_HI_NORTH(NT)-1
!
          IF(CHILDTASK_H_SAVE(N)%I_HI_NORTH(NT)==ITE_CHILD_X)            &  !<-- We do not reduce the area for H data
            CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NT)=                 &  !    transfer if the bndry segment reaches
            CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NT)+1                   !    the physical limit of that bndry
!
          NBASE=CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NT)              &
               -CHILDTASK_H_SAVE(N)%I_LO_NORTH(NT)+1
!
          NBASE_3D=LM*NBASE*N_BLEND_H_CHILD(N)
          NWORDS  =(3*LM+1)*NBASE*N_BLEND_H_CHILD(N)                        !<-- # of Nbndry words to transfer from parent to child
          WORDS_BOUND_H_NORTH(N)%TASKS(NT)=NWORDS                           !<-- Save total number of words
!
          CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA    &
                        ,'CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA')
!
          ALLOCATE(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA(1:NWORDS) &  !<-- 1-D bndry data string for child tasks
                                                             ,stat=ISTAT)   !     with north boundary H points.
          NLOC_1=1
          NLOC_2=NLOC_1+NBASE*N_BLEND_H_CHILD(N)-1
!
          NBASE_EXP=CHILDTASK_H_SAVE(N)%I_HI_NORTH(NT)                   &
                   -CHILDTASK_H_SAVE(N)%I_LO_NORTH(NT)+1
!
          NLOC_2_EXP=NLOC_1+NBASE_EXP*(N_BLEND_H_CHILD(N)+1)-1              !<-- Extend PD_B_* by one row to allow 4-pt averaging to V pts
!
          PD_B_NORTH(N)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(PD_B_NORTH(N)%TASKS(NT)%DATA                   &
                        ,'PD_B_NORTH(N)%TASKS(NT)%DATA')
!
          ALLOCATE(PD_B_NORTH(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2_EXP)       &
                                                           ,stat=ISTAT)
!
          NLOC_1=NLOC_2+1                                                   !<-- Start at NLOC_2, NOT NLOC_2_EXPAND
          NLOC_2=NLOC_1+NBASE_3D-1
          T_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- T_B_NORTH storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          Q_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- Q_B_NORTH storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          CW_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- CW_B_NORTH storage location
!
        ENDDO  nt_north_h
!
      ELSE  north_h                                                         !<-- Dummy nonzero length
!
        ALLOCATE(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(1:1))       
        CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(1)%DATA       &
                      ,'CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(1)%DATA')
        ALLOCATE(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(1)%DATA(1:1)         &
                                                           ,stat=ISTAT)
        ALLOCATE(PD_B_NORTH(N)%TASKS(1:1))
        PD_B_NORTH(N)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(PD_B_NORTH(N)%TASKS(1)%DATA                      &
                      ,'PD_B_NORTH(N)%TASKS(1)%DATA')
        ALLOCATE(PD_B_NORTH(N)%TASKS(1)%DATA(1:1),stat=ISTAT)
        ALLOCATE(T_B_NORTH (N)%TASKS(1:1))
        ALLOCATE(Q_B_NORTH (N)%TASKS(1:1))
        ALLOCATE(CW_B_NORTH(N)%TASKS(1:1))
        ALLOCATE(WORDS_BOUND_H_NORTH(N)%TASKS(1:1))           
!
      ENDIF north_h
!
!------------------------------------------------------------------------
!
      north_v: IF(NUM_TASKS_SEND_V_N(N)>0)THEN                              !<-- Parent task has child north boundary V points?
!
        NCHILD_TASKS=NUM_TASKS_SEND_V_N(N)
        ALLOCATE(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(1:NCHILD_TASKS))        !<-- 1-D bndry data string for child tasks with Nbndry V points
        ALLOCATE(WORDS_BOUND_V_NORTH(N)%TASKS(1:NCHILD_TASKS))              !<-- # of words in Nbndry V point 1-D data string
!
        ALLOCATE(PD_B_NORTH_V(N)%TASKS(1:NCHILD_TASKS))                     !<-- PD_B_NORTH_V for each child task
        ALLOCATE(U_B_NORTH(N)%TASKS(1:NCHILD_TASKS))                        !<-- U_B_NORTH for each child task
        ALLOCATE(V_B_NORTH(N)%TASKS(1:NCHILD_TASKS))                        !<-- V_B_NORTH for each child task
!
        DO NT=1,NCHILD_TASKS
          NBASE=CHILDTASK_V_SAVE(N)%I_HI_NORTH(NT)                       &
               -CHILDTASK_V_SAVE(N)%I_LO_NORTH(NT)+1
          NBASE_3D=LM*NBASE*N_BLEND_V_CHILD(N)
          NWORDS=2*NBASE_3D                                                 !<-- Total number of V boundary words for child task's segment
          WORDS_BOUND_V_NORTH(N)%TASKS(NT)=NWORDS                           !<-- Save total number of words
!
          CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA    &
                        ,'CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA')
          ALLOCATE(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA(1:NWORDS) &  !<-- 1-D bndry data string for child tasks with Nbndry V points
                                                             ,stat=ISTAT)   
!
          PD_B_NORTH_V(N)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(PD_B_NORTH_V(N)%TASKS(NT)%DATA                 &
                        ,'PD_B_NORTH_V(N)%TASKS(NT)%DATA')
!
          ALLOCATE(PD_B_NORTH_V(N)%TASKS(NT)%DATA(1:NBASE*N_BLEND_V_CHILD(N)) &
                                                                 ,stat=ISTAT)
!
          NLOC_1=1
          NLOC_2=NLOC_1+NBASE_3D-1
          U_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- U_B_NORTH storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          V_B_NORTH(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- V_B_NORTH storage location
!
        ENDDO
!
      ELSE  north_v                                                         !<-- Dummy nonzero length
!
        ALLOCATE(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(1:1))       
        CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(1)%DATA       &
                      ,'CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(1)%DATA')
        ALLOCATE(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(1)%DATA(1:1)         &
                                                           ,stat=ISTAT)
        ALLOCATE(PD_B_NORTH_V(N)%TASKS(1:1))
        PD_B_NORTH_V(N)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(PD_B_NORTH_V(N)%TASKS(1)%DATA                    &
                      ,'PD_B_NORTH_V(N)%TASKS(1)%DATA')
        ALLOCATE(PD_B_NORTH_V(N)%TASKS(1)%DATA(1:1),stat=ISTAT)
        ALLOCATE(U_B_NORTH(N)%TASKS(1:1))
        ALLOCATE(V_B_NORTH(N)%TASKS(1:1))
        ALLOCATE(WORDS_BOUND_V_NORTH(N)%TASKS(1:1))           
!
      ENDIF north_v
!
!------------------------------------------------------------------------
!
!----------
!***  West
!----------
!
      west_h: IF(NUM_TASKS_SEND_H_W(N)>0)THEN                               !<-- Parent task has child west boundary H points?
!
        NCHILD_TASKS=NUM_TASKS_SEND_H_W(N)
        ALLOCATE(CHILD_BOUND_H_WEST(N,INDX2)%TASKS(1:NCHILD_TASKS))         !<-- 1-D bndry data string for child tasks with Wbndry H points
        ALLOCATE(WORDS_BOUND_H_WEST(N)%TASKS(1:NCHILD_TASKS))               !<-- # of words in Wbndry H point 1-D data string
!
        ALLOCATE(PD_B_WEST(N)%TASKS(1:NCHILD_TASKS))                        !<-- PD_B_WEST for each child task
        ALLOCATE(T_B_WEST (N)%TASKS(1:NCHILD_TASKS))                        !<-- T_B_WEST for each child task
        ALLOCATE(Q_B_WEST (N)%TASKS(1:NCHILD_TASKS))                        !<-- Q_B_WEST for each child task
        ALLOCATE(CW_B_WEST(N)%TASKS(1:NCHILD_TASKS))                        !<-- CW_B_WEST for each child task
!
        nt_west_h: DO NT=1,NCHILD_TASKS
!
          N_TASK=CHILDTASK_BNDRY_H_RANKS(N)%WEST(NT)+1
          JTE_CHILD_X=CTASK_LIMITS(N)%LIMITS(4,N_TASK)
          J_END_TRANSFER=MIN(JTE_CHILD_X+2,JM_CHILD(N))
!
!------------------------------------------------------------------------
!***  The child J extent of words to be transferred from this parent task
!***  to child task NT is one less than the limit used for saving values
!***  of PDB on the child boundary.  We needed to save PDB at one point
!***  further north than the northernmost V in the segment to be able 
!***  to do 4-pt averaging of PDB onto the V points in order to do
!***  hydrostatic updating of V by the parent.  Now indicate that
!***  reduction in the number of points to be transferred.
!------------------------------------------------------------------------
!
          CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NT)=                    &  !<-- Wbndry J limit for transfer to child
            CHILDTASK_H_SAVE(N)%J_HI_WEST(NT)-1
!
          IF(CHILDTASK_H_SAVE(N)%J_HI_WEST(NT)==JTE_CHILD_X)             &  !<-- We do not reduce the area for H data
            CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NT)=                  &  !    transfer if the bndry segment reaches
            CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NT)+1                    !    the physical limit of that bndry
!
          NBASE=CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NT)               &
               -CHILDTASK_H_SAVE(N)%J_LO_WEST(NT)+1
!
          NBASE_3D=LM*NBASE*N_BLEND_H_CHILD(N)
          NWORDS  =(3*LM+1)*NBASE*N_BLEND_H_CHILD(N)                        !<-- # of Wbndry words to transfer from parent to child
          WORDS_BOUND_H_WEST(N)%TASKS(NT)=NWORDS                            !<-- Save total number of words
!
          CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA     &
                        ,'CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA')
          ALLOCATE(CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA(1:NWORDS))    !<-- 1-D bndry data string for child tasks
!                                                                                with west boundary H points
          NLOC_1=1
          NLOC_2=NLOC_1+NBASE*N_BLEND_H_CHILD(N)-1
!
          NBASE_EXP=CHILDTASK_H_SAVE(N)%J_HI_WEST(NT)                    &
                   -CHILDTASK_H_SAVE(N)%J_LO_WEST(NT)+1
!
          NLOC_2_EXP=NLOC_1+NBASE_EXP*(N_BLEND_H_CHILD(N)+1)-1              !<-- Extend PD_B_* by one row to allow 4-pt averaging to V pts
!
          PD_B_WEST(N)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(PD_B_WEST(N)%TASKS(NT)%DATA                    &
                        ,'PD_B_WEST(N)%TASKS(NT)%DATA')
          ALLOCATE(PD_B_WEST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2_EXP),stat=ISTAT)
!
          NLOC_1=NLOC_2+1                                                   !<-- Start at NLOC_2, NOT NLOC_2_EXPAND
          NLOC_2=NLOC_1+NBASE_3D-1
          T_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- T_B_WEST storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          Q_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- Q_B_WEST storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          CW_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- CW_B_WEST storage location
!
        ENDDO  nt_west_h
!
      ELSE  west_h                                                          !<-- Dummy nonzero length
!
        ALLOCATE(CHILD_BOUND_H_WEST(N,INDX2)%TASKS(1:1))       
        CHILD_BOUND_H_WEST(N,INDX2)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(CHILD_BOUND_H_WEST(N,INDX2)%TASKS(1)%DATA        &
                      ,'CHILD_BOUND_H_WEST(N,INDX2)%TASKS(1)%DATA')
        ALLOCATE(CHILD_BOUND_H_WEST(N,INDX2)%TASKS(1)%DATA(1:1))
        ALLOCATE(PD_B_WEST(N)%TASKS(1:1))
        PD_B_WEST(N)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(PD_B_WEST(N)%TASKS(1)%DATA                       &
                      ,'PD_B_WEST(N)%TASKS(1)%DATA')
        ALLOCATE(PD_B_WEST(N)%TASKS(1)%DATA(1:1),stat=ISTAT)
        ALLOCATE(T_B_WEST (N)%TASKS(1:1))
        ALLOCATE(Q_B_WEST (N)%TASKS(1:1))
        ALLOCATE(CW_B_WEST(N)%TASKS(1:1))
        ALLOCATE(WORDS_BOUND_H_WEST(N)%TASKS(1:1))           
!
      ENDIF west_h
!
!------------------------------------------------------------------------
!
      west_v: IF(NUM_TASKS_SEND_V_W(N)>0)THEN                               !<-- Parent task has child west boundary V points?
!
        NCHILD_TASKS=NUM_TASKS_SEND_V_W(N)
        ALLOCATE(CHILD_BOUND_V_WEST(N,INDX2)%TASKS(1:NCHILD_TASKS))         !<-- 1-D bndry data string for child tasks with Wbndry V points
        ALLOCATE(WORDS_BOUND_V_WEST(N)%TASKS(1:NCHILD_TASKS))               !<-- # of words in Wbndry V point 1-D data string
!
        ALLOCATE(PD_B_WEST_V(N)%TASKS(1:NCHILD_TASKS))                      !<-- PD_B_WEST_V for each child task
        ALLOCATE(U_B_WEST(N)%TASKS(1:NCHILD_TASKS))                         !<-- U_B_WEST for each child task
        ALLOCATE(V_B_WEST(N)%TASKS(1:NCHILD_TASKS))                         !<-- V_B_WEST for each child task
!
        DO NT=1,NCHILD_TASKS
          NBASE=CHILDTASK_V_SAVE(N)%J_HI_WEST(NT)                        &  
               -CHILDTASK_V_SAVE(N)%J_LO_WEST(NT)+1
          NBASE_3D=LM*NBASE*N_BLEND_V_CHILD(N)
          NWORDS=2*NBASE_3D                                                 !<-- Total number of V boundary words for child task's segment
          WORDS_BOUND_V_WEST(N)%TASKS(NT)=NWORDS                            !<-- Save total number of words
!
          CHILD_BOUND_V_WEST(N,INDX2)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(CHILD_BOUND_V_WEST(N,INDX2)%TASKS(NT)%DATA     &
                        ,'CHILD_BOUND_V_WEST(N,INDX2)%TASKS(NT)%DATA')
          ALLOCATE(CHILD_BOUND_V_WEST(N,INDX2)%TASKS(NT)%DATA(1:NWORDS))    !<-- 1-D bndry data string for child tasks with Wbndry V points
!
          PD_B_WEST_V(N)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(PD_B_WEST_V(N)%TASKS(NT)%DATA                  &
                        ,'PD_B_WEST_V(N)%TASKS(NT)%DATA')
          ALLOCATE(PD_B_WEST_V(N)%TASKS(NT)%DATA(1:NBASE*N_BLEND_V_CHILD(N)),stat=ISTAT)
!
          NLOC_1=1
          NLOC_2=NLOC_1+NBASE_3D-1
          U_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_WEST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- U_B_WEST storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          V_B_WEST(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_WEST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- V_B_WEST storage location
!
        ENDDO
!
      ELSE  west_v                                                          !<-- Dummy nonzero length
!
        ALLOCATE(CHILD_BOUND_V_WEST(N,INDX2)%TASKS(1:1))       
        CHILD_BOUND_V_WEST(N,INDX2)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(CHILD_BOUND_V_WEST(N,INDX2)%TASKS(1)%DATA        &
                      ,'CHILD_BOUND_V_WEST(N,INDX2)%TASKS(1)%DATA')
        ALLOCATE(CHILD_BOUND_V_WEST(N,INDX2)%TASKS(1)%DATA(1:1))
        ALLOCATE(PD_B_WEST_V(N)%TASKS(1:1))
        PD_B_WEST_V(N)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(PD_B_WEST_V(N)%TASKS(1)%DATA                     &
                      ,'PD_B_WEST_V(N)%TASKS(1)%DATA')
        ALLOCATE(PD_B_WEST_V(N)%TASKS(1)%DATA(1:1),stat=ISTAT)
        ALLOCATE(U_B_WEST(N)%TASKS(1:1))
        ALLOCATE(V_B_WEST(N)%TASKS(1:1))
        ALLOCATE(WORDS_BOUND_V_WEST(N)%TASKS(1:1))           
!
      ENDIF west_v
!
!------------------------------------------------------------------------
!
!----------
!***  East
!----------
!
      east_h: IF(NUM_TASKS_SEND_H_E(N)>0)THEN                               !<-- Parent task has child east boundary H points?
!
        NCHILD_TASKS=NUM_TASKS_SEND_H_E(N)
        ALLOCATE(CHILD_BOUND_H_EAST(N,INDX2)%TASKS(1:NCHILD_TASKS))         !<-- 1-D bndry data string for child tasks with Ebndry H points
        ALLOCATE(WORDS_BOUND_H_EAST(N)%TASKS(1:NCHILD_TASKS))               !<-- # of words in Ebndry H point 1-D data string
!
        ALLOCATE(PD_B_EAST(N)%TASKS(1:NCHILD_TASKS))                        !<-- PD_B_EAST for each child task
        ALLOCATE(T_B_EAST (N)%TASKS(1:NCHILD_TASKS))                        !<-- T_B_EAST for each child task
        ALLOCATE(Q_B_EAST (N)%TASKS(1:NCHILD_TASKS))                        !<-- Q_B_EAST for each child task
        ALLOCATE(CW_B_EAST(N)%TASKS(1:NCHILD_TASKS))                        !<-- CW_B_EAST for each child task
!
        nt_east_h: DO NT=1,NCHILD_TASKS
!
          N_TASK=CHILDTASK_BNDRY_H_RANKS(N)%EAST(NT)+1
          JTE_CHILD_X=CTASK_LIMITS(N)%LIMITS(4,N_TASK)
          J_END_TRANSFER=MIN(JTE_CHILD_X+2,JM_CHILD(N))
!
!------------------------------------------------------------------------
!***  The child J extent of words to be transferred from this parent task
!***  to child task NT is one less than the limit used for saving values
!***  of PDB on the child boundary.  We needed to save PDB at one point
!***  further north than the northernmost V in the segment to be able 
!***  to do 4-pt averaging of PDB onto the V points in order to do
!***  hydrostatic updating of V by the parent.  Now indicate that
!***  reduction in the number of points to be transferred.
!------------------------------------------------------------------------
!
          CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NT)=                    &  !<-- Ebndry J limit for transfer to child
            CHILDTASK_H_SAVE(N)%J_HI_EAST(NT)-1
!
          IF(CHILDTASK_H_SAVE(N)%J_HI_EAST(NT)==JTE_CHILD_X)             &  !<-- We do not reduce the area for H data
            CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NT)=                  &  !    transfer if the bndry segment reaches
            CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NT)+1                    !    the physical limit of that bndry
!
          NBASE=CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NT)               &
               -CHILDTASK_H_SAVE(N)%J_LO_EAST(NT)+1
!
          NBASE_3D=LM*NBASE*N_BLEND_H_CHILD(N)
          NWORDS  =(3*LM+1)*NBASE*N_BLEND_H_CHILD(N)                        !<-- # of Ebndry words to transfer from parent to child
          WORDS_BOUND_H_EAST(N)%TASKS(NT)=NWORDS                            !<-- Save total number of words
!
          CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA     &
                        ,'CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA')
          ALLOCATE(CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA(1:NWORDS))    !<-- 1-D bndry data string for child tasks
!                                                                                with east boundary H points
          NLOC_1=1
          NLOC_2=NLOC_1+NBASE*N_BLEND_H_CHILD(N)-1
!
          NBASE_EXP=CHILDTASK_H_SAVE(N)%J_HI_EAST(NT)                    &
                   -CHILDTASK_H_SAVE(N)%J_LO_EAST(NT)+1
!
          NLOC_2_EXP=NLOC_1+NBASE_EXP*(N_BLEND_H_CHILD(N)+1)-1              !<-- Extend PD_B_* by one row to allow 4-pt averaging to V pts
!
          PD_B_EAST(N)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(PD_B_EAST(N)%TASKS(NT)%DATA                    &
                        ,'PD_B_EAST(N)%TASKS(NT)%DATA')
          ALLOCATE(PD_B_EAST(N)%TASKS(NT)%DATA(NLOC_1:NLOC_2_EXP),stat=ISTAT)
!
          NLOC_1=NLOC_2+1                                                   !<-- Start at NLOC_2, NOT NLOC_2_EXPAND
          NLOC_2=NLOC_1+NBASE_3D-1
          T_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- T_B_EAST storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          Q_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- Q_B_EAST storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          CW_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)  !<-- CW_B_EAST storage location
!
        ENDDO  nt_east_h
!
      ELSE east_h                                                           !<-- Dummy nonzero length
!
        ALLOCATE(CHILD_BOUND_H_EAST(N,INDX2)%TASKS(1:1))       
        CHILD_BOUND_H_EAST(N,INDX2)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(CHILD_BOUND_H_EAST(N,INDX2)%TASKS(1)%DATA        &
                      ,'CHILD_BOUND_H_EAST(N,INDX2)%TASKS(1)%DATA')
        ALLOCATE(CHILD_BOUND_H_EAST(N,INDX2)%TASKS(1)%DATA(1:1))
        ALLOCATE(PD_B_EAST(N)%TASKS(1:1))
        PD_B_EAST(N)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(PD_B_EAST(N)%TASKS(1)%DATA                       &
                      ,'PD_B_EAST(N)%TASKS(1)%DATA')
        ALLOCATE(PD_B_EAST(N)%TASKS(1)%DATA(1:1))
        ALLOCATE(T_B_EAST (N)%TASKS(1:1))
        ALLOCATE(Q_B_EAST (N)%TASKS(1:1))
        ALLOCATE(CW_B_EAST(N)%TASKS(1:1))
        ALLOCATE(WORDS_BOUND_H_EAST(N)%TASKS(1:1))           
!
      ENDIF east_h
!
!-----------------------------------------------------------------------
!
       east_v: IF(NUM_TASKS_SEND_V_E(N)>0)THEN                             !<-- Parent task has child east boundary V points?
!
        NCHILD_TASKS=NUM_TASKS_SEND_V_E(N)
        ALLOCATE(CHILD_BOUND_V_EAST(N,INDX2)%TASKS(1:NCHILD_TASKS))        !<-- 1-D bndry data string for child tasks with Ebndry V points
        ALLOCATE(WORDS_BOUND_V_EAST(N)%TASKS(1:NCHILD_TASKS))              !<-- # of words in Ebndry V point 1-D data string
!
        ALLOCATE(PD_B_EAST_V(N)%TASKS(1:NCHILD_TASKS))                     !<-- PD_B_EAST_V for each child task
        ALLOCATE(U_B_EAST(N)%TASKS(1:NCHILD_TASKS))                        !<-- U_B_EAST for each child task
        ALLOCATE(V_B_EAST(N)%TASKS(1:NCHILD_TASKS))                        !<-- V_B_EAST for each child task
!
        DO NT=1,NCHILD_TASKS
          NBASE=CHILDTASK_V_SAVE(N)%J_HI_EAST(NT)                       &
               -CHILDTASK_V_SAVE(N)%J_LO_EAST(NT)+1
          NBASE_3D=LM*NBASE*N_BLEND_V_CHILD(N)
          NWORDS=2*NBASE_3D                                                !<-- Total number of V boundary words for child task's segment
          WORDS_BOUND_V_EAST(N)%TASKS(NT)=NWORDS                           !<-- Save total number of words
!
          CHILD_BOUND_V_EAST(N,INDX2)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(CHILD_BOUND_V_EAST(N,INDX2)%TASKS(NT)%DATA    &
                        ,'CHILD_BOUND_V_EAST(N,INDX2)%TASKS(NT)%DATA')
          ALLOCATE(CHILD_BOUND_V_EAST(N,INDX2)%TASKS(NT)%DATA(1:NWORDS))   !<-- 1-D bndry data string for child tasks with Ebndry V points
!
          PD_B_EAST_V(N)%TASKS(NT)%DATA=>NULL()
          CALL CHECK_REAL(PD_B_EAST_V(N)%TASKS(NT)%DATA                 &
                        ,'PD_B_EAST_V(N)%TASKS(NT)%DATA')
          ALLOCATE(PD_B_EAST_V(N)%TASKS(NT)%DATA(1:NBASE*N_BLEND_V_CHILD(N)),stat=ISTAT)
!
          NLOC_1=1
          NLOC_2=NLOC_1+NBASE_3D-1
          U_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_EAST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- U_B_EAST storage location
!
          NLOC_1=NLOC_2+1
          NLOC_2=NLOC_1+NBASE_3D-1
          V_B_EAST(N)%TASKS(NT)%DATA=>CHILD_BOUND_V_EAST(N,INDX2)%TASKS(NT)%DATA(NLOC_1:NLOC_2)   !<-- V_B_EAST storage location
!
        ENDDO
!
      ELSE  east_v                                                         !<-- Dummy nonzero length
!
        ALLOCATE(CHILD_BOUND_V_EAST(N,INDX2)%TASKS(1:1))       
        CHILD_BOUND_V_EAST(N,INDX2)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(CHILD_BOUND_V_EAST(N,INDX2)%TASKS(1)%DATA       &
                      ,'CHILD_BOUND_V_EAST(N,INDX2)%TASKS(1)%DATA')
        ALLOCATE(CHILD_BOUND_V_EAST(N,INDX2)%TASKS(1)%DATA(1:1))
        ALLOCATE(PD_B_EAST_V(N)%TASKS(1:1))
        PD_B_EAST_V(N)%TASKS(1)%DATA=>NULL()
        CALL CHECK_REAL(PD_B_EAST_V(N)%TASKS(1)%DATA                    &
                      ,'PD_B_EAST_V(N)%TASKS(1)%DATA')
        ALLOCATE(PD_B_EAST_V(N)%TASKS(1)%DATA(1:1),stat=ISTAT)
        ALLOCATE(U_B_EAST(N)%TASKS(1:1))
        ALLOCATE(V_B_EAST(N)%TASKS(1:1))
        ALLOCATE(WORDS_BOUND_V_EAST(N)%TASKS(1:1))           
!
      ENDIF east_v
!
!-----------------------------------------------------------------------
!***  Here we set logical flags so each parent tasks knows whether or
!***  not it must send data to any side of child N's boundary.
!-----------------------------------------------------------------------
!
      IF(NUM_TASKS_SEND_H_S(N)>0.OR.                                    &
         NUM_TASKS_SEND_H_N(N)>0.OR.                                    &
         NUM_TASKS_SEND_H_W(N)>0.OR.                                    &
         NUM_TASKS_SEND_H_E(N)>0.OR.                                    &
         NUM_TASKS_SEND_V_S(N)>0.OR.                                    &
         NUM_TASKS_SEND_V_N(N)>0.OR.                                    &
         NUM_TASKS_SEND_V_W(N)>0.OR.                                    &
         NUM_TASKS_SEND_V_E(N)>0)THEN
!
        SEND_CHILD_DATA(N)=.TRUE.
!
      ELSE
        SEND_CHILD_DATA(N)=.FALSE.
      ENDIF
!
!-----------------------------------------------------------------------
!***  Allocate and nullify the new handles created for ISends between
!***  parent tasks and nest boundary tasks.  That association of tasks
!***  obviously changes each time the nests move.
!-----------------------------------------------------------------------
!
!-------------------------------
!***  For child boundary, south
!-------------------------------
!
      IF(NUM_TASKS_SEND_H_S(N)>0)THEN
        ALLOCATE(HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_H_S(N)))
        DO NN=1,NUM_TASKS_SEND_H_S(N)
          HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
        ENDDO
      ELSE
        ALLOCATE(HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV(1:1))
        HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
      ENDIF
!
      IF(NUM_TASKS_SEND_V_S(N)>0)THEN
        ALLOCATE(HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_V_S(N)))
        DO NN=1,NUM_TASKS_SEND_V_S(N)
          HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
        ENDDO
      ELSE
        ALLOCATE(HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV(1:1))
        HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
      ENDIF
!
!-------------------------------
!***  For child boundary, north
!-------------------------------
!
      IF(NUM_TASKS_SEND_H_N(N)>0)THEN
        ALLOCATE(HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_H_N(N)))
        DO NN=1,NUM_TASKS_SEND_H_N(N)
          HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
        ENDDO
      ELSE
        ALLOCATE(HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV(1:1))
        HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
      ENDIF
!
      IF(NUM_TASKS_SEND_V_N(N)>0)THEN
        ALLOCATE(HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_V_N(N)))
        DO NN=1,NUM_TASKS_SEND_V_N(N)
          HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
        ENDDO
      ELSE
        ALLOCATE(HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV(1:1))
        HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
      ENDIF
!
!------------------------------
!***  For child boundary, west
!------------------------------
!
      IF(NUM_TASKS_SEND_H_W(N)>0)THEN
        ALLOCATE(HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_H_W(N)))
        DO NN=1,NUM_TASKS_SEND_H_W(N)
          HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
        ENDDO
      ELSE
        ALLOCATE(HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV(1:1))
        HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
      ENDIF
!
      IF(NUM_TASKS_SEND_V_W(N)>0)THEN
        ALLOCATE(HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_V_W(N)))
        DO NN=1,NUM_TASKS_SEND_V_W(N)
          HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
        ENDDO
      ELSE
        ALLOCATE(HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV(1:1))
        HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
      ENDIF
!
!------------------------------
!***  For child boundary, east
!------------------------------
!
      IF(NUM_TASKS_SEND_H_E(N)>0)THEN
        ALLOCATE(HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_H_E(N)))
        DO NN=1,NUM_TASKS_SEND_H_E(N)
          HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
        ENDDO
      ELSE
        ALLOCATE(HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV(1:1))
        HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
      ENDIF
!
      IF(NUM_TASKS_SEND_V_E(N)>0)THEN
        ALLOCATE(HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV(1:NUM_TASKS_SEND_V_E(N)))
        DO NN=1,NUM_TASKS_SEND_V_E(N)
          HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV(NN)=MPI_REQUEST_NULL
        ENDDO
      ELSE
        ALLOCATE(HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV(1:1))
        HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV(1)=MPI_REQUEST_NULL
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE POINT_INTERP_DATA_TO_MEMORY
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_SENDS_CHILD_DATA_LIMITS(N_CHILD)
!
!-----------------------------------------------------------------------
!***  Parents send children basic bookkeeping information needed
!***  for the exchange of boundary data during the integration.
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument Variables
!-----------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: N_CHILD                              !<-- The child being considered by this parent
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: ID_CHILDTASK,ID_DOM,IERR                    &
                           ,MYPE,N,NT,NTX
!
      INTEGER(kind=KINT),DIMENSION(4) :: INFO=(/-9999,0,0,0/)
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: STATUS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Parent tasks send three pieces of information to child tasks 
!***  on child domain boundaries so those child tasks will be
!***  able to receive boundary data and use it properly:
!
!     (1) The local rank of the sending parent.
!     (2) The child boundary tasks' starting (I,J) on the parent task.
!     (3) The child boundary tasks' ending (I,J) on the parent task.
!
!***  The child task must be able to know if the data it receives
!***  pertains to south boundary H or V points, north boundary
!***  H or V, points, etc.  Thus the MPI tag will indicate
!***  the boundary's side and variable type.
!
!      11111 --> South H
!      22222 --> South V
!      33333 --> North H
!      44444 --> North V
!      55555 --> West H  
!      66666 --> West V
!      77777 --> East H  
!      88888 --> East V
!
!*** (The child tasks know which side of their domain's boundary they
!***  are on of course but since a corner task is on more than one side,
!***  the above tags indicating the side are used for all child tasks.)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Reinitialize the first word of the data packet.  It will be 
!***  changed to the valid task ID (i.e., a non-negative integer)
!***  of a parent task that has boundary data to send to a child. 
!-----------------------------------------------------------------------
!
      INFO(1)=-9999
!
!-----------------------------------------------------------------------
!
      N=N_CHILD
!
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK(COMM_TO_MY_CHILDREN(N),MYPE,IERR)                 !<-- Obtain rank of parent task
      ID_DOM=MY_CHILDREN_ID(N)
!
!-------------
!***  South H
!-------------
!
      sh_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                             !<-- Each parent task loops through all tasks on each child
!
        INFO(1)=-1
!
        IF(NUM_TASKS_SEND_H_S(N)>0)THEN                                
!
          DO NTX=1,NUM_TASKS_SEND_H_S(N)                                   !<-- Look for a child task with south boundary H points
            ID_CHILDTASK=CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NTX)
!              
            IF(NT==ID_CHILDTASK)THEN                                       !<-- If yes, we found a child task w/ south boundary H points
              INFO(1)=MYPE                                                 !<-- Save the parent task rank
              INFO(2)=CHILDTASK_H_SAVE(N)%I_LO_SOUTH(NTX)                  !<-- Save the starting index of boundary segment on child
              INFO(3)=CHILDTASK_H_SAVE(N)%I_HI_SOUTH_TRANSFER(NTX)         !<-- Save the ending index of boundary segment on child
              INFO(4)=CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NTX)                  !<-- Save the ending index of expanded boundary segment on child
!
              CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK             &  !<-- Parent task sends the key data to the child Sbndry task
                           ,11111                                       &  !<-- Tag for south boundary H points
                           ,COMM_TO_MY_CHILDREN(N),IERR)
!
              CYCLE sh_loop                                                !<-- Move on to next child task
            ENDIF
          ENDDO
!
        ENDIF
!
        ID_CHILDTASK=NT                                                    !<-- This child task has no south boundary H points
        CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK                   &  !<-- So send this child task some dummy information
                     ,11111                                             &  !<-- Tag for south boundary H points
                     ,COMM_TO_MY_CHILDREN(N),IERR)
!
      ENDDO sh_loop
!
!-------------
!***  South V
!-------------
!
      sv_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                             !<-- Each parent task loops through all tasks on each child
!
        INFO(1)=-1
!
        IF(NUM_TASKS_SEND_V_S(N)>0)THEN                                
!
          DO NTX=1,NUM_TASKS_SEND_V_S(N)                                   !<-- Look for a child task with south boundary V points
            ID_CHILDTASK=CHILDTASK_BNDRY_V_RANKS(N)%SOUTH(NTX)
!              
            IF(NT==ID_CHILDTASK)THEN                                       !<-- If yes, we found a child task w/ south boundary V points
              INFO(1)=MYPE                                                 !<-- Save the parent task rank
              INFO(2)=CHILDTASK_V_SAVE(N)%I_LO_SOUTH(NTX)                  !<-- Save the starting index of boundary segment on child
              INFO(3)=CHILDTASK_V_SAVE(N)%I_HI_SOUTH(NTX)                  !<-- Save the ending index of boundary segment on child
!
              CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK             &  !<-- Parent task sends the key data to the child Sbndry task
                           ,22222                                       &  !<-- Tag for south boundary V points
                           ,COMM_TO_MY_CHILDREN(N),IERR)
!
              CYCLE sv_loop                                                !<-- Move on to next child task
            ENDIF
          ENDDO
!
        ENDIF
!
        ID_CHILDTASK=NT                                                    !<-- This child task has no south boundary H points
        CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK                   &  !<-- So send this child task some dummy information
                     ,22222                                             &  !<-- Tag for south boundary V points
                     ,COMM_TO_MY_CHILDREN(N),IERR)
!
      ENDDO sv_loop
!
!-------------
!***  North H
!-------------
!
      nh_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                             !<-- Each parent task loops through all tasks on each child
!
        INFO(1)=-1
!
        IF(NUM_TASKS_SEND_H_N(N)>0)THEN                                
!
          DO NTX=1,NUM_TASKS_SEND_H_N(N)                                   !<-- Look for a child task with north boundary H points
            ID_CHILDTASK=CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NTX)
!              
            IF(NT==ID_CHILDTASK)THEN                                       !<-- If yes, we found a child task w/ north boundary H points
              INFO(1)=MYPE                                                 !<-- Save the parent task rank
              INFO(2)=CHILDTASK_H_SAVE(N)%I_LO_NORTH(NTX)                  !<-- Save the starting index of boundary segment on child
              INFO(3)=CHILDTASK_H_SAVE(N)%I_HI_NORTH_TRANSFER(NTX)         !<-- Save the ending index of boundary segment on child
              INFO(4)=CHILDTASK_H_SAVE(N)%I_HI_NORTH(NTX)                  !<-- Save the ending index of expanded boundary segment on child
!
              CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK             &  !<-- Parent task sends the key data to the child Nbndry task
                           ,33333                                       &  !<-- Tag for north boundary H points
                           ,COMM_TO_MY_CHILDREN(N),IERR)
!
              CYCLE nh_loop                                                !<-- Move on to next child task
            ENDIF
          ENDDO
!
        ENDIF
!
        ID_CHILDTASK=NT                                                    !<-- This child task has no north boundary H points
        CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK                   &  !<-- So send this child task some dummy information
                     ,33333                                             &  !<-- Tag for north boundary H points
                     ,COMM_TO_MY_CHILDREN(N),IERR)
!
      ENDDO nh_loop
!
!-------------
!***  North V
!-------------
!
      nv_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                             !<-- Each parent task loops through all tasks on each child
!
        INFO(1)=-1
!
        IF(NUM_TASKS_SEND_V_N(N)>0)THEN                                
!
          DO NTX=1,NUM_TASKS_SEND_V_N(N)                                   !<-- Look for a child task with north boundary V points
            ID_CHILDTASK=CHILDTASK_BNDRY_V_RANKS(N)%NORTH(NTX)
!              
            IF(NT==ID_CHILDTASK)THEN                                       !<-- If yes, we found a child task w/ north boundary V points
              INFO(1)=MYPE                                                 !<-- Save the parent task rank
              INFO(2)=CHILDTASK_V_SAVE(N)%I_LO_NORTH(NTX)                  !<-- Save the starting index of boundary segment on child
              INFO(3)=CHILDTASK_V_SAVE(N)%I_HI_NORTH(NTX)                  !<-- Save the ending index of boundary segment on child
!
              CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK             &  !<-- Parent task sends the key data to the child Nbndry task
                           ,44444                                       &  !<-- Tag for north boundary V points
                           ,COMM_TO_MY_CHILDREN(N),IERR)
!
              CYCLE nv_loop                                                !<-- Move on to next child task
            ENDIF
          ENDDO
!
        ENDIF
!
        ID_CHILDTASK=NT                                                    !<-- This child task has no north boundary V points
        CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK                   &  !<-- So send this child task some dummy information
                     ,44444                                             &  !<-- Tag for north boundary V points
                     ,COMM_TO_MY_CHILDREN(N),IERR)
!
      ENDDO nv_loop
!
!------------
!***  West H
!------------
!
      wh_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                             !<-- Each parent task loops through all tasks on each child
!
        INFO(1)=-1
!
        IF(NUM_TASKS_SEND_H_W(N)>0)THEN                                
!
          DO NTX=1,NUM_TASKS_SEND_H_W(N)                                   !<-- Look for a child task with west boundary H points
            ID_CHILDTASK=CHILDTASK_BNDRY_H_RANKS(N)%WEST(NTX)
!              
            IF(NT==ID_CHILDTASK)THEN                                       !<-- If yes, we found a child task w/ west boundary H points
              INFO(1)=MYPE                                                 !<-- Save the parent task rank
              INFO(2)=CHILDTASK_H_SAVE(N)%J_LO_WEST(NTX)                   !<-- Save the starting index of boundary segment on child
              INFO(3)=CHILDTASK_H_SAVE(N)%J_HI_WEST_TRANSFER(NTX)          !<-- Save the ending index of boundary segment on child
              INFO(4)=CHILDTASK_H_SAVE(N)%J_HI_WEST(NTX)                   !<-- Save the ending index of expanded boundary segment on child
!
              CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK             &  !<-- Parent task sends the key data to the child Wbndry task
                           ,55555                                       &  !<-- Tag for west boundary H points
                           ,COMM_TO_MY_CHILDREN(N),IERR)
!
              CYCLE wh_loop                                                !<-- Move on to next child task
            ENDIF
          ENDDO
!
        ENDIF
!
        ID_CHILDTASK=NT                                                    !<-- This child task has no west boundary H points
        CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK                   &  !<-- So send this child task some dummy information
                     ,55555                                             &  !<-- Tag for west boundary H points
                     ,COMM_TO_MY_CHILDREN(N),IERR)
!
      ENDDO wh_loop
!
!------------
!***  West V
!------------
!
      wv_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                             !<-- Each parent task loops through all tasks on each child
!
        INFO(1)=-1
!
        IF(NUM_TASKS_SEND_V_W(N)>0)THEN                                
!
          DO NTX=1,NUM_TASKS_SEND_V_W(N)                                   !<-- Look for a child task with west boundary V points
            ID_CHILDTASK=CHILDTASK_BNDRY_V_RANKS(N)%WEST(NTX)
!              
            IF(NT==ID_CHILDTASK)THEN                                       !<-- If yes, we found a child task w/ west boundary V points
              INFO(1)=MYPE                                                 !<-- Save the parent task rank
              INFO(2)=CHILDTASK_V_SAVE(N)%J_LO_WEST(NTX)                   !<-- Save the starting index of boundary segment on child
              INFO(3)=CHILDTASK_V_SAVE(N)%J_HI_WEST(NTX)                   !<-- Save the ending index of boundary segment on child
!
              CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK             &  !<-- Parent task sends the key data to the child Wbndry task
                           ,66666                                       &  !<-- Tag for west boundary V points
                           ,COMM_TO_MY_CHILDREN(N),IERR)
!
              CYCLE wv_loop                                                !<-- Move on to next child task
            ENDIF
          ENDDO
!
        ENDIF
!
        ID_CHILDTASK=NT                                                    !<-- This child task has no west boundary V points
        CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK                   &  !<-- So send this child task some dummy information
                     ,66666                                             &  !<-- Tag for west boundary V points
                     ,COMM_TO_MY_CHILDREN(N),IERR)
!
      ENDDO wv_loop
!
!------------
!***  East H
!------------
!
      eh_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                             !<-- Each parent task loops through all tasks on each child
!
        INFO(1)=-1
!
        IF(NUM_TASKS_SEND_H_E(N)>0)THEN
!
          DO NTX=1,NUM_TASKS_SEND_H_E(N)                                   !<-- Look for a child task with east boundary H points
            ID_CHILDTASK=CHILDTASK_BNDRY_H_RANKS(N)%EAST(NTX)
!              
            IF(NT==ID_CHILDTASK)THEN                                       !<-- If yes, we found a child task w/ east boundary H points
              INFO(1)=MYPE                                                 !<-- Save the parent task rank
              INFO(2)=CHILDTASK_H_SAVE(N)%J_LO_EAST(NTX)                   !<-- Save the starting index of boundary segment on child
              INFO(3)=CHILDTASK_H_SAVE(N)%J_HI_EAST_TRANSFER(NTX)          !<-- Save the ending index of boundary segment on child
              INFO(4)=CHILDTASK_H_SAVE(N)%J_HI_EAST(NTX)                   !<-- Save the ending index of expanded boundary segment on child
!
              CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK             &  !<-- Parent task sends the key data to the child Ebndry task
                           ,77777                                       &  !<-- Tag for east boundary H points
                           ,COMM_TO_MY_CHILDREN(N),IERR)
!
              CYCLE eh_loop                                                !<-- Move on to next child task
            ENDIF
          ENDDO
!
        ENDIF
!
        ID_CHILDTASK=NT                                                    !<-- This child task has no east boundary H points
        CALL MPI_SEND(INFO,4,MPI_INTEGER,ID_CHILDTASK                   &  !<-- So send this child task some dummy information
                     ,77777                                             &  !<-- Tag for east boundary H points
                     ,COMM_TO_MY_CHILDREN(N),IERR)
!
      ENDDO eh_loop
!
!------------
!***  East V
!------------
!
      ev_loop: DO NT=0,FTASKS_DOMAIN(ID_DOM)-1                             !<-- Each parent task loops through all tasks on each child
!
        INFO(1)=-1
!
        IF(NUM_TASKS_SEND_V_E(N)>0)THEN                                
!
          DO NTX=1,NUM_TASKS_SEND_V_E(N)                                   !<-- Look for a child task with east boundary V points
            ID_CHILDTASK=CHILDTASK_BNDRY_V_RANKS(N)%EAST(NTX)
!              
            IF(NT==ID_CHILDTASK)THEN                                       !<-- If yes, we found a child task w/ east boundary V points
              INFO(1)=MYPE                                                 !<-- Save the parent task rank
              INFO(2)=CHILDTASK_V_SAVE(N)%J_LO_EAST(NTX)                   !<-- Save the starting index of boundary segment on child
              INFO(3)=CHILDTASK_V_SAVE(N)%J_HI_EAST(NTX)                   !<-- Save the ending index of boundary segment on child
!
              CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK             &  !<-- Parent task sends the key data to the child Ebndry task
                           ,88888                                       &  !<-- Tag for east boundary V points
                           ,COMM_TO_MY_CHILDREN(N),IERR)
!
              CYCLE ev_loop                                                !<-- Move on to next child task
            ENDIF
          ENDDO
!
        ENDIF
!
        ID_CHILDTASK=NT                                                    !<-- This child task has no east boundary V points
        CALL MPI_SEND(INFO,3,MPI_INTEGER,ID_CHILDTASK                   &  !<-- So send this child task some dummy information
                     ,88888                                             &  !<-- Tag for east boundary V points
                     ,COMM_TO_MY_CHILDREN(N),IERR)
!
      ENDDO ev_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_SENDS_CHILD_DATA_LIMITS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CHILD_RECVS_CHILD_DATA_LIMITS(EXP_STATE)
!
!-----------------------------------------------------------------------
!***  Children receive from their parents basic bookkeeping information
!***  needed for the exchange of boundary data during the integration.
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument Variables
!-----------------------
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE                           !<-- Parent-Child Coupler export state
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: ID_DOM,KOUNT,LENGTH,N,N1,NBASE
!
      INTEGER(kind=KINT) :: ILIM_HI,ILIM_LO,JLIM_HI,JLIM_LO
!
      INTEGER(kind=KINT) :: IERR,ISTAT,RC,RC_LIMITS 
!
      INTEGER(kind=KINT),DIMENSION(4) :: INFO=(/-9999,0,0,0/)
      INTEGER(kind=KINT),DIMENSION(4) :: TEMP
      INTEGER(kind=KINT),DIMENSION(4,2) :: PARENT_INFO
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: STATUS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  The child tasks receive the key information from their
!***  parent tasks.  At this point the child tasks do not know
!***  the local ranks of the parent tasks that will be sending
!***  information to them thus MPI_ANY_SOURCE is used in the
!***  receive.  However this means that when there are two
!***  parent tasks sending to a nest task (rather than only one)
!***  then the two overlap points on the nest boundary segment
!***  computed by the parent tasks will ultimately have values
!***  depending on which parent task's preliminary information
!***  is received last in the MPI_ANY_SOURCE Recv below.  Since
!***  the values in those overlap points are not bit identical
!***  then any successive runs can have slightly different answers.
!***  To avoid that happening when two parent tasks are sending,
!***  the child task will receive their ranks and then put them
!***  in ascending order so that all subsequent updates of the
!***  nest boundary overlap points are always done in the same
!***  way regardless of the order the preliminary information
!***  is received with MPI_ANY_SOURCE.
!
!***  All nests execute this routine once during the Init step and
!***  then again on those parent timesteps during the Run step when
!***  the child has moved. 
!-----------------------------------------------------------------------
!
      call mpi_comm_rank(comm_to_my_parent,mype,ierr)      !<-- Obtain local rank of child task
!     write(0,*)' PARENT_TO_CHILD_DATA_LIMITS child ready to recv from parent mype=',mype
      ID_DOM=ID_PARENTS(MY_DOMAIN_ID)                                      !<-- Domain ID of this child's parent
!
!-------------
!***  South H
!-------------
!
      KOUNT=0
      INDX_MIN_H%SOUTH= 1000000
      INDX_MAX_H%SOUTH=-1000000
!
      DO N=1,FTASKS_DOMAIN(ID_DOM)                                         !<-- Child task loops through its parent's tasks
!
        CALL MPI_RECV(INFO                                              &  !<-- Receive data packet from each parent task
                     ,4                                                 &  !<-- # of words in data packet
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,MPI_ANY_SOURCE                                    &  !<-- Accept data from any parent task that is sending
                     ,11111                                             &  !<-- Tag used for south boundary H points
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator between child and parent
                     ,STATUS                                            &  !<-- Status of Recv
                     ,IERR)
!
        IF(INFO(1)>=0)THEN                                                 !<-- If yes, this parent task sent key preliminary bndry info
          KOUNT=KOUNT+1
          DO N1=1,4
            PARENT_INFO(N1,KOUNT)=INFO(N1)                                 !<-- Save the data from that parent task
          ENDDO
        ENDIF
!
      ENDDO
!
      IF(KOUNT==2)THEN                                                     !<-- Nest task recvs data from two parent tasks
        IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                          !<-- Data recvd from 'out of order' parent tasks
!
          DO N1=1,4                                                        !<-- Save parent data in order of ascending task IDs
            TEMP(N1)         =PARENT_INFO(N1,1)                            !
            PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                            !
            PARENT_INFO(N1,2)=TEMP(N1)                                     !<--
          ENDDO
!
        ENDIF
      ENDIF
!
      IF(KOUNT>0)THEN
        DO N=1,KOUNT
          PARENT_TASK(N)%SOUTH_H%ID_SOURCE   =PARENT_INFO(1,N)             !<-- Rank of parent task that will send Sboundary H data segment
          PARENT_TASK(N)%SOUTH_H%INDX_START  =PARENT_INFO(2,N)             !<-- Istart on child grid of the boundary data segment
          PARENT_TASK(N)%SOUTH_H%INDX_END    =PARENT_INFO(3,N)             !<-- Iend on child grid of the boundary data segment
          PARENT_TASK(N)%SOUTH_H%INDX_END_EXP=PARENT_INFO(4,N)             !<-- Iend on child grid of the expanded boundary data segment
!
          NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_H
          LENGTH=(3*LM+1)*NBASE
!
          PARENT_TASK(N)%SOUTH_H%LENGTH=LENGTH                             !<-- # of words in this parent task's datastring of child bndry
!
          ALLOCATE(PARENT_TASK(N)%SOUTH_H%STRING(1:LENGTH))                !<-- Sboundary H datastring to be received from parent task
!
          INDX_MIN_H%SOUTH=MIN(INDX_MIN_H%SOUTH,PARENT_INFO(2,N))          !<-- Starting child I for union of parent task segments sent  
          INDX_MAX_H%SOUTH=MAX(INDX_MAX_H%SOUTH,PARENT_INFO(3,N))          !<-- Ending child I for union of parent task segments sent
        ENDDO
      ENDIF
!
      NUM_PARENT_TASKS_SENDING_H%SOUTH=KOUNT
!!!   LENGTH_BND_SEG_H%SOUTH=INDX_MAX_H%SOUTH-INDX_MIN_H%SOUTH+1
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
        ALLOCATE(BOUND_1D_SOUTH_H(1:LENGTH),stat=ISTAT)                    !<-- 1-D combined H-point data on child task's Sbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Lower I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_LO_SOUTH_H'                  &  !<-- Name of the boundary array's lower I limit
                              ,value=ILIM_LO                            &  !<-- The boundary array's lower I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Upper I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_HI_SOUTH_H'                  &  !<-- Name of the boundary array's upper I limit
                              ,value=ILIM_HI                            &  !<-- The boundary array's upper I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Lower J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_LO_SOUTH_H'                  &  !<-- Name of the boundary array's lower J limit
                              ,value=JLIM_LO                            &  !<-- The boundary array's lower J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Upper J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_HI_SOUTH_H'                  &  !<-- Name of the boundary array's upper J limit
                              ,value=JLIM_HI                            &  !<-- The boundary array's upper J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF south_h
!
!-----------------------------------------------------------------------
!
!-------------
!***  South V
!-------------
!
      KOUNT=0
      INDX_MIN_V%SOUTH= 1000000
      INDX_MAX_V%SOUTH=-1000000
!
      DO N=1,FTASKS_DOMAIN(ID_DOM)                                         !<-- Child task loops through its parent's tasks
!
        CALL MPI_RECV(INFO                                              &  !<-- Receive data packet from each parent task
                     ,3                                                 &  !<-- # of words in data packet
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,MPI_ANY_SOURCE                                    &  !<-- Accept data from any parent task that is sending
                     ,22222                                             &  !<-- Tag used for south boundary V points
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator between child and parent
                     ,STATUS                                            &  !<-- Status of Recv
                     ,IERR)
!
        IF(INFO(1)>=0)THEN                                                 !<-- If yes, this parent task has key preliminary bndry info
          KOUNT=KOUNT+1
          DO N1=1,3
            PARENT_INFO(N1,KOUNT)=INFO(N1)                                 !<-- Save the data from that parent task
          ENDDO
        ENDIF
!
      ENDDO
!
      IF(KOUNT==2)THEN                                                     !<-- Nest task recvs data from two parent tasks
        IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                          !<-- Data recvd from 'out of order' parent tasks
!
          DO N1=1,3                                                        !<-- Save parent data in order of ascending task IDs
            TEMP(N1)         =PARENT_INFO(N1,1)                            !
            PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                            !
            PARENT_INFO(N1,2)=TEMP(N1)                                     !<--
          ENDDO
!
        ENDIF
      ENDIF
!
      IF(KOUNT>0)THEN
        DO N=1,KOUNT
          PARENT_TASK(N)%SOUTH_V%ID_SOURCE =PARENT_INFO(1,N)               !<-- Rank of parent task that will send Sboundary V data segment
          PARENT_TASK(N)%SOUTH_V%INDX_START=PARENT_INFO(2,N)               !<-- Istart on child grid of the boundary data segment
          PARENT_TASK(N)%SOUTH_V%INDX_END  =PARENT_INFO(3,N)               !<-- Iend on child grid of the boundary data segment
!
          NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_V
          LENGTH=2*LM*NBASE
!
          PARENT_TASK(N)%SOUTH_V%LENGTH=LENGTH                             !<-- # of words in this parent task's datastring of child bndry
!
          ALLOCATE(PARENT_TASK(N)%SOUTH_V%STRING(1:LENGTH))                !<-- Sboundary V datastring to be received from parent task
!
          INDX_MIN_V%SOUTH=MIN(INDX_MIN_V%SOUTH,PARENT_INFO(2,N))          !<-- Starting child I for union of parent task segments sent
          INDX_MAX_V%SOUTH=MAX(INDX_MAX_V%SOUTH,PARENT_INFO(3,N))          !<-- Ending child I for union of parent task segments sent
        ENDDO
      ENDIF
!
      NUM_PARENT_TASKS_SENDING_V%SOUTH=KOUNT
!!!   LENGTH_BND_SEG_V%SOUTH=INDX_MAX_V%SOUTH-INDX_MIN_V%SOUTH+1
!
      south_v: IF(NUM_PARENT_TASKS_SENDING_V%SOUTH>0)THEN                  !<-- Does this child task recv any Sboundary V data from parent?
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
        ALLOCATE(BOUND_1D_SOUTH_V(1:LENGTH),stat=ISTAT)                    !<-- 1-D combined V-point data on child task's Sbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Lower I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_LO_SOUTH_V'                  &  !<-- Name of the boundary array's lower I limit
                              ,value=ILIM_LO                            &  !<-- The boundary array's lower I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Upper I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_HI_SOUTH_V'                  &  !<-- Name of the boundary array's upper I limit
                              ,value=ILIM_HI                            &  !<-- The boundary array's upper I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Lower J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_LO_SOUTH_V'                  &  !<-- Name of the boundary array's lower J limit
                              ,value=JLIM_LO                            &  !<-- The boundary array's lower J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Upper J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_HI_SOUTH_V'                  &  !<-- Name of the boundary array's upper J limit
                              ,value=JLIM_HI                            &  !<-- The boundary array's upper J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF south_v
!
!-----------------------------------------------------------------------
!
!-------------
!***  North H
!-------------
!
      KOUNT=0
      INDX_MIN_H%NORTH= 1000000
      INDX_MAX_H%NORTH=-1000000
!
      DO N=1,FTASKS_DOMAIN(ID_DOM)                                         !<-- Child task loops through its parent's tasks
!
        CALL MPI_RECV(INFO                                              &  !<-- Receive data packet from each parent task
                     ,4                                                 &  !<-- # of words in data packet
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,MPI_ANY_SOURCE                                    &  !<-- Accept data from any parent task that is sending
                     ,33333                                             &  !<-- Tag used for north boundary H points
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator between child and parent
                     ,STATUS                                            &  !<-- Status of Recv
                     ,IERR)
!
        IF(INFO(1)>=0)THEN                                                 !<-- If yes, this parent task has key preliminary bndry info
          KOUNT=KOUNT+1
          DO N1=1,4
            PARENT_INFO(N1,KOUNT)=INFO(N1)                                 !<-- Save the data from that parent task
          ENDDO
        ENDIF
!
      ENDDO
!
      IF(KOUNT==2)THEN                                                     !<-- Nest task recvs data from two parent tasks
        IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                          !<-- Data recvd from 'out of order' parent tasks
!
          DO N1=1,4                                                        !<-- Save parent data in order of ascending task IDs
            TEMP(N1)         =PARENT_INFO(N1,1)                            !
            PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                            !
            PARENT_INFO(N1,2)=TEMP(N1)                                     !<--
          ENDDO
!
        ENDIF
      ENDIF
!
      IF(KOUNT>0)THEN
        DO N=1,KOUNT
          PARENT_TASK(N)%NORTH_H%ID_SOURCE   =PARENT_INFO(1,N)             !<-- Rank of parent task that will send Nboundary H data segment
          PARENT_TASK(N)%NORTH_H%INDX_START  =PARENT_INFO(2,N)             !<-- Istart on child grid of the boundary data segment
          PARENT_TASK(N)%NORTH_H%INDX_END    =PARENT_INFO(3,N)             !<-- Iend on child grid of the boundary data segment
          PARENT_TASK(N)%NORTH_H%INDX_END_EXP=PARENT_INFO(4,N)             !<-- Iend on child grid of the expanded boundary data segment
!
          NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_H
          LENGTH=(3*LM+1)*NBASE
!
          PARENT_TASK(N)%NORTH_H%LENGTH=LENGTH                             !<-- # of words in this parent task's datastring of child bndry
!
          ALLOCATE(PARENT_TASK(N)%NORTH_H%STRING(1:LENGTH))                !<-- Nboundary H datastring to be received from parent task
!
          INDX_MIN_H%NORTH=MIN(INDX_MIN_H%NORTH,PARENT_INFO(2,N))          !<-- Starting child I for union of parent task segments sent
          INDX_MAX_H%NORTH=MAX(INDX_MAX_H%NORTH,PARENT_INFO(3,N))          !<-- Ending child I for union of parent task segments sent
        ENDDO
      ENDIF
!
      NUM_PARENT_TASKS_SENDING_H%NORTH=KOUNT
!!!   LENGTH_BND_SEG_H%NORTH=INDX_MAX_H%NORTH-INDX_MIN_H%NORTH+1
!
      north_h: IF(NUM_PARENT_TASKS_SENDING_H%NORTH>0)THEN                  !<-- Does this child task recv Nboundary H data from parent?
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
        ALLOCATE(BOUND_1D_NORTH_H(1:LENGTH),stat=ISTAT)                    !<-- 1-D combined H-point data on child task's Nbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Lower I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_LO_NORTH_H'                  &  !<-- Name of the boundary array's lower I limit
                              ,value=ILIM_LO                            &  !<-- The boundary array's lower I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Upper I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_HI_NORTH_H'                  &  !<-- Name of the boundary array's upper I limit
                              ,value=ILIM_HI                            &  !<-- The boundary array's upper I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Lower J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_LO_NORTH_H'                  &  !<-- Name of the boundary array's lower J limit
                              ,value=JLIM_LO                            &  !<-- The boundary array's lower J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Upper J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_HI_NORTH_H'                  &  !<-- Name of the boundary array's upper J limit
                              ,value=JLIM_HI                            &  !<-- The boundary array's upper J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF north_h
!
!-----------------------------------------------------------------------
!
!-------------
!***  North V
!-------------
!
      KOUNT=0
      INDX_MIN_V%NORTH= 1000000
      INDX_MAX_V%NORTH=-1000000
!
      DO N=1,FTASKS_DOMAIN(ID_DOM)                                         !<-- Child task loops through its parent's tasks
!
        CALL MPI_RECV(INFO                                              &  !<-- Receive data packet from each parent task
                     ,3                                                 &  !<-- # of words in data packet
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,MPI_ANY_SOURCE                                    &  !<-- Accept data from any parent task that is sending
                     ,44444                                             &  !<-- Tag used for north boundary V points
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator between child and parent
                     ,STATUS                                            &  !<-- Status of Recv
                     ,IERR)
!
        IF(INFO(1)>=0)THEN                                                 !<-- If yes, this parent task has key preliminary bndry info
          KOUNT=KOUNT+1
!   if(kount>2)then
!     write(0,*)' BUG: exceeded two parent tasks sending to this child bndry task'
!   endif
          DO N1=1,3
            PARENT_INFO(N1,KOUNT)=INFO(N1)                                 !<-- Save the data from that parent task
          ENDDO
!!!!!!!!!!!!!!!!!!!!!!!debug
!                    else
!   write(0,*)' PARENT_TO_CHILD_DATA_LIMITS child recvd dummy north V from parent task #n=',n,' with id=',-1*info(1)
!!!!!!!!!!!!!!!!!!!!!!!debug
        ENDIF
!
      ENDDO
!
      IF(KOUNT==2)THEN                                                     !<-- Nest task recvs data from two parent tasks
        IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                          !<-- Data recvd from 'out of order' parent tasks
!
          DO N1=1,3                                                        !<-- Save parent data in order of ascending task IDs
            TEMP(N1)         =PARENT_INFO(N1,1)                            !
            PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                            !
            PARENT_INFO(N1,2)=TEMP(N1)                                     !<--
          ENDDO
!
        ENDIF
      ENDIF
!
      IF(KOUNT>0)THEN
        DO N=1,KOUNT
          PARENT_TASK(N)%NORTH_V%ID_SOURCE =PARENT_INFO(1,N)               !<-- Rank of parent task that will send Nboundary V data segment
          PARENT_TASK(N)%NORTH_V%INDX_START=PARENT_INFO(2,N)               !<-- Istart on child grid of the boundary data segment
          PARENT_TASK(N)%NORTH_V%INDX_END  =PARENT_INFO(3,N)               !<-- Iend on child grid of the boundary data segment
!
          NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_V
          LENGTH=2*LM*NBASE
!
          PARENT_TASK(N)%NORTH_V%LENGTH=LENGTH                             !<-- # of words in this parent task's datastring of child bndry
!
          ALLOCATE(PARENT_TASK(N)%NORTH_V%STRING(1:LENGTH))                !<-- Nboundary V datastring to be received from parent task
!
          INDX_MIN_V%NORTH=MIN(INDX_MIN_V%NORTH,PARENT_INFO(2,N))          !<-- Starting child I for union of parent task segments sent
          INDX_MAX_V%NORTH=MAX(INDX_MAX_V%NORTH,PARENT_INFO(3,N))          !<-- Ending child I for union of parent task segments sent
        ENDDO
      ENDIF
!
      NUM_PARENT_TASKS_SENDING_V%NORTH=KOUNT
!!!   LENGTH_BND_SEG_V%NORTH=INDX_MAX_V%NORTH-INDX_MIN_V%NORTH+1
!
      north_v: IF(NUM_PARENT_TASKS_SENDING_V%NORTH>0)THEN                  !<-- Does this child task recv any Nboundary V data from parent?
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
        ALLOCATE(BOUND_1D_NORTH_V(1:LENGTH),stat=ISTAT)                    !<-- 1-D combined V-point data on child task's Nbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Lower I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_LO_NORTH_V'                  &  !<-- Name of the boundary array's lower I limit
                              ,value=ILIM_LO                            &  !<-- The boundary array's lower I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Upper I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_HI_NORTH_V'                  &  !<-- Name of the boundary array's upper I limit
                              ,value=ILIM_HI                            &  !<-- The boundary array's upper I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Lower J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_LO_NORTH_V'                  &  !<-- Name of the boundary array's lower J limit
                              ,value=JLIM_LO                            &  !<-- The boundary array's lower J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Upper J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_HI_NORTH_V'                  &  !<-- Name of the boundary array's upper J limit
                              ,value=JLIM_HI                            &  !<-- The boundary array's upper J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF north_v
!
!-----------------------------------------------------------------------
!
!------------
!***  West H
!------------
!
      KOUNT=0
      INDX_MIN_H%WEST= 1000000
      INDX_MAX_H%WEST=-1000000
!
      DO N=1,FTASKS_DOMAIN(ID_DOM)                                         !<-- Child task loops through its parent's tasks
!
        CALL MPI_RECV(INFO                                              &  !<-- Receive data packet from each parent task
                     ,4                                                 &  !<-- # of words in data packet
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,MPI_ANY_SOURCE                                    &  !<-- Accept data from any parent task that is sending
                     ,55555                                             &  !<-- Tag used for west boundary H points
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator between child and parent
                     ,STATUS                                            &  !<-- Status of Recv
                     ,IERR)
!
        IF(INFO(1)>=0)THEN                                                 !<-- If yes, this parent task has key preliminary bndry info
          KOUNT=KOUNT+1
          DO N1=1,4
            PARENT_INFO(N1,KOUNT)=INFO(N1)                                 !<-- Save the data from that parent task
          ENDDO
        ENDIF
!
      ENDDO
!
      IF(KOUNT==2)THEN                                                     !<-- Nest task recvs data from two parent tasks
        IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                          !<-- Data recvd from 'out of order' parent tasks
!
          DO N1=1,4                                                        !<-- Save parent data in order of ascending task IDs
            TEMP(N1)         =PARENT_INFO(N1,1)                            !
            PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                            !
            PARENT_INFO(N1,2)=TEMP(N1)                                     !<--
          ENDDO
!
        ENDIF
      ENDIF
!
      IF(KOUNT>0)THEN
        DO N=1,KOUNT
          PARENT_TASK(N)%WEST_H%ID_SOURCE   =PARENT_INFO(1,N)              !<-- Rank of parent task that will send Wboundary H data segment
          PARENT_TASK(N)%WEST_H%INDX_START  =PARENT_INFO(2,N)              !<-- Jstart on child grid of the boundary data segment
          PARENT_TASK(N)%WEST_H%INDX_END    =PARENT_INFO(3,N)              !<-- Jend on child grid of the boundary data segment
          PARENT_TASK(N)%WEST_H%INDX_END_EXP=PARENT_INFO(4,N)              !<-- Jend on child grid of the expanded boundary data segment
!
          NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_H
          LENGTH=(3*LM+1)*NBASE
!
          PARENT_TASK(N)%WEST_H%LENGTH=LENGTH                              !<-- # of words in this parent task's datastring of child bndry
!
          ALLOCATE(PARENT_TASK(N)%WEST_H%STRING(1:LENGTH))                 !<-- Wboundary H datastring to be received from parent task
!
          INDX_MIN_H%WEST=MIN(INDX_MIN_H%WEST,PARENT_INFO(2,N))            !<-- Starting child J for union of parent task segments sent
          INDX_MAX_H%WEST=MAX(INDX_MAX_H%WEST,PARENT_INFO(3,N))            !<-- Ending child J for union of parent task segments sent
        ENDDO
      ENDIF
!
      NUM_PARENT_TASKS_SENDING_H%WEST=KOUNT
!!!   LENGTH_BND_SEG_H%WEST=INDX_MAX_H%WEST-INDX_MIN_H%WEST+1
!
      west_h: IF(NUM_PARENT_TASKS_SENDING_H%WEST>0)THEN                    !<-- Does this child task recv Wboundary H data from parent?
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
        ALLOCATE(BOUND_1D_WEST_H(1:LENGTH),stat=ISTAT)                     !<-- 1-D combined H-point data on child task's Wbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Lower I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_LO_WEST_H'                   &  !<-- Name of the boundary array's lower I limit
                              ,value=ILIM_LO                            &  !<-- The boundary array's lower I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Upper I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_HI_WEST_H'                   &  !<-- Name of the boundary array's upper I limit
                              ,value=ILIM_HI                            &  !<-- The boundary array's upper I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Lower J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_LO_WEST_H'                   &  !<-- Name of the boundary array's lower J limit
                              ,value=JLIM_LO                            &  !<-- The boundary array's lower J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Upper J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_HI_WEST_H'                   &  !<-- Name of the boundary array's upper J limit
                              ,value=JLIM_HI                            &  !<-- The boundary array's upper J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF west_h
!
!-----------------------------------------------------------------------
!
!------------
!***  West V
!------------
!
      KOUNT=0
      INDX_MIN_V%WEST= 1000000
      INDX_MAX_V%WEST=-1000000
!
      DO N=1,FTASKS_DOMAIN(ID_DOM)                                         !<-- Child task loops through its parent's tasks
!
        CALL MPI_RECV(INFO                                              &  !<-- Receive data packet from each parent task
                     ,3                                                 &  !<-- # of words in data packet
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,MPI_ANY_SOURCE                                    &  !<-- Accept data from any parent task that is sending
                     ,66666                                             &  !<-- Tag used for west boundary V points
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator between child and parent
                     ,STATUS                                            &  !<-- Status of Recv
                     ,IERR)
!
        IF(INFO(1)>=0)THEN                                                 !<-- If yes, this parent task has key preliminary bndry info
          KOUNT=KOUNT+1
          DO N1=1,3
            PARENT_INFO(N1,KOUNT)=INFO(N1)                                 !<-- Save the data from that parent task
          ENDDO
        ENDIF
!
      ENDDO
!
      IF(KOUNT==2)THEN                                                     !<-- Nest task recvs data from two parent tasks
        IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                          !<-- Data recvd from 'out of order' parent tasks
!
          DO N1=1,3                                                        !<-- Save parent data in order of ascending task IDs
            TEMP(N1)         =PARENT_INFO(N1,1)                            !
            PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                            !
            PARENT_INFO(N1,2)=TEMP(N1)                                     !<--
          ENDDO
!
        ENDIF
      ENDIF
!
      IF(KOUNT>0)THEN
        DO N=1,KOUNT
          PARENT_TASK(N)%WEST_V%ID_SOURCE =PARENT_INFO(1,N)                !<-- Rank of parent task that will send Wboundary V data segment
          PARENT_TASK(N)%WEST_V%INDX_START=PARENT_INFO(2,N)                !<-- Jstart on child grid of the boundary data segment
          PARENT_TASK(N)%WEST_V%INDX_END  =PARENT_INFO(3,N)                !<-- Jend on child grid of the boundary data segment
!
          NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_V
          LENGTH=2*LM*NBASE
!
          PARENT_TASK(N)%WEST_V%LENGTH=LENGTH                              !<-- # of words in this parent task's datastring of child bndry
!
          ALLOCATE(PARENT_TASK(N)%WEST_V%STRING(1:LENGTH))                 !<-- Wboundary V datastring to be received from parent task
!
          INDX_MIN_V%WEST=MIN(INDX_MIN_V%WEST,PARENT_INFO(2,N))            !<-- Starting child J for union of parent task segments sent
          INDX_MAX_V%WEST=MAX(INDX_MAX_V%WEST,PARENT_INFO(3,N))            !<-- Ending child J for union of parent task segments sent
        ENDDO
      ENDIF
!
      NUM_PARENT_TASKS_SENDING_V%WEST=KOUNT
!!!   LENGTH_BND_SEG_V%WEST=INDX_MAX_V%WEST-INDX_MIN_V%WEST+1
!
      west_v: IF(NUM_PARENT_TASKS_SENDING_V%WEST>0)THEN                    !<-- Does this child task recv Wboundary V data from parent?
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
        ALLOCATE(BOUND_1D_WEST_V(1:LENGTH),stat=ISTAT)                     !<-- 1-D combined V-point data on child task's Wbndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Lower I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_LO_WEST_V'                   &  !<-- Name of the boundary array's lower I limit
                              ,value=ILIM_LO                            &  !<-- The boundary array's lower I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Upper I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_HI_WEST_V'                   &  !<-- Name of the boundary array's upper I limit
                              ,value=ILIM_HI                            &  !<-- The boundary array's upper I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Lower J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_LO_WEST_V'                   &  !<-- Name of the boundary array's lower J limit
                              ,value=JLIM_LO                            &  !<-- The boundary array's lower J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Upper J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_HI_WEST_V'                   &  !<-- Name of the boundary array's upper J limit
                              ,value=JLIM_HI                            &  !<-- The boundary array's upper J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF west_v
!
!-----------------------------------------------------------------------
!
!------------
!***  East H
!------------
!
      KOUNT=0
      INDX_MIN_H%EAST= 1000000
      INDX_MAX_H%EAST=-1000000
!
      DO N=1,FTASKS_DOMAIN(ID_DOM)                                         !<-- Child task loops through its parent's tasks
!
        CALL MPI_RECV(INFO                                              &  !<-- Receive data packet from each parent task
                     ,4                                                 &  !<-- # of words in data packet
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,MPI_ANY_SOURCE                                    &  !<-- Accept data from any parent task that is sending
                     ,77777                                             &  !<-- Tag used for east boundary H points
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator between child and parent
                     ,STATUS                                            &  !<-- Status of Recv
                     ,IERR)
!
        IF(INFO(1)>=0)THEN                                                 !<-- If yes, this parent task has key preliminary bndry info
          KOUNT=KOUNT+1
          DO N1=1,4
            PARENT_INFO(N1,KOUNT)=INFO(N1)                                 !<-- Save the data from that parent task
          ENDDO
        ENDIF
!
      ENDDO
!
      IF(KOUNT==2)THEN                                                     !<-- Nest task recvs data from two parent tasks
        IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                          !<-- Data recvd from 'out of order' parent tasks
!
          DO N1=1,4                                                        !<-- Save parent data in order of ascending task IDs
            TEMP(N1)         =PARENT_INFO(N1,1)                            !
            PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                            !
            PARENT_INFO(N1,2)=TEMP(N1)                                     !<--
          ENDDO
!
        ENDIF
      ENDIF
!
      IF(KOUNT>0)THEN
        DO N=1,KOUNT
          PARENT_TASK(N)%EAST_H%ID_SOURCE   =PARENT_INFO(1,N)              !<-- Rank of parent task that will send Eboundary H data segment
          PARENT_TASK(N)%EAST_H%INDX_START  =PARENT_INFO(2,N)              !<-- Jstart on child grid of the boundary data segment
          PARENT_TASK(N)%EAST_H%INDX_END    =PARENT_INFO(3,N)              !<-- Jend on child grid of the boundary data segment
          PARENT_TASK(N)%EAST_H%INDX_END_EXP=PARENT_INFO(4,N)              !<-- Jend on child grid of the expanded boundary data segment
!
          NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_H
          LENGTH=(3*LM+1)*NBASE
!
          PARENT_TASK(N)%EAST_H%LENGTH=LENGTH                              !<-- # of words in this parent task's datastring of child bndry
!
          ALLOCATE(PARENT_TASK(N)%EAST_H%STRING(1:LENGTH))                 !<-- Eboundary H datastring to be received from parent task
!
          INDX_MIN_H%EAST=MIN(INDX_MIN_H%EAST,PARENT_INFO(2,N))            !<-- Starting child J for union of parent task segments sent
          INDX_MAX_H%EAST=MAX(INDX_MAX_H%EAST,PARENT_INFO(3,N))            !<-- Ending child J for union of parent task segments sent
        ENDDO
      ENDIF
!
      NUM_PARENT_TASKS_SENDING_H%EAST=KOUNT
!!!   LENGTH_BND_SEG_H%EAST=INDX_MAX_H%EAST-INDX_MIN_H%EAST+1
!
      east_h: IF(NUM_PARENT_TASKS_SENDING_H%EAST>0)THEN                    !<-- Does this child task recv Eboundary H data from parent?
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
        ALLOCATE(BOUND_1D_EAST_H(1:LENGTH),stat=ISTAT)                     !<-- 1-D combined H-point data on child task's Ebndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Lower I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_LO_EAST_H'                   &  !<-- Name of the boundary array's lower I limit
                              ,value=ILIM_LO                            &  !<-- The boundary array's lower I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Upper I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_HI_EAST_H'                   &  !<-- Name of the boundary array's upper I limit
                              ,value=ILIM_HI                            &  !<-- The boundary array's upper I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Lower J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_LO_EAST_H'                   &  !<-- Name of the boundary array's lower J limit
                              ,value=JLIM_LO                            &  !<-- The boundary array's lower J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry H Data Upper J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_HI_EAST_H'                   &  !<-- Name of the boundary array's upper J limit
                              ,value=JLIM_HI                            &  !<-- The boundary array's upper J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF east_h
!
!-----------------------------------------------------------------------
!
!------------
!***  East V
!------------
!
      KOUNT=0
      INDX_MIN_V%EAST= 1000000
      INDX_MAX_V%EAST=-1000000
!
      DO N=1,FTASKS_DOMAIN(ID_DOM)                                         !<-- Child task loops through its parent's tasks
!
        CALL MPI_RECV(INFO                                              &  !<-- Receive data packet from each parent task
                     ,3                                                 &  !<-- # of words in data packet
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,MPI_ANY_SOURCE                                    &  !<-- Accept data from any parent task that is sending
                     ,88888                                             &  !<-- Tag used for east boundary V points
                     ,COMM_TO_MY_PARENT                                 &  !<-- MPI communicator between child and parent
                     ,STATUS                                            &  !<-- Status of Recv
                     ,IERR)
!
        IF(INFO(1)>=0)THEN                                                 !<-- If yes, this parent task has key preliminary bndry info
          KOUNT=KOUNT+1
          DO N1=1,3
            PARENT_INFO(N1,KOUNT)=INFO(N1)                                 !<-- Save the data from that parent task
          ENDDO
        ENDIF
!
      ENDDO
!
      IF(KOUNT==2)THEN                                                     !<-- Nest task recvs data from two parent tasks
        IF(PARENT_INFO(1,1)>PARENT_INFO(1,2))THEN                          !<-- Data recvd from 'out of order' parent tasks
!
          DO N1=1,3                                                        !<-- Save parent data in order of ascending task IDs
            TEMP(N1)         =PARENT_INFO(N1,1)                            !
            PARENT_INFO(N1,1)=PARENT_INFO(N1,2)                            !
            PARENT_INFO(N1,2)=TEMP(N1)                                     !<--
          ENDDO
!
        ENDIF
      ENDIF
!
      IF(KOUNT>0)THEN
        DO N=1,KOUNT
!
          PARENT_TASK(N)%EAST_V%ID_SOURCE =PARENT_INFO(1,N)                !<-- Rank of parent task that will send Eboundary V data segment
          PARENT_TASK(N)%EAST_V%INDX_START=PARENT_INFO(2,N)                !<-- Jstart on child grid of the boundary data segment
          PARENT_TASK(N)%EAST_V%INDX_END  =PARENT_INFO(3,N)                !<-- Jend on child grid of the boundary data segment
!
          NBASE =(PARENT_INFO(3,N)-PARENT_INFO(2,N)+1)*N_BLEND_V
          LENGTH=2*LM*NBASE
!
          PARENT_TASK(N)%EAST_V%LENGTH=LENGTH                              !<-- # of words in this parent task's datastring of child bndry
!
          ALLOCATE(PARENT_TASK(N)%EAST_V%STRING(1:LENGTH))                 !<-- Eboundary V datastring to be received from parent task
!
          INDX_MIN_V%EAST=MIN(INDX_MIN_V%EAST,PARENT_INFO(2,N))            !<-- Starting child J for union of parent task segments sent
          INDX_MAX_V%EAST=MAX(INDX_MAX_V%EAST,PARENT_INFO(3,N))            !<-- Ending child J for union of parent task segments sent
        ENDDO
      ENDIF
!
      NUM_PARENT_TASKS_SENDING_V%EAST=KOUNT
!!!   LENGTH_BND_SEG_V%EAST=INDX_MAX_V%EAST-INDX_MIN_V%EAST+1
!
      east_v: IF(NUM_PARENT_TASKS_SENDING_V%EAST>0)THEN                    !<-- Does this child task recv Eboundary V data from parent?
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
        ALLOCATE(BOUND_1D_EAST_V(1:LENGTH),stat=ISTAT)                     !<-- 1-D combined V-point data on child task's Ebndry segment
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Lower I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_LO_EAST_V'                   &  !<-- Name of the boundary array's lower I limit
                              ,value=ILIM_LO                            &  !<-- The boundary array's lower I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Upper I Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='ILIM_HI_EAST_V'                   &  !<-- Name of the boundary array's upper I limit
                              ,value=ILIM_HI                            &  !<-- The boundary array's upper I limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Lower J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_LO_EAST_V'                   &  !<-- Name of the boundary array's lower J limit
                              ,value=JLIM_LO                            &  !<-- The boundary array's lower J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Load Child Bndry V Data Upper J Limit"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state=EXP_STATE                          &  !<-- The Parent_Child Coupler export state
                              ,name ='JLIM_HI_EAST_V'                   &  !<-- Name of the boundary array's upper J limit
                              ,value=JLIM_HI                            &  !<-- The boundary array's upper J limit
                              ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_LIMITS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF east_v
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CHILD_RECVS_CHILD_DATA_LIMITS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE CHILD_SENDS_TOPO_TO_PARENT(IMP_STATE)
!
!-----------------------------------------------------------------------
!***  The children send their boundary surface geopotential to their 
!***  parents so the parents can properly balance their own data that 
!***  they interpolate to child boundary gridpoints. 
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument Variables
!-----------------------
!
      TYPE(ESMF_State),INTENT(IN) :: IMP_STATE                             !<-- Parent-Child Coupler import state
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IERR,KOUNT,MYPE                             &
                           ,N,N1,N2,NUM_WORDS                           &
                           ,RC,RC_FIS
!
      INTEGER(kind=KINT) :: I_LO,I_HI,I_OFFSET                          &
                           ,J_LO,J_HI,J_OFFSET
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: STATUS
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: FIS_SEND
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!     I_OFFSET=IMS-1+NHALO                                                 !<-- Offset of I in unloaded FIS vs. original FIS
!     J_OFFSET=JMS-1+NHALO                                                 !<-- Offset of J in unloaded FIS vs. original FIS
      I_OFFSET=0                                                           !<-- ESMF_INDEX now GLOBAL
      J_OFFSET=0                                                           !<-- ESMF_INDEX now GLOBAL
!
!-----------------------------------------------------------------------
!***  Extract the Sfc Geopotential from the Coupler's import state.
!***  If this child domain is also a parent then it already extracted
!***  its FIS in 'parent_block' of PARENT_CHILD_CPL_INITIALIZE but we
!***  now extract FIS again in case this child domain is not a parent.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract FIS Field from Parent-Child Coupler Import State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =IMP_STATE                             &  !<-- The parent-child coupler import state
                        ,itemName='FIS'                                 &  !<-- Extract FIS Field
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FIS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract FIS from ESMF Field in Parent-Child Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#ifdef ESMF_3
      CALL ESMF_FieldGet(field  =HOLD_FIELD                             &  !<-- Field that holds the data pointer
                        ,localDe=0                                      &
                        ,farray =FIS                                    &  !<-- Put the pointer here
                        ,rc     =RC)
#else
      CALL ESMF_FieldGet(field    =HOLD_FIELD                           &  !<-- Field that holds the data pointer
                        ,localDe  =0                                    &
                        ,farrayPtr=FIS                                  &  !<-- Put the pointer here
                        ,rc       =RC)
#endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FIS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_RANK(COMM_TO_MY_PARENT,MYPE,IERR)                      !<-- Obtain rank of this child task
!
!-----------------------------------------------------------------------
!
!------------------------------
!***  Child South Boundary FIS
!------------------------------
!
      IF(NUM_PARENT_TASKS_SENDING_H%SOUTH>0)THEN                           !<-- Child tasks know which parent tasks compute their BCs
!
        DO N=1,NUM_PARENT_TASKS_SENDING_H%SOUTH                            !<-- Child sends its FIS to parent tasks that will be
!                                                                               computing its BCs.
          I_LO=PARENT_TASK(N)%SOUTH_H%INDX_START                           !<-- Starting I of child covered by parent task N
          I_HI=PARENT_TASK(N)%SOUTH_H%INDX_END_EXP                         !<-- Ending I of child for expanded area covered by parent task N
          NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H+1)                            !<-- # of child points covered by parent task N
          ALLOCATE(FIS_SEND(1:NUM_WORDS))                                  !<-- Array to hold child FIS values to go to parent task N
                                                                           !    with extra row of values for 4-pt interp to V pts
          KOUNT=0
!
          DO N2=1,N_BLEND_H+1                                              !<-- Extra row for 4-pt interp of PD to V pts
          DO N1=I_LO,I_HI
            KOUNT=KOUNT+1
            FIS_SEND(KOUNT)=FIS(N1-I_OFFSET,N2-J_OFFSET)                   !<-- 1-D FIS of child points covered by parent task N
          ENDDO
          ENDDO
!
          CALL MPI_SEND(FIS_SEND                                        &  !<-- Send FIS data to parent task N
                       ,NUM_WORDS                                       &  !<-- There are NUM_WORDS words in the data
                       ,MPI_REAL                                        &  !<-- Data is type Real
                       ,PARENT_TASK(N)%SOUTH_H%ID_SOURCE                &  !<-- Data sent to this parent task
                       ,MYPE                                            &  !<-- Use child task ID as the tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between this child and its parent
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
      IF(NUM_PARENT_TASKS_SENDING_H%NORTH>0)THEN                           !<-- Child tasks know which parent tasks compute their BCs
!                                                                               
        DO N=1,NUM_PARENT_TASKS_SENDING_H%NORTH                            !<-- Child sends its FIS to parent tasks that will be
!                                                                               computing its BCs.
          I_LO=PARENT_TASK(N)%NORTH_H%INDX_START                           !<-- Starting I of child covered by parent task N
          I_HI=PARENT_TASK(N)%NORTH_H%INDX_END_EXP                         !<-- Ending I of child for expanded area covered by parent task N
          NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H+1)                            !<-- # of child points covered by parent task N
          ALLOCATE(FIS_SEND(1:NUM_WORDS))                                  !<-- Array to hold child FIS values to go to parent task N
                                                                           !    with extra row of values for 4-pt interp to V pts
          KOUNT=0
!
          DO N2=JTE-N_BLEND_H  ,JTE                                        !<-- Extra row for 4-pt interp of PD to V pts
          DO N1=I_LO,I_HI
            KOUNT=KOUNT+1
            FIS_SEND(KOUNT)=FIS(N1-I_OFFSET,N2-J_OFFSET)                   !<-- 1-D FIS of child points covered by parent task N
          ENDDO
          ENDDO
!
          CALL MPI_SEND(FIS_SEND                                        &  !<-- Send FIS data to parent task N
                       ,NUM_WORDS                                       &  !<-- There are NUM_WORDS words in the data
                       ,MPI_REAL                                        &  !<-- Data is type Real
                       ,PARENT_TASK(N)%NORTH_H%ID_SOURCE                &  !<-- Data sent to this parent task
                       ,MYPE                                            &  !<-- Use child task ID as the tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between this child and its parent
                       ,IERR)
!
          DEALLOCATE(FIS_SEND)
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
      IF(NUM_PARENT_TASKS_SENDING_H%WEST>0)THEN                            !<-- Child tasks know which parent tasks compute their BCs
!                                                                          
        DO N=1,NUM_PARENT_TASKS_SENDING_H%WEST                             !<-- Child sends its FIS to parent tasks that will be
!                                                                               computing its BCs.
          J_LO=PARENT_TASK(N)%WEST_H%INDX_START                            !<-- Starting J of child covered by parent task N
          J_HI=PARENT_TASK(N)%WEST_H%INDX_END_EXP                          !<-- Ending J of child for expanded area covered by parent task N
          NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H+1)                            !<-- # of child points covered by parent task N
          ALLOCATE(FIS_SEND(1:NUM_WORDS))                                  !<-- Array to hold child FIS values to go to parent task N
                                                                           !    with extra row of values for 4-pt interp to V pts
          KOUNT=0
!
          DO N2=J_LO,J_HI
          DO N1=1,N_BLEND_H+1                                              !<-- Extra row for 4-pt interp of PD to V pts
            KOUNT=KOUNT+1
            FIS_SEND(KOUNT)=FIS(N1-I_OFFSET,N2-J_OFFSET)                   !<-- 1-D FIS of child points covered by parent task N
          ENDDO
          ENDDO
!
          CALL MPI_SEND(FIS_SEND                                        &  !<-- Send FIS data to parent task N
                       ,NUM_WORDS                                       &  !<-- There are NUM_WORDS words in the data
                       ,MPI_REAL                                        &  !<-- Data is type Real
                       ,PARENT_TASK(N)%WEST_H%ID_SOURCE                 &  !<-- Data sent to this parent task
                       ,MYPE                                            &  !<-- Use child task ID as the tag
                       ,COMM_TO_MY_PARENT                               &  !<-- MPI communicator between this child and its parent
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
          NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H+1)                          !<-- # of child points covered by parent task N
          ALLOCATE(FIS_SEND(1:NUM_WORDS))                                !<-- Array to hold child FIS values to go to parent task N
                                                                         !    with extra row of values for 4-pt interp to V pts
          KOUNT=0
!
          DO N2=J_LO,J_HI
          DO N1=ITE-N_BLEND_H  ,ITE                                      !<-- Extra row for 4-pt interp of PD to V pts
            KOUNT=KOUNT+1
            FIS_SEND(KOUNT)=FIS(N1-I_OFFSET,N2-J_OFFSET)                 !<-- 1-D FIS of child points covered by parent task N
          ENDDO
          ENDDO
!
          CALL MPI_SEND(FIS_SEND                                      &  !<-- Send FIS data to parent task N
                       ,NUM_WORDS                                     &  !<-- There are NUM_WORDS words in the data
                       ,MPI_REAL                                      &  !<-- Data is type Real
                       ,PARENT_TASK(N)%EAST_H%ID_SOURCE               &  !<-- Data sent to this parent task
                       ,MYPE                                          &  !<-- Use child task ID as the tag
                       ,COMM_TO_MY_PARENT                             &  !<-- MPI communicator between this child and its parent
                       ,IERR)
!
          DEALLOCATE(FIS_SEND)
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CHILD_SENDS_TOPO_TO_PARENT
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_RECVS_CHILD_TOPO(N_CHILD)
!
!-----------------------------------------------------------------------
!***  Parents receive the boundary topography from their children
!***  so the parents can properly balance their own data that they
!***  interpolate to child gridpoints.
!-----------------------------------------------------------------------
!
!-----------------------
!*** Argument Variables
!-----------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: N_CHILD                              !<-- The child who is sending
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IERR,N,NTX,NUM_WORDS
!
      INTEGER(kind=KINT) :: I_LO,I_HI,J_LO,J_HI
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: STATUS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      N=N_CHILD
!
!-----------------------------------------------------------------------
!***  The parent tasks receive sfc geopotentials from the
!***  child boundary tasks.
!-----------------------------------------------------------------------
!
!------------------------------
!***  Child South Boundary FIS
!------------------------------
!
      IF(NUM_TASKS_SEND_H_S(N)>0)THEN                                      !<-- This parent task covers some child Sboundary H points
!
        ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(1:NUM_TASKS_SEND_H_S(N)))        !<-- FIS data slot for each Sbndry child task on parent task
!
        DO NTX=1,NUM_TASKS_SEND_H_S(N)                                     !<-- Loop through those particular child tasks 
!                                                                          
          I_LO=CHILDTASK_H_SAVE(N)%I_LO_SOUTH(NTX)                         !<-- Starting I of child point bndry segment on this parent task
          I_HI=CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NTX)                         !<-- Ending I of child point bndry segment on this parent task
          NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H_CHILD(N)+1)                   !<-- # of child points in its bndry segment on parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
          ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(NTX)%DATA(1:NUM_WORDS))        !<-- FIS on child S bndry covered by this parent task
!
          CALL MPI_RECV(FIS_CHILD_SOUTH(N)%TASKS(NTX)%DATA              &  !<-- Recv FIS values from Sbndry child task NTX
                       ,NUM_WORDS                                       &  !<-- # of FIS values
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NTX)           &  !<-- The child task sending
                       ,CHILDTASK_BNDRY_H_RANKS(N)%SOUTH(NTX)           &  !<-- Tag
                       ,COMM_TO_MY_CHILDREN(N)                          &  !<-- MPI communicator
                       ,STATUS                                          &
                       ,IERR)
!
        ENDDO
!
      ELSE
!
        ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(1:1))
        ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(1)%DATA(1:1))
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!------------------------------
!***  Child North Boundary FIS
!------------------------------
!
      IF(NUM_TASKS_SEND_H_N(N)>0)THEN                                      !<-- This parent task covers some child Nboundary H points
!
        ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(1:NUM_TASKS_SEND_H_N(N)))        !<-- FIS data slot for each Nbndry child task on parent task
!
        DO NTX=1,NUM_TASKS_SEND_H_N(N)                                     !<-- Loop through those particular child tasks 
!                                                                          
          I_LO=CHILDTASK_H_SAVE(N)%I_LO_NORTH(NTX)                         !<-- Starting I of child point bndry segment on this parent task
          I_HI=CHILDTASK_H_SAVE(N)%I_HI_NORTH(NTX)                         !<-- Ending I of child point bndry segment on this parent task
          NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H_CHILD(N)+1)                   !<-- # of child points in its bndry segment on parent task
          ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(NTX)%DATA(1:NUM_WORDS))        !<-- FIS on child N bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
!
          CALL MPI_RECV(FIS_CHILD_NORTH(N)%TASKS(NTX)%DATA              &  !<-- Recv FIS values from Nbndry child task NTX
                       ,NUM_WORDS                                       &  !<-- # of FIS values
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NTX)           &  !<-- The child task sending
                       ,CHILDTASK_BNDRY_H_RANKS(N)%NORTH(NTX)           &  !<-- Tag
                       ,COMM_TO_MY_CHILDREN(N)                          &  !<-- MPI communicator
                       ,STATUS                                          &
                       ,IERR)
!
        ENDDO
!
      ELSE
!
        ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(1:1))
        ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(1)%DATA(1:1))
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!-----------------------------
!***  Child West Boundary FIS
!-----------------------------
!
      IF(NUM_TASKS_SEND_H_W(N)>0)THEN                                      !<-- This parent task covers some child Wboundary H points
!
        ALLOCATE(FIS_CHILD_WEST(N)%TASKS(1:NUM_TASKS_SEND_H_W(N)))         !<-- FIS data slot for each Wbndry child task on parent task
!
        DO NTX=1,NUM_TASKS_SEND_H_W(N)                                     !<-- Loop through those particular child tasks 
!                                                                          
          J_LO=CHILDTASK_H_SAVE(N)%J_LO_WEST(NTX)                          !<-- Starting J of child point bndry segment on this parent task
          J_HI=CHILDTASK_H_SAVE(N)%J_HI_WEST(NTX)                          !<-- Ending J of child point bndry segment on this parent task
          NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H_CHILD(N)+1)                   !<-- # of child points in its bndry segment on parent task
          ALLOCATE(FIS_CHILD_WEST(N)%TASKS(NTX)%DATA(1:NUM_WORDS))         !<-- FIS on child W bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
!
          CALL MPI_RECV(FIS_CHILD_WEST(N)%TASKS(NTX)%DATA               &  !<-- Recv FIS values from Wbndry child task NTX
                       ,NUM_WORDS                                       &  !<-- # of FIS values
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,CHILDTASK_BNDRY_H_RANKS(N)%WEST(NTX)            &  !<-- The child task sending
                       ,CHILDTASK_BNDRY_H_RANKS(N)%WEST(NTX)            &  !<-- Tag
                       ,COMM_TO_MY_CHILDREN(N)                          &  !<-- MPI communicator
                       ,STATUS                                          &
                       ,IERR)
!
        ENDDO
!
      ELSE
!
        ALLOCATE(FIS_CHILD_WEST(N)%TASKS(1:1))
        ALLOCATE(FIS_CHILD_WEST(N)%TASKS(1)%DATA(1:1))
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!-----------------------------
!***  Child East Boundary FIS
!-----------------------------
!
      IF(NUM_TASKS_SEND_H_E(N)>0)THEN                                      !<-- This parent task covers some child Eboundary H points
!
        ALLOCATE(FIS_CHILD_EAST(N)%TASKS(1:NUM_TASKS_SEND_H_E(N)))         !<-- FIS data slot for each Ebndry child task on parent task
!
        DO NTX=1,NUM_TASKS_SEND_H_E(N)                                     !<-- Loop through those particular child tasks 
!                                                                          
          J_LO=CHILDTASK_H_SAVE(N)%J_LO_EAST(NTX)                          !<-- Starting J of child point bndry segment on this parent task
          J_HI=CHILDTASK_H_SAVE(N)%J_HI_EAST(NTX)                          !<-- Ending J of child point bndry segment on this parent task
          NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H_CHILD(N)+1)                   !<-- # of child points in its bndry segment on parent task
          ALLOCATE(FIS_CHILD_EAST(N)%TASKS(NTX)%DATA(1:NUM_WORDS))         !<-- FIS on child E bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
!
          CALL MPI_RECV(FIS_CHILD_EAST(N)%TASKS(NTX)%DATA               &  !<-- Recv FIS values from Ebndry child task NTX
                       ,NUM_WORDS                                       &  !<-- # of FIS values
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,CHILDTASK_BNDRY_H_RANKS(N)%EAST(NTX)            &  !<-- The child task sending
                       ,CHILDTASK_BNDRY_H_RANKS(N)%EAST(NTX)            &  !<-- Tag
                       ,COMM_TO_MY_CHILDREN(N)                          &  !<-- MPI communicator
                       ,STATUS                                          &
                       ,IERR)
!
        ENDDO
!
      ELSE
!
        ALLOCATE(FIS_CHILD_EAST(N)%TASKS(1:1))
        ALLOCATE(FIS_CHILD_EAST(N)%TASKS(1)%DATA(1:1))
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_RECVS_CHILD_TOPO
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_COMPUTES_CHILD_TOPO(N_CHILD                     &
                                           ,I_PARENT_SW                 &
                                           ,J_PARENT_SW                 &
                                           ,IM_CHILD                    &
                                           ,JM_CHILD                    &
                                           ,N_BLEND_H_CHILD             &
                                           ,LBND1,UBND1,LBND2,UBND2     &
                                           ,MOVING_NEST_TOPO            &
                                             )
!
!-----------------------------------------------------------------------
!***  Parents fill the working arrays of their moving children's
!***  boundary topography.  The parents carry full arrays of topography
!***  at each of their moving children's resolutions so the data only
!***  needs to be lifted from those arrays.  This avoids Sends and
!***  Recvs of that data which could force the parents to Wait. 
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: N_CHILD                          &  !<-- Which child are we considering?
                                      ,I_PARENT_SW                      &  !<-- Parent I index of child's SW corner 
                                      ,J_PARENT_SW                      &  !<-- Parent J index of child's SW corner 
                                      ,IM_CHILD                         &  !<-- I dimension of child domain
                                      ,JM_CHILD                         &  !<-- J dimension of child domain
                                      ,N_BLEND_H_CHILD                  &  !<-- Width of child's boundary blending region
!
                                      ,LBND1,LBND2                      &  !<-- I limits of nest-res FIS on this parentt ask
                                      ,UBND1,UBND2                     
!
      REAL(kind=KFPT),DIMENSION(LBND1:UBND1,LBND2:UBND2),INTENT(IN) ::  &
                                                      MOVING_NEST_TOPO     !<-- Nest-resolution topography on the parent task
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,I_HI,I_LO,I_OFFSET,ISTART                 &
                           ,J,J_HI,J_LO,J_OFFSET,JSTART                 &
                           ,KOUNT,N,NTX,NUM_WORDS
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: FIS_X
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      N=N_CHILD
!
      ISTART=MAX(IMS,IDS)                                                  !<-- The SW corner of this parent domain
      JSTART=MAX(JMS,JDS)                                                  !<--
!
!------------------------------
!***  Child South Boundary FIS
!------------------------------
!
      IF(NUM_TASKS_SEND_H_S(N)>0)THEN                                      !<-- This parent task covers some child Sboundary H points
!
        ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(1:NUM_TASKS_SEND_H_S(N)))        !<-- FIS data slot for each Sbndry child task on parent task
!
        DO NTX=1,NUM_TASKS_SEND_H_S(N)                                     !<-- Loop through those particular child tasks 
!                                                                          
          I_LO=CHILDTASK_H_SAVE(N)%I_LO_SOUTH(NTX)                         !<-- Starting I of child point bndry segment on this parent task
          I_HI=CHILDTASK_H_SAVE(N)%I_HI_SOUTH(NTX)                         !<-- Ending I of child point bndry segment on this parent task
          NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H_CHILD+1)                      !<-- # of child points in its bndry segment on parent task
!
          ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(NTX)%DATA(1:NUM_WORDS))        !<-- FIS on child S bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
          FIS_X=>FIS_CHILD_SOUTH(N)%TASKS(NTX)%DATA
!
          I_OFFSET=(I_PARENT_SW-ISTART)*PARENT_CHILD_SPACE_RATIO(N)     &  !<-- I offset of child SW corner in full topo array on parent
                   +LBND1-1                                         
          J_OFFSET=(J_PARENT_SW-JSTART)*PARENT_CHILD_SPACE_RATIO(N)     &  !<-- J offset of child SW corner in full topo array on parent
                   +LBND2-1
          KOUNT=0
!
          DO J=1,N_BLEND_H_CHILD+1
          DO I=I_LO,I_HI
            KOUNT=KOUNT+1
            FIS_X(KOUNT)=MOVING_NEST_TOPO(I+I_OFFSET,J+J_OFFSET)           !<-- Lift child topography from full parent array
          ENDDO
          ENDDO
!
        ENDDO
!
      ELSE
!
        ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(1:1))
        ALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(1)%DATA(1:1))
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!------------------------------
!***  Child North Boundary FIS
!------------------------------
!
      IF(NUM_TASKS_SEND_H_N(N)>0)THEN                                      !<-- This parent task covers some child Nboundary H points
!
        ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(1:NUM_TASKS_SEND_H_N(N)))        !<-- FIS data slot for each Nbndry child task on parent task
!
        DO NTX=1,NUM_TASKS_SEND_H_N(N)                                     !<-- Loop through those particular child tasks 
!                                                                          
          I_LO=CHILDTASK_H_SAVE(N)%I_LO_NORTH(NTX)                         !<-- Starting I of child point bndry segment on this parent task
          I_HI=CHILDTASK_H_SAVE(N)%I_HI_NORTH(NTX)                         !<-- Ending I of child point bndry segment on this parent task
          NUM_WORDS=(I_HI-I_LO+1)*(N_BLEND_H_CHILD+1)                      !<-- # of child points in its bndry segment on parent task
!
          ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(NTX)%DATA(1:NUM_WORDS))        !<-- FIS on child N bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
          FIS_X=>FIS_CHILD_NORTH(N)%TASKS(NTX)%DATA
!
          I_OFFSET=(I_PARENT_SW-ISTART)*PARENT_CHILD_SPACE_RATIO(N)     &  !<-- I offset of child SW corner in full topo array on parent
                   +LBND1-1                                          
          J_OFFSET=(J_PARENT_SW-JSTART)*PARENT_CHILD_SPACE_RATIO(N)     &  !<-- J offset of child SW corner in full topo array on parent
                   +LBND2-1                                          
          KOUNT=0
!
!!!       DO J=JM_CHILD-N_BLEND_H_CHILD+1,JM_CHILD
          DO J=JM_CHILD-N_BLEND_H_CHILD  ,JM_CHILD
          DO I=I_LO,I_HI
            KOUNT=KOUNT+1
            FIS_X(KOUNT)=MOVING_NEST_TOPO(I+I_OFFSET,J+J_OFFSET)           !<-- Lift child topography from full parent array
          ENDDO
          ENDDO
!
        ENDDO
!
      ELSE
!
        ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(1:1))
        ALLOCATE(FIS_CHILD_NORTH(N)%TASKS(1)%DATA(1:1))
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!-----------------------------
!***  Child West Boundary FIS
!-----------------------------
!
      IF(NUM_TASKS_SEND_H_W(N)>0)THEN                                      !<-- This parent task covers some child Wboundary H points
!
        ALLOCATE(FIS_CHILD_WEST(N)%TASKS(1:NUM_TASKS_SEND_H_W(N)))         !<-- FIS data slot for each Wbndry child task on parent task
!
        DO NTX=1,NUM_TASKS_SEND_H_W(N)                                     !<-- Loop through those particular child tasks 
!                                                                          
          J_LO=CHILDTASK_H_SAVE(N)%J_LO_WEST(NTX)                          !<-- Starting J of child point bndry segment on this parent task
          J_HI=CHILDTASK_H_SAVE(N)%J_HI_WEST(NTX)                          !<-- Ending J of child point bndry segment on this parent task
          NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H_CHILD+1)                      !<-- # of child points in its bndry segment on parent task
!
          ALLOCATE(FIS_CHILD_WEST(N)%TASKS(NTX)%DATA(1:NUM_WORDS))         !<-- FIS on child W bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
          FIS_X=>FIS_CHILD_WEST(N)%TASKS(NTX)%DATA
!
          I_OFFSET=(I_PARENT_SW-ISTART)*PARENT_CHILD_SPACE_RATIO(N)     &  !<-- I offset of child SW corner in full topo array on parent
                   +LBND1-1                                          
          J_OFFSET=(J_PARENT_SW-JSTART)*PARENT_CHILD_SPACE_RATIO(N)     &  !<-- J offset of child SW corner in full topo array on parent
                   +LBND2-1                                          
          KOUNT=0
!
          DO J=J_LO,J_HI
          DO I=1,N_BLEND_H_CHILD+1
            KOUNT=KOUNT+1
            FIS_X(KOUNT)=MOVING_NEST_TOPO(I+I_OFFSET,J+J_OFFSET)           !<-- Lift child topography from full parent array
          ENDDO
          ENDDO
!
        ENDDO
!
      ELSE
!
        ALLOCATE(FIS_CHILD_WEST(N)%TASKS(1:1))
        ALLOCATE(FIS_CHILD_WEST(N)%TASKS(1)%DATA(1:1))
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!-----------------------------
!***  Child East Boundary FIS
!-----------------------------
!
      IF(NUM_TASKS_SEND_H_E(N)>0)THEN                                      !<-- This parent task covers some child Eboundary H points
!
        ALLOCATE(FIS_CHILD_EAST(N)%TASKS(1:NUM_TASKS_SEND_H_E(N)))         !<-- FIS data slot for each Ebndry child task on parent task
!
        DO NTX=1,NUM_TASKS_SEND_H_E(N)                                     !<-- Loop through those particular child tasks 
!                                                                          
          J_LO=CHILDTASK_H_SAVE(N)%J_LO_EAST(NTX)                          !<-- Starting J of child point bndry segment on this parent task
          J_HI=CHILDTASK_H_SAVE(N)%J_HI_EAST(NTX)                          !<-- Ending J of child point bndry segment on this parent task
          NUM_WORDS=(J_HI-J_LO+1)*(N_BLEND_H_CHILD+1)                      !<-- # of child points in its bndry segment on parent task
!
          ALLOCATE(FIS_CHILD_EAST(N)%TASKS(NTX)%DATA(1:NUM_WORDS))         !<-- FIS on child E bndry covered by this parent task
                                                                           !    with extra row for 4-pt interpolation of PD to V pts
          FIS_X=>FIS_CHILD_EAST(N)%TASKS(NTX)%DATA
!
          I_OFFSET=(I_PARENT_SW-ISTART)*PARENT_CHILD_SPACE_RATIO(N)     &  !<-- I offset of child SW corner in full topo array on parent
                   +LBND1-1                                          
          J_OFFSET=(J_PARENT_SW-JSTART)*PARENT_CHILD_SPACE_RATIO(N)     &  !<-- J offset of child SW corner in full topo array on parent 
                   +LBND2-1                                          
          KOUNT=0
!
          DO J=J_LO,J_HI
!!!       DO I=IM_CHILD-N_BLEND_H_CHILD+1,IM_CHILD
          DO I=IM_CHILD-N_BLEND_H_CHILD  ,IM_CHILD
            KOUNT=KOUNT+1
            FIS_X(KOUNT)=MOVING_NEST_TOPO(I+I_OFFSET,J+J_OFFSET)           !<-- Lift child topography from full parent array
          ENDDO
          ENDDO
!
        ENDDO
!
      ELSE
!
        ALLOCATE(FIS_CHILD_EAST(N)%TASKS(1:1))
        ALLOCATE(FIS_CHILD_EAST(N)%TASKS(1)%DATA(1:1))
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_COMPUTES_CHILD_TOPO
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
                                         ,JM_CHILD_X                    & !    Input
!                                                                           -----------
                                         ,CHILD_H_SBND                  & !    Output
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
!***  Given a child's actual surface geopotential, generate a new 
!***  value of PD for the child boundary points based on the 
!***  surrounding parent points.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: IMS,IME,JMS,JME                  &  !<-- Parent task's memory limits
                                      ,IM_CHILD_X,JM_CHILD_X            &  !<-- Index limits of the nest domain
                                      ,N_BLEND                          &  !<-- # of domain boundary blending rows 
                                      ,NLEV                             &  !<-- # of vertical levels in parent array
                                      ,NUM_CHILD_TASKS                  &  !<-- # of fcst tasks on this child
                                      ,NUM_TASKS_SEND_SBND              &  !<-- # of child tasks with Sboundary regions on parent task
                                      ,NUM_TASKS_SEND_NBND              &  !<-- # of child tasks with Nboundary regions on parent task
                                      ,NUM_TASKS_SEND_WBND              &  !<-- # of child tasks with Wboundary regions on parent task
                                      ,NUM_TASKS_SEND_EBND                 !<-- # of child tasks with Eboundary regions on parent task
!
!!!   INTEGER(kind=KINT),DIMENSION(1:NUM_CHILD_TASKS),INTENT(IN) ::     &
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(IN) ::             &
                                                 LOCAL_TASK_RANKS_S     &  !<-- Ranks of this child's Sbndry fcst tasks
                                                ,LOCAL_TASK_RANKS_N     &  !<-- Ranks of this child's Nbndry fcst tasks
                                                ,LOCAL_TASK_RANKS_W     &  !<-- Ranks of this child's Wbndry fcst tasks
                                                ,LOCAL_TASK_RANKS_E        !<-- Ranks of this child's Ebndry fcst tasks
!
      INTEGER(kind=KINT),DIMENSION(1:4,1:NUM_CHILD_TASKS),INTENT(IN) :: &
                                                                 LIMITS    !<-- ITS,ITE,JTS,JTE on each task of the child
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(IN) ::             &
                                                   I_LO_SOUTH           &  !<-- Starting I of Sbndry region on child tasks
                                                  ,I_HI_SOUTH           &  !<-- Ending I of Sbndry region on child tasks
                                                  ,I_HI_SOUTH_TRANSFER  &  !<-- Ending I of Sbndry region for transfer to child
                                                  ,I_LO_NORTH           &  !<-- Starting I of Nbndry region on child tasks
                                                  ,I_HI_NORTH           &  !<-- Ending I of Nbndry region on child tasks
                                                  ,I_HI_NORTH_TRANSFER  &  !<-- Ending I of Nbndry region for transfer to child
                                                  ,J_LO_WEST            &  !<-- Starting J of Wbndry region on child tasks
                                                  ,J_HI_WEST            &  !<-- Ending J of Wbndry region on child tasks
                                                  ,J_HI_WEST_TRANSFER   &  !<-- Ending J of Wbndry region for transfer to child
                                                  ,J_LO_EAST            &  !<-- Starting J of Ebndry region on child tasks
                                                  ,J_HI_EAST            &  !<-- Ending J of Ebndry region on child tasks
                                                  ,J_HI_EAST_TRANSFER      !<-- Ending J of Ebndry region for transfer to child
!
      INTEGER(kind=KINT),DIMENSION(:,:,:),POINTER,INTENT(IN) ::         &
                                                   I_INDX_PARENT_SBND   &  !<-- I indices of parent points west/east of child Sbndry point
                                                  ,I_INDX_PARENT_NBND   &  !<-- I indices of parent points west/east of child Nbndry point
                                                  ,I_INDX_PARENT_WBND   &  !<-- I indices of parent points west/east of child Wbndry point
                                                  ,I_INDX_PARENT_EBND   &  !<-- I indices of parent points west/east of child Ebndry point
                                                  ,J_INDX_PARENT_SBND   &  !<-- J indices of parent points south/north of child Sbndry point
                                                  ,J_INDX_PARENT_NBND   &  !<-- J indices of parent points south/north of child Nbndry point
                                                  ,J_INDX_PARENT_WBND   &  !<-- J indices of parent points south/north of child Wbndry point
                                                  ,J_INDX_PARENT_EBND      !<-- J indices of parent points south/north of child Ebndry point
!
      REAL(kind=KFPT),INTENT(IN) :: PT                                  &  !<-- Top pressure of model domain (Pa)
                                   ,PDTOP                                  !<-- Pressure at top of sigma domain (Pa)
!
      REAL(kind=KFPT),DIMENSION(:),POINTER,INTENT(IN) :: SG1            &  !<-- Interface sigmas, pressure domain
                                                        ,SG2               !<-- Interface sigmas, sigma domain
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS      &  !<-- Parent FIS
                                                              ,PD          !<-- Parent PD
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(IN) :: T &  !<-- Parent sensible temperature (K)
                                                                     ,Q    !<-- Parent specific humidity (kg/kg)
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER,INTENT(IN) ::            &
                                                          WEIGHT_SBND   &  !<-- Bilinear interp weights of the 4 parent points around 
                                                         ,WEIGHT_NBND   &  !    each point on child's boundary sides (S,N,W,E).
                                                         ,WEIGHT_WBND   &  !    
                                                         ,WEIGHT_EBND      !<--
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(IN) :: FIS_CHILD_SBND &  !<-- Sfc geopot on Sbndry points of each child task
                                                        ,FIS_CHILD_NBND &  !<-- Sfc geopot on Nbndry points of each child task
                                                        ,FIS_CHILD_WBND &  !<-- Sfc geopot on Wbndry points of each child task
                                                        ,FIS_CHILD_EBND    !<-- Sfc geopot on Ebndry points of each child task
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
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J,L                                       &
                           ,I_EAST,I_WEST,J_SOUTH,J_NORTH               &
                           ,I_START,I_END,J_START,J_END                 &
                           ,I_START_TRANSFER,I_END_TRANSFER             &
                           ,J_START_TRANSFER,J_END_TRANSFER             &
                           ,KOUNT_PTS                                   &
                           ,KOUNT_TRANSFER                              &
                           ,N_SIDE,NUM_TASKS_SEND,NTX                   &
                           ,RC
!
      INTEGER(kind=KINT),DIMENSION(:,:,:),POINTER :: I_INDX_PARENT_BND  &
                                                    ,J_INDX_PARENT_BND
!
      REAL(kind=KFPT) :: COEFF_1,COEFF_2,D_LNP_DFI,FIS_CHILD            &
                        ,LOG_P1_PARENT,PDTOP_PT,PHI_DIFF,PSFC_CHILD     &
                        ,Q_INTERP,T_INTERP
!
      REAL(kind=KFPT) :: PX_NE,PX_NW,PX_SE,PX_SW                        &
                        ,WGHT_NE,WGHT_NW,WGHT_SE,WGHT_SW
!
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE :: LOG_PBOT            &
                                                   ,LOG_PTOP
!
      REAL(kind=KFPT),DIMENSION(:,:,:),ALLOCATABLE :: PHI_INTERP        &
                                                     ,PINT_INTERP
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: WEIGHT_BND
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER :: CHILD_BOUND_H             &
                                             ,FIS_CHILD_BND             &
                                             ,PDB
!
      INTEGER(kind=KINT) :: LOG_LENGTH
!
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE :: TMP
!
      integer,dimension(8) :: values
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Loop through the four sides of the nest domain boundary (S,N,W,E).
!***  We use some dummy variables/pointers generically for all four
!***  of the sides.
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
!.......................................................................
!$omp parallel do private(                                              &
!$omp         coeff_1,coeff_2,d_lnp_dfi,fis_child,                      &
!$omp         i,i_east,i_end,i_end_transfer,                            &
!$omp         i_start,i_start_transfer,i_west,                          &
!$omp         j,j_end,j_end_transfer,j_north,j_south,                   &
!$omp         j_start,j_start_transfer,                                 &
!$omp         kount_pts,kount_transfer,                                 &
!$omp         l,log_length,log_p1_parent,log_pbot,log_ptop,             &
!$omp         ntx,pdtop_pt,phi_diff,phi_interp,pint_interp,psfc_child,  &
!$omp         px_ne,px_nw,px_se,px_sw,                                  &
!$omp         q_interp,t_interp,tmp,wght_ne,wght_nw,wght_se,wght_sw)
!.......................................................................
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
!***  Allocate the nest working arrays valid for the current child task
!***  on the child's grid.
!-----------------------------------------------------------------------
!
          ALLOCATE(PINT_INTERP(I_START:I_END,J_START:J_END,1:LM+1))
          ALLOCATE( PHI_INTERP(I_START:I_END,J_START:J_END,1:LM+1))
          ALLOCATE(   LOG_PTOP(I_START:I_END,J_START:J_END))
          ALLOCATE(   LOG_PBOT(I_START:I_END,J_START:J_END))
!
#ifdef IBM
          ALLOCATE(TMP(I_START:I_END,J_START:J_END))
          LOG_LENGTH=(I_END-I_START+1)*(J_END-J_START+1)
#endif
!
!-----------------------------------------------------------------------
!***  Compute parent heights of layer interfaces at the four points
!***  surrounding each child boundary point.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First the bottom layer (L=NLEV).
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
            PX_SW=PD(I_WEST,J_SOUTH)+PT                                    !<-- Sfc pressure on parent point SW of nest point
            PX_SE=PD(I_EAST,J_SOUTH)+PT                                    !<-- Sfc pressure on parent point SE of nest point
            PX_NW=PD(I_WEST,J_NORTH)+PT                                    !<-- Sfc pressure on parent point NW of nest point
            PX_NE=PD(I_EAST,J_NORTH)+PT                                    !<-- Sfc pressure on parent point NE of nest point
!
            PINT_INTERP(I,J,LM+1)=WGHT_SW*PX_SW                         &  !<-- Parent's surface pressure interp'd to this child's 
                                 +WGHT_SE*PX_SE                         &  !    gridpoint (I,J) along child's boundary 
                                 +WGHT_NW*PX_NW                         &  !    on child task NTX.
                                 +WGHT_NE*PX_NE               
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
!***  Now that we have the parent's sfc pressure and sfc geopotential
!***  at the child boundary points (I,J), compute the interface heights
!***  based on the horizontally interpolated interface pressure and
!***  the T and Q.
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
              PX_SW=SG2(L)*PD(I_WEST,J_SOUTH)+PDTOP_PT                     !<-- Pressure, top of layer L, parent point SW of nest point
              PX_SE=SG2(L)*PD(I_EAST,J_SOUTH)+PDTOP_PT                     !<-- Pressure, top of layer L, parent point SE of nest point
              PX_NW=SG2(L)*PD(I_WEST,J_NORTH)+PDTOP_PT                     !<-- Pressure, top of layer L, parent point NW of nest point
              PX_NE=SG2(L)*PD(I_EAST,J_NORTH)+PDTOP_PT                     !<-- Pressure, top of layer L, parent point NE of nest point
!
              WGHT_SW=WEIGHT_BND(I,J,INDX_SW)                              !<-- Bilinear weight for parent's point SW of child's point
              WGHT_SE=WEIGHT_BND(I,J,INDX_SE)                              !<-- Bilinear weight for parent's point SE of child's point
              WGHT_NW=WEIGHT_BND(I,J,INDX_NW)                              !<-- Bilinear weight for parent's point NW of child's point
              WGHT_NE=WEIGHT_BND(I,J,INDX_NE)                              !<-- Bilinear weight for parent's point NE of child's point
!
              PINT_INTERP(I,J,L)=WGHT_SW*PX_SW                          &  !<-- Top interface pressure interp'd to child gridpoint
                                +WGHT_SE*PX_SE                          &  !    along child's boundary for child task NTX
                                +WGHT_NW*PX_NW                          &
                                +WGHT_NE*PX_NE              
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
#ifndef IBM
              LOG_PTOP(I,J)=LOG(PINT_INTERP(I,J,L))                        !<-- Log of parent (top) interface pressure at child bndry point
!
              PHI_INTERP(I,J,L)=PHI_INTERP(I,J,L+1)                     &  !<-- Top interface geopotl of parent at child gridpoint (I,J)
                               +R_D*T_INTERP*(1.+P608*Q_INTERP)         &
                               *(LOG_PBOT(I,J)-LOG_PTOP(I,J))
!
              LOG_PBOT(I,J)=LOG_PTOP(I,J)                                  !<--- Move Log(Ptop) to bottom of next model layer up
!
#else
              TMP(I,J)=R_D*T_INTERP*(1.+P608*Q_INTERP)
#endif
            ENDDO
            ENDDO
!
#ifdef IBM
            CALL VSLOG(LOG_PTOP,PINT_INTERP(:,:,L),LOG_LENGTH)             !<-- Log of parent (top) interface pressure at child bndry point
            DO J=J_START,J_END                                             !<-- J limits of child task bndry region on parent task
            DO I=I_START,I_END                                             !<-- I limits of child task bndry region on parent task
              PHI_INTERP(I,J,L)=PHI_INTERP(I,J,L+1)                     &  !<-- Top interface geopotl of parent at child gridpoint (I,J)
                               +TMP(I,J)*(LOG_PBOT(I,J)-LOG_PTOP(I,J))
              LOG_PBOT(I,J)=LOG_PTOP(I,J)                                  !<--- Move Log(Ptop) to bottom of next model layer up
            ENDDO
            ENDDO
#endif
!
          ENDDO
!
!-----------------------------------------------------------------------
!***  Use the child's actual sfc geopotential to derive the value of
!***  PD at the child boundary points based on the parent's heights
!***  and pressures at its (the parent's) layer interfaces over the
!***  child's boundary points.
!
!***  If the child's terrain is lower than the value of the parent's
!***  terrain interpolated to the child point then extrapolate the
!***  parent's interpolated Sfc Pressure down to the child's terrain
!***  quadratically.
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
!***  Save only the PD's that are on the nests' boundary region's
!***  mass points into the datastring that will be transferred from
!***  the parent to the nests.  The extra points in the PDB pointer
!***  are for 4-PT averaging onto V points.
!***  Recall that the CHILD_BOUND_* pointer is the 1-D string of all
!***  data sent from parent tasks to child boundary tasks and PDB is
!***  the first data in that string.
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
#ifdef IBM
          DEALLOCATE(TMP)
#endif
!
!-----------------------------------------------------------------------
!
        ENDDO child_task_loop
!
!.......................................................................
!$omp end parallel do
!.......................................................................
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
                                          ,VBL_FLAG                     &
!
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
                                          ,JM_CHILD_X                   & !   Input
!                                                                           ----------
                                          ,VBL_CHILD_SBND               & !   Output
                                          ,VBL_CHILD_NBND               & !     |
                                          ,VBL_CHILD_WBND               & !     |
                                          ,VBL_CHILD_EBND )               !     v
!
!-----------------------------------------------------------------------
!***  Parent tasks interpolate their values of variables
!***  to child grid points. 
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: IMS,IME,JMS,JME                  &  !<-- Parent task's memory limits
                                      ,IM_CHILD_X,JM_CHILD_X            &  !<-- Index limits of the nest domain
                                      ,N_BLEND                          &  !<-- # of domain boundary blending rows
                                      ,NLEV                             &  !<-- # of vertical levels in parent array
                                      ,N_REMOVE                         &  !<-- # of rows to ignore on north/east sides (H=>0;V=>1)
!
                                      ,NUM_TASKS_SEND_SBND              &  !<-- # of child tasks with Sbndry points on parent task
                                      ,NUM_TASKS_SEND_NBND              &  !<-- # of child tasks with Nbndry points on parent task
                                      ,NUM_TASKS_SEND_WBND              &  !<-- # of child tasks with Wbndry points on parent task
                                      ,NUM_TASKS_SEND_EBND                 !<-- # of child tasks with Ebndry points on parent task
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(IN) ::             &
                                                 I_LO_SOUTH             &  !<-- Starting I of Sbndry region on child tasks
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
      INTEGER(kind=KINT),DIMENSION(:,:,:),POINTER,INTENT(IN) ::         &
                                                 I_INDX_PARENT_SBND     &  !<-- I indices of parent points west/east of child Sbndry point
                                                ,I_INDX_PARENT_NBND     &  !<-- I indices of parent points west/east of child Nbndry point
                                                ,I_INDX_PARENT_WBND     &  !<-- I indices of parent points west/east of child Wbndry point
                                                ,I_INDX_PARENT_EBND     &  !<-- I indices of parent points west/east of child Ebndry point
                                                ,J_INDX_PARENT_SBND     &  !<-- J indices of parent points south/north of child Sbndry point
                                                ,J_INDX_PARENT_NBND     &  !<-- J indices of parent points south/north of child Nbndry point
                                                ,J_INDX_PARENT_WBND     &  !<-- J indices of parent points south/north of child Wbndry point
                                                ,J_INDX_PARENT_EBND        !<-- J indices of parent points south/north of child Ebndry point
!
      REAL(kind=KFPT),INTENT(IN) :: PT                                  &  !<-- Top pressure of model domain (Pa)
                                   ,PDTOP                                  !<-- Pressure at top of sigma domain (Pa)
!
      REAL(kind=KFPT),DIMENSION(:),POINTER,INTENT(IN) :: PSGML1         &  !<-- Midlayer pressures, pressure domain
                                                        ,SGML2          &  !<-- Midlayer sigmas, sigma domain
                                                        ,SG1            &  !<-- Interface sigmas, pressure domain
                                                        ,SG2               !<-- Interface sigmas, sigma domain
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PD          !<-- Parent PD
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(IN) ::   &
                                                           VBL_PARENT      !<-- Parent mass point variable
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER,INTENT(IN) ::            &
                                                         WEIGHT_SBND    &  !<-- Bilinear interp weights of the 4 parent points around 
                                                        ,WEIGHT_NBND    &  !    each point on child's boundaries.
                                                        ,WEIGHT_WBND    &  !   
                                                        ,WEIGHT_EBND       !  
!
      CHARACTER(len=*),INTENT(IN) :: VBL_FLAG                              !<-- Which variable is the parent interpolating?
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(IN) :: PD_SBND        &  !<-- Boundary region PD (Pa) (column mass in sigma domain)
                                                        ,PD_NBND        &  !    on the four sides of the child boundary.
                                                        ,PD_WBND        &  !
                                                        ,PD_EBND           !
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(OUT) :: VBL_CHILD_SBND &  !<-- Mass variable in child bndry region as computed
                                                         ,VBL_CHILD_NBND &  !<-- by parent.
                                                         ,VBL_CHILD_WBND &  !
                                                         ,VBL_CHILD_EBND    !
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J,L                                       &
                           ,I_EAST,I_WEST,J_SOUTH,J_NORTH               &
                           ,I_START,I_START_EXPAND                      &
                           ,I_END,I_END_EXPAND                          &
                           ,J_START,J_START_EXPAND                      &
                           ,J_END,J_END_EXPAND                          &
                           ,KNT_PTS,KNT_PTS_X                           &
                           ,LOC_1,LOC_2                                 &
                           ,N_ADD,N_EXP,N_SIDE,N_STRIDE,NTX             &
                           ,NUM_LEVS_SPLINE,NUM_TASKS_SEND              &
                           ,RC
!
      INTEGER(kind=KINT),DIMENSION(:,:,:),POINTER :: I_INDX_PARENT_BND  &
                                                    ,J_INDX_PARENT_BND
!
      REAL(kind=KFPT) :: COEFF_1,DELP_EXTRAP,DP1,DP2,DP3                &
                        ,PDTOP_PT,PROD1,PROD2,PROD3,R_DELP
!
      REAL(kind=KFPT) :: PX_NE,PX_NW,PX_SE,PX_SW                        &
                        ,WGHT_NE,WGHT_NW,WGHT_SE,WGHT_SW
!
      REAL(kind=KFPT),DIMENSION(1:LM) :: PMID_CHILD
!
      REAL(kind=KFPT),DIMENSION(1:LM+1) :: P_INPUT                      &
                                          ,SEC_DERIV                    &
                                          ,VBL_INPUT
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: T_LOWEST                 
!
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE :: PINT_INTERP_HI      &
                                                   ,PMID_INTERP         &
                                                   ,VBL_INTERP
!
      REAL(kind=KFPT),DIMENSION(:,:,:),ALLOCATABLE :: PINT_INTERP_LO
!
      REAL(kind=KFPT),DIMENSION(1:LM+1,1:4) :: C_TMP                       !<-- Working array for ESSL spline call
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: VBL_COL_CHILD
!
      REAL(kind=KFPT),DIMENSION(:,:,:),POINTER :: WEIGHT_BND
!
      LOGICAL(kind=KLOG),DIMENSION(:),ALLOCATABLE :: INVERSION
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
!***  Loop through the four sides of the nest domain boundary (S,N,W,E).
!***  We use some dummy variables/pointers generically for all four
!***  of the sides.
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
!.......................................................................
!$omp parallel do                                                       &
!$omp private(c_tmp,coeff_1,delp_extrap,dp1,dp2,dp3,                    &
!$omp         i,i_east,i_end,i_end_expand,i_start,i_start_expand,       &
!$omp         i_west,inversion,                                         &
!$omp         j,j_end,j_end_expand,j_north,j_south,                     &
!$omp         j_start,j_start_expand,                                   &
!$omp         knt_pts,knt_pts_x,loc_1,l,loc_2,                          &
!$omp         n_add,n_stride,ntx,num_levs_spline,                       &
!$omp         p_input,pdtop_pt,pint_interp_hi,pint_interp_lo,           &
!$omp         pmid_child,pmid_interp,prod1,prod2,prod3,                 &
!$omp         px_ne,px_nw,px_se,px_sw,r_delp,                           &
!$omp         t_lowest,vbl_col_child,vbl_input,vbl_interp,              &
!$omp         wght_ne,wght_nw,wght_se,wght_sw)
!.......................................................................
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
          ALLOCATE(T_LOWEST   (1:N_STRIDE))
          ALLOCATE(INVERSION  (1:N_STRIDE))
!
          ALLOCATE(PINT_INTERP_HI(I_START:I_END,J_START:J_END))
          ALLOCATE(PINT_INTERP_LO(I_START:I_END,J_START:J_END,1:NLEV+1))
!
!-----------------------------------------------------------------------
!***  We need the mid-layer pressure values in the parent layers
!***  over the child boundary point locations since those are
!***  required for the vertical interpolation of variables
!***  to the mid-layers in the child.
!***  Compute the interface pressures of the parent layers
!***  then take the means to get the mid-layer values.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Start with the bottom layer (L=NLEV).
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
            PX_SW=PD(I_WEST,J_SOUTH)+PT                                    !<-- Sfc pressure on parent point SW of nest point
            PX_SE=PD(I_EAST,J_SOUTH)+PT                                    !<-- Sfc pressure on parent point SE of nest point
            PX_NW=PD(I_WEST,J_NORTH)+PT                                    !<-- Sfc pressure on parent point NW of nest point
            PX_NE=PD(I_EAST,J_NORTH)+PT                                    !<-- Sfc pressure on parent point NE of nest point
!
            PINT_INTERP_LO(I,J,NLEV+1)=WGHT_SW*PX_SW                    &  !<-- Parent's surface pressure interp'd to this child's
                                      +WGHT_SE*PX_SE                    &  !    gridpoint (I,J) along child's boundary for
                                      +WGHT_NW*PX_NW                    &  !    child task NTX.
                                      +WGHT_NE*PX_NE
!
          ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!***  Now compute those mid-layer pressures in the parent layers
!***  as well as the values of the parent's variable at those
!***  pressure levels.
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
              PX_SW=SG2(L)*PD(I_WEST,J_SOUTH)+PDTOP_PT                     !<-- Pressure, top of layer L, parent point SW of nest point
              PX_SE=SG2(L)*PD(I_EAST,J_SOUTH)+PDTOP_PT                     !<-- Pressure, top of layer L, parent point SE of nest point
              PX_NW=SG2(L)*PD(I_WEST,J_NORTH)+PDTOP_PT                     !<-- Pressure, top of layer L, parent point NW of nest point
              PX_NE=SG2(L)*PD(I_EAST,J_NORTH)+PDTOP_PT                     !<-- Pressure, top of layer L, parent point NE of nest point
!
              WGHT_SW=WEIGHT_BND(I,J,INDX_SW)                              !<-- Bilinear weight for parent's point SW of child's point
              WGHT_SE=WEIGHT_BND(I,J,INDX_SE)                              !<-- Bilinear weight for parent's point SE of child's point
              WGHT_NW=WEIGHT_BND(I,J,INDX_NW)                              !<-- Bilinear weight for parent's point NW of child's point
              WGHT_NE=WEIGHT_BND(I,J,INDX_NE)                              !<-- Bilinear weight for parent's point NE of child's point
!
              PINT_INTERP_HI(I,J)=WGHT_SW*PX_SW                         &  !<-- Top interface pressure interp'd to child gridpoint
                                 +WGHT_SE*PX_SE                         &  !    in child's boundary region on child task NTX
                                 +WGHT_NW*PX_NW                         &
                                 +WGHT_NE*PX_NE
!
              KNT_PTS=KNT_PTS+1
!
              PMID_INTERP(KNT_PTS,L)=0.5*(PINT_INTERP_HI(I,J)           &  !<-- Parent midlayer pressure interp'd to child gridpoint
                                         +PINT_INTERP_LO(I,J,L+1))         !    in child's boundary region of child task NTX
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
              PINT_INTERP_LO(I,J,L)=PINT_INTERP_HI(I,J)
!
            ENDDO
            ENDDO
!
          ENDDO
!
!-----------------------------------------------------------------------
!***  The parent uses mass-weighted averages of the lowest and 
!***  2nd lowest layer temperatures to avoid extrapolation problems
!***  when there is an inversion present.
!-----------------------------------------------------------------------
!
          IF(TRIM(VBL_FLAG)=='T')THEN
            KNT_PTS=0
!
            DO J=J_START,J_END   
            DO I=I_START,I_END  
              KNT_PTS=KNT_PTS+1
              INVERSION(KNT_PTS)=.FALSE.
              IF(VBL_INTERP(KNT_PTS,NLEV)<VBL_INTERP(KNT_PTS,NLEV-1))THEN
                DP1=PINT_INTERP_LO(I,J,NLEV+1)-PINT_INTERP_LO(I,J,NLEV)
                DP2=PINT_INTERP_LO(I,J,NLEV)-PINT_INTERP_LO(I,J,NLEV-1)
                DP3=PINT_INTERP_LO(I,J,NLEV-1)-PINT_INTERP_LO(I,J,NLEV-2)
 
                PROD1=VBL_INTERP(KNT_PTS,NLEV)*DP1
                PROD2=VBL_INTERP(KNT_PTS,NLEV-1)*DP2
                PROD3=VBL_INTERP(KNT_PTS,NLEV-2)*DP3
!
                INVERSION(KNT_PTS)=.TRUE.
                T_LOWEST(KNT_PTS)=VBL_INTERP(KNT_PTS,NLEV)
                VBL_INTERP(KNT_PTS,NLEV)=(PROD1+PROD2+PROD3)/(DP1+DP2+DP3)
                VBL_INTERP(KNT_PTS,NLEV-1)=(PROD2+PROD3)/(DP2+DP3)
              ENDIF
            ENDDO
            ENDDO
          ENDIF
!
!-----------------------------------------------------------------------
!***  Compute values of the variable at mid-layers of the child
!***  over nest boundary region points based on the parent's
!***  interpolated values.
!-----------------------------------------------------------------------
!
          KNT_PTS  =0
          KNT_PTS_X=0
          LOC_1=0
          N_ADD=(LM-1)*N_STRIDE
!
!-----------------------------------------------------------------------
!***  Recall that the PDB_H pointer contains data for one extra row
!***  beyond the child boundary on all sides since we need the
!***  ability to do 4-PT averages of PD to V points.   The PDB_V
!***  pointer did not need the extra row and so does not have it.
!***  Here we must use PDB but the values for the prognostic variables
!***  are only generated on the true boundary points.  To address
!***  PDB correctly we loop over the potentially expanded boundary
!***  region but CYCLE when we are at points not within the true
!***  boundary region of the child.
!-----------------------------------------------------------------------
!
          main_loop: DO J=J_START_EXPAND,J_END_EXPAND                      !<-- J limits of expanded child task bndry region on parent task
!
          DO I=I_START_EXPAND,I_END_EXPAND                                 !<-- Child's expanded boundary I values on this parent task
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
!***  Midlayer pressures for the nest boundary points.
!-----------------------------------------------------------------------
!
            DO L=1,LM                                 
              PMID_CHILD(L)=PSGML1(L)+SGML2(L)*PDB(NTX)%DATA(KNT_PTS_X)    !<-- Nest midlayer pressure
            ENDDO
! 
!-----------------------------------------------------------------------
!***  Use spline interpolation to move variables from their
!***  vertical locations in the column after horizontal interpolation
!***  from the surrounding parent points to child boundary point levels.
!***  The target locations are the new midlayer pressures in the
!***  nest boundary point columns based on the new surface pressure
!***  for the nest's terrain.
!
!***  If the target location lies below the middle of the lowest parent
!***  layer in the newly created child column then extrapolate linearly
!***  in pressure to obtain a value at the lowest child mid-layer and
!***  fill in the remaining 'underground' levels using the call to
!***  'SPLINE' just as with all the other higher levels.
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
              NUM_LEVS_SPLINE=LM+1                                         !<-- Add another input level at nest's lowest
              P_INPUT(LM+1)=PMID_CHILD(LM)                                 !    mid-layer.
!
!             IF(TRIM(VBL_FLAG)=='T'.AND.INVERSION(KNT_PTS))THEN           !<-- For temperature inversions place the parent's
!               VBL_INPUT(LM+1)=T_LOWEST(KNT_PTS)                          !    original cold sfc lyr into the new bottom lyr.
!             ELSE
                R_DELP=1./(P_INPUT(LM)-P_INPUT(LM-1))
                DELP_EXTRAP=PMID_CHILD(LM)-P_INPUT(LM)
!
                COEFF_1=(VBL_INPUT(LM)-VBL_INPUT(LM-1))*R_DELP
                VBL_INPUT(LM+1)=VBL_INPUT(LM)+COEFF_1*DELP_EXTRAP          !<-- Extrapolated value at nest's new bottom layer.
!             ENDIF
!
            ENDIF
!
#ifdef IBM
            CALL SCSINT(P_INPUT                                         &  !<-- Input mid-layer pressure
                       ,VBL_INPUT                                       &  !<-- Input mid-layer mass variable value
                       ,C_TMP                                           &  !<-- Auxiliary matrix, C(1:num_levs_spline,1:4)
                       ,NUM_LEVS_SPLINE                                 &  !<-- # of input levels
                       ,0                                               &  !<-- Natural boundary conditions, nothing precomputed (zero or negative integer)
                       ,PMID_CHILD                                      &  !<-- Child mid-layer pressures to interpolate to
                       ,VBL_COL_CHILD                                   &  !<-- Child mid-layer variable value returned
                       ,LM)                                                !<-- # of child mid-layers to interpolate to

#else
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
#endif
!
!-----------------------------------------------------------------------
!
          ENDDO
          ENDDO main_loop
!
!-----------------------------------------------------------------------
!
          DEALLOCATE(PMID_INTERP)
          DEALLOCATE(VBL_INTERP)
          DEALLOCATE(INVERSION)
          DEALLOCATE(T_LOWEST)
!
          DEALLOCATE(PINT_INTERP_HI)
          DEALLOCATE(PINT_INTERP_LO)
!
!-----------------------------------------------------------------------
!
        ENDDO child_task_loop
!
!.......................................................................
!$omp end parallel do
!.......................................................................
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
                                         ,N_BLEND_H_CHILD               & !     |
                                         ,N_BLEND_V_CHILD               & !     |
                                         ,IM_CHILD                      & !     |
                                         ,JM_CHILD                      & !     |
                                                                          !     |
                                         ,INC_FIX_N                     & !   Input
!                                                                           ---------
                                         ,PD_V                          & !   Output
                                         ,PDB_S_V                       & !     |
                                         ,PDB_N_V                       & !     |
                                         ,PDB_W_V                       & !     v
                                         ,PDB_E_V )                       !  
!
!-----------------------------------------------------------------------
!***  Use 4-pt horizontal interpolation to compute PD on V points
!***  of the parent domain and of the nest boundary given those 
!***  values on H points.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: IMS,IME,JMS,JME                  &  !<-- Parent task's memory limits
                                      ,IM_CHILD,JM_CHILD                &  !<-- Nest domain limits
                                      ,INC_FIX_N                        &  !<-- Increment for selecting nest tasks for averaging H to V
                                      ,N_BLEND_H_CHILD                  &  !<-- H rows in nest's boundary region
                                      ,N_BLEND_V_CHILD                  &  !<-- V rows in nest's boundary region
                                      ,NUM_TASKS_SEND_SBND              &  !<-- # of child tasks with Sbndry V points on parent task
                                      ,NUM_TASKS_SEND_NBND              &  !<-- # of child tasks with Nbndry V points on parent task
                                      ,NUM_TASKS_SEND_WBND              &  !<-- # of child tasks with Wbndry V points on parent task
                                      ,NUM_TASKS_SEND_EBND                 !<-- # of child tasks with Ebndry V points on parent task
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(IN) ::             &
                                                        I_LO_SOUTH_H    &  !<-- Starting I of Sbndry region H points on child tasks
                                                       ,I_HI_SOUTH_H    &  !<-- Ending I of Sbndry region H points on child tasks
                                                       ,I_LO_NORTH_H    &  !<-- Starting I of Nbndry region H points on child tasks
                                                       ,I_HI_NORTH_H    &  !<-- Ending I of Nbndry region H points on child tasks
                                                       ,J_LO_WEST_H     &  !<-- Starting J of Wbndry region H points on child tasks
                                                       ,J_HI_WEST_H     &  !<-- Ending J of Wbndry region H points on child tasks
                                                       ,J_LO_EAST_H     &  !<-- Starting J of Ebndry region H points on child tasks
                                                       ,J_HI_EAST_H        !<-- Ending J of Ebndry region H points on child tasks
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(IN) ::             &
                                                        I_LO_SOUTH_V    &  !<-- Starting I of Sbndry region V points on child tasks
                                                       ,I_HI_SOUTH_V    &  !<-- Ending I of Sbndry region V points on child tasks
                                                       ,I_LO_NORTH_V    &  !<-- Starting I of Nbndry region V points on child tasks
                                                       ,I_HI_NORTH_V    &  !<-- Ending I of Nbndry region V points on child tasks
                                                       ,J_LO_WEST_V     &  !<-- Starting J of Wbndry region V points on child tasks
                                                       ,J_HI_WEST_V     &  !<-- Ending J of Wbndry region V points on child tasks
                                                       ,J_LO_EAST_V     &  !<-- Starting J of Ebndry region V points on child tasks
                                                       ,J_HI_EAST_V        !<-- Ending J of Ebndry region V points on child tasks
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PD_H        !<-- Parent PD (Pa) (column mass in sigma domain) on H points
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(IN) :: PDB_S_H        &  !<-- Boundary region PD (Pa) (column mass in sigma domain)
                                                        ,PDB_N_H        &  !    on mass points on the four sides of the child boundary.
                                                        ,PDB_W_H        &  !
                                                        ,PDB_E_H           !<--
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: PD_V       !<-- Parent PD (Pa) (column mass in sigma domain) on V points
!
      TYPE(REAL_DATA),DIMENSION(:),POINTER,INTENT(OUT) :: PDB_S_V       &  !<-- Child boundary PD (Pa) on child domain Sbndry V points
                                                         ,PDB_N_V       &  !<-- Child boundary PD (Pa) on child domain Nbndry V points
                                                         ,PDB_W_V       &  !<-- Child boundary PD (Pa) on child domain Wbndry V points
                                                         ,PDB_E_V          !<-- Child boundary PD (Pa) on child domain Ebndry V points
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: DIFF_START,DIFF_START_PTS                   &
                           ,I,J                                         &
                           ,I_ADD_H,I_INC_H                             &
                           ,I_START_H,I_START_V                         &
                           ,I_END_H,I_END_V                             &
                           ,I_HI_H,I_LO_H                               &
                           ,J_START_H,J_START_V                         &
                           ,J_END_H,J_END_V                             &
                           ,KNT_H,KNT_V                                 &
                           ,N_ADD_H,N_OFFSET_H,N_SIDE                   &
                           ,NTX,NTX_H,NUM_TASKS_SEND
!
      INTEGER(kind=KINT),DIMENSION(4) :: NTX_ADD
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
!***  First obtain PD on the parent's V points.
!-----------------------------------------------------------------------
!
      DO J=JMS,JME-1
      DO I=IMS,IME-1
        PD_V(I,J)=0.25*(PD_H(I,J)+PD_H(I+1,J)+PD_H(I,J+1)+PD_H(I+1,J+1))
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Now average values of PD on and adjacent to the nests' boundary
!***  to obtain PD on the nests' boundary V points.
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
            J_END_V  =N_BLEND_V_CHILD                                      !<-- J index of last Sbndry V point on this child task
!
            I_END_H  =I_HI_SOUTH_H(NTX)                                    !<-- I index of last Sbndry H point east of last Sbndry V point
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
            J_START_V=JM_CHILD-N_BLEND_V_CHILD                             !<-- J index of first Nbndry V point on this child task
            J_END_V  =JM_CHILD-1                                           !<-- J index of last Nbndry V point on this child task
!
            I_END_H  =I_HI_NORTH_H(NTX)                                    !<-- I index of last Nbndry H point east of last Nbndry V point
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
            I_END_V  =N_BLEND_V_CHILD                                      !<-- I index of last Wbndry V point on this child task
            J_START_V=J_LO_WEST_V(NTX)                                     !<-- J index of first Wbndry V point on this child task
            J_END_V  =J_HI_WEST_V(NTX)                                     !<-- J index of last Wbndry V point on this child task
!
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
            DIFF_START_PTS=DIFF_START*(N_BLEND_H_CHILD+1)
            N_ADD_H=N_BLEND_H_CHILD+1
!
!---------------------------------------------
!***  East boundary limits on this child task
!---------------------------------------------
!
          ELSEIF(N_SIDE==4)THEN
            I_START_V=IM_CHILD-N_BLEND_V_CHILD                             !<-- I index of first Ebndry V point on this child task
            I_END_V  =IM_CHILD-1                                           !<-- I index of last Ebndry V point on this child task
            J_START_V=J_LO_EAST_V(NTX)                                     !<-- J index of first Ebndry V point on this child task
            J_END_V  =J_HI_EAST_V(NTX)                                     !<-- J index of last Ebndry V point on this child task
!
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
            DIFF_START_PTS=DIFF_START*(N_BLEND_H_CHILD+1)
            N_ADD_H=N_BLEND_H_CHILD+1
!
          ENDIF
!
!-----------------------------------------------------------------------
!***  Ready to average nest boundary PD on H to PD on V.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Recall that the PDB_H pointer contains data for one extra row
!***  beyond the child boundary on all sides since we need the
!***  ability to do 4-pt averages of PD to V points.
!***  Here we must use PDB_H but the values for PDB_V (and thus U,V)
!***  are only generated on the true boundary points.  To address PDB_H
!***  correctly we must take into account those extra points in PDB_H
!***  as we march through the V point locations.
!
!***  Also we must be aware that the nest boundary V segments that are
!***  being considered for each nest task must correspond to the same
!***  nest task H point values of PDB.  Because of overlap that can
!***  occur due to the segments of PDB on H being larger than PBD on
!***  V we must explicitly check to be sure that nest PDB on V is
!***  being filled by nest PDB on H for the same nest task.
!-----------------------------------------------------------------------
!
          KNT_V=0
!
          DO J=J_START_V,J_END_V
            IF(N_SIDE<=2)THEN
              N_OFFSET_H=(DIFF_START+1)*(J-J_START_V+1)-1
            ELSEIF(N_SIDE>=3)THEN
              N_OFFSET_H=J-J_START_V+DIFF_START_PTS
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
!
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
!***  Extract variables for a nest's boundary from the datastring
!***  received by the child from its parent.  A child task might 
!***  receive segments of boundary data from two parent tasks in 
!***  which case the pieces are be combined.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: LENGTH_DATA                      &  !<-- # of words in datastring
                                      ,I_START                          &  !<-- Starting I of data segment in string on child's grid
                                      ,I_END                            &  !<-- Ending I of data segment in string on child's grid
                                      ,J_START                          &  !<-- Starting J of data segment in string on child's grid
                                      ,J_END                            &  !<-- Ending J of data segment in string on child's grid
                                      ,ILIM_LO                          &  !<-- Lower I limit of full boundary segment 
                                      ,ILIM_HI                          &  !<-- Upper I limit of full boundary segment 
                                      ,JLIM_LO                          &  !<-- Lower J limit of full boundary segment 
                                      ,JLIM_HI                             !<-- Upper J limit of full boundary segment 
!
      REAL(kind=KFPT),DIMENSION(:),INTENT(IN) :: DATASTRING                !<-- The string of boundary data from the parent
!      
      REAL(kind=KFPT),DIMENSION(ILIM_LO:ILIM_HI,JLIM_LO:JLIM_HI)        &
                                          ,INTENT(OUT),OPTIONAL :: PDB     !<-- PD for segment of the child boundary
!      
      REAL(kind=KFPT),DIMENSION(ILIM_LO:ILIM_HI,JLIM_LO:JLIM_HI,1:LM)   &
                                           ,INTENT(OUT),OPTIONAL :: TB  &  !<-- Temperature for segment of the child boundary
                                                                   ,QB  &  !<-- Specific humidity for segment of the child boundary
                                                                  ,CWB  &  !<-- Cloud condensate for segment of the child boundary
                                                                   ,UB  &  !<-- U wind for segment of the child boundary
                                                                   ,VB     !<-- V wind for segment of the child boundary
!      
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J,K,N
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the appropriate boundary variables from the datastring
!***  provided by the parent.
!***  Here we transfer from the 1-D string to the 2-D and 3-D arrays.
!***  This easily allows for proper combining of separate strings
!***  whose ends may overlap that are arriving from different parent
!***  tasks.  To export these boundary arrays out of the coupler
!***  though the data must be put back into 1-D.
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
!***  Load the child boundary values received from the parent 
!***  into the Parent-Child coupler export state.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: ILIM_LO                          &  !<-- Lower I limit of full boundary segment
                                      ,ILIM_HI                          &  !<-- Upper I limit of full boundary segment
                                      ,JLIM_LO                          &  !<-- Lower J limit of full boundary segment
                                      ,JLIM_HI                             !<-- Upper J limit of full boundary segment
!      
      REAL(kind=KFPT),DIMENSION(ILIM_LO:ILIM_HI,JLIM_LO:JLIM_HI)        &
                                           ,INTENT(IN),OPTIONAL :: PDB     !<-- PD for segment of the child boundary
!      
      REAL(kind=KFPT),DIMENSION(ILIM_LO:ILIM_HI,JLIM_LO:JLIM_HI,1:LM)   &
                                            ,INTENT(IN),OPTIONAL :: TB  &  !<-- Temperature for segment of the child boundary
                                                                   ,QB  &  !<-- Specific humidity for segment of the child boundary
                                                                  ,CWB  &  !<-- Cloud condensate for segment of the child boundary
                                                                   ,UB  &  !<-- U wind for segment of the child boundary
                                                                   ,VB     !<-- V wind for segment of the child boundary
!
      CHARACTER(*),INTENT(IN) :: DATA_NAME                                 !<-- Name used for each child task's boundary segment
!
      REAL(kind=KFPT),DIMENSION(:),POINTER,INTENT(INOUT) :: DATA_EXP       !<-- Combined boundary segment data on child task
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXPORT_STATE                       !<-- The Parent-Child Coupler export state
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,J,K                                       &
                           ,ISTAT,RC,RC_EXP_BNDRY
!
      INTEGER(kind=KINT),SAVE :: NN
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC          =ESMF_SUCCESS
      RC_EXP_BNDRY=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  The children's final boundary data will be loaded into the
!***  Parent-Child Coupler's export state as 1-D arrays (Attributes)
!***  since they are not spanning the childrens' ESMF Grids (as Fields).
!***  But because they are Attributes and not Fields, they will need
!***  to be reset every time.
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
#ifdef ESMF_3
        CALL ESMF_AttributeSet(state    =EXPORT_STATE                   &  !<-- This Parent_child Coupler export state
                              ,name     =DATA_NAME                      &  !<-- Name of the children's new boundary H data
                              ,count    =NN                             &  !<-- # of words in the data
                              ,valueList=DATA_EXP                       &  !<-- The children's new boundary H data
                              ,rc       =RC)
#else
        CALL ESMF_AttributeSet(state    =EXPORT_STATE                   &  !<-- This Parent_child Coupler export state
                              ,name     =DATA_NAME                      &  !<-- Name of the children's new boundary H data
                              ,itemCount=NN                             &  !<-- # of words in the data
                              ,valueList=DATA_EXP                       &  !<-- The children's new boundary H data
                              ,rc       =RC)
#endif
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
#ifdef ESMF_3
        CALL ESMF_AttributeSet(state    =EXPORT_STATE                    &  !<-- This Parent_child Coupler export state
                              ,name     =DATA_NAME                       &  !<-- Name of the children's new boundary V data
                              ,count    =NN                              &  !<-- # of words in the data
                              ,valueList=DATA_EXP                        &  !<-- The children's new boundary V data
                              ,rc       =RC)
#else
        CALL ESMF_AttributeSet(state    =EXPORT_STATE                    &  !<-- This Parent_child Coupler export state
                              ,name     =DATA_NAME                       &  !<-- Name of the children's new boundary V data
                              ,itemCount=NN                              &  !<-- # of words in the data
                              ,valueList=DATA_EXP                        &  !<-- The children's new boundary V data
                              ,rc       =RC)
#endif
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
!
      SUBROUTINE DEALLOC_WORK_PARENTS(N,TIME_FLAG)
!
!-----------------------------------------------------------------------
!***  Parents deallocate all working pointers that needed to be
!***  allocated with unique dimensions at the outset of the forecast 
!***  and again each time a nest moves.  Those allocations took place
!***  in subroutine POINT_INTERP_DATA_TO_MEMORY.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: N                                   !<-- Which child of this parent
!
      CHARACTER(*),INTENT(IN) :: TIME_FLAG                                 !<-- Current or future boundary data for children?
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: IERR,INDX2,ISTAT,NCHILD_TASKS,NMAX,NT
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Select the appropriate part of the working array depending on
!***  whether we are now concerned with children's boundaries for
!***  their current time or for their future.
!-----------------------------------------------------------------------
!
      IF(TIME_FLAG=='Future')THEN
        INDX2=1
      ELSEIF(TIME_FLAG=='Current')THEN
        INDX2=2
      ENDIF
!
!-----------------------------------------------------------------------
!***  The parent must not deallocate the actual BC update buffers
!***  until it knows that they have been Recvd by the children.
!***  That is the reason for the Waits below.
!-----------------------------------------------------------------------
!
!-----------
!***  South
!-----------
!
      NMAX=UBOUND(HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV,1)
      IF(NMAX>0)THEN
        DO NT=1,NMAX
          CALL MPI_WAIT(HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
      south_h: IF(ASSOCIATED(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS))THEN
!
        NCHILD_TASKS=UBOUND(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS,1)          !<-- # of Sbndry tasks on child N that recvd H point data
!
        DO NT=1,NCHILD_TASKS
          IF(ASSOCIATED(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA))THEN
            DEALLOCATE(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS(NT)%DATA)
          ENDIF
        ENDDO
!
        DEALLOCATE(CHILD_BOUND_H_SOUTH(N,INDX2)%TASKS,stat=ISTAT)
!
        IF(ASSOCIATED(HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV))THEN
          DEALLOCATE(HANDLE_H_SOUTH(N,INDX2)%NTASKS_TO_RECV)
        ENDIF
!
      ENDIF south_h
!
      NCHILD_TASKS=UBOUND(PD_B_SOUTH(N)%TASKS,1)                           !<-- # of Sbndry tasks on child N that recvd H point data
      DO NT=1,NCHILD_TASKS
        IF(ASSOCIATED(PD_B_SOUTH(N)%TASKS(NT)%DATA))THEN
          DEALLOCATE(PD_B_SOUTH(N)%TASKS(NT)%DATA)
        ENDIF
        DEALLOCATE(FIS_CHILD_SOUTH(N)%TASKS(NT)%DATA)
      ENDDO
!
      DEALLOCATE(PD_B_SOUTH(N)%TASKS)
      DEALLOCATE( T_B_SOUTH(N)%TASKS)
      DEALLOCATE( Q_B_SOUTH(N)%TASKS)
      DEALLOCATE(CW_B_SOUTH(N)%TASKS)
      DEALLOCATE(FIS_CHILD_SOUTH(N)%TASKS)
      DEALLOCATE(WORDS_BOUND_H_SOUTH(N)%TASKS)
!
      NMAX=UBOUND(HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV,1)
      IF(NMAX>0)THEN
        DO NT=1,NMAX
          CALL MPI_WAIT(HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
      south_v: IF(ASSOCIATED(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS))THEN
!
        NCHILD_TASKS=UBOUND(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS,1)          !<-- # of Sbndry tasks on child N that recvd V point data
!
        DO NT=1,NCHILD_TASKS
          IF(ASSOCIATED(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(NT)%DATA))THEN
            DEALLOCATE(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS(NT)%DATA)
          ENDIF
        ENDDO
!
        DEALLOCATE(CHILD_BOUND_V_SOUTH(N,INDX2)%TASKS)
!
        IF(ASSOCIATED(HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV))THEN
          DEALLOCATE(HANDLE_V_SOUTH(N,INDX2)%NTASKS_TO_RECV)
        ENDIF
!
      ENDIF south_v
!
      NCHILD_TASKS=UBOUND(PD_B_SOUTH_V(N)%TASKS,1)                         !<-- # of Sbndry tasks on child N that recvd V point data
      DO NT=1,NCHILD_TASKS
        IF(ASSOCIATED(PD_B_SOUTH_V(N)%TASKS(NT)%DATA))THEN
          DEALLOCATE(PD_B_SOUTH_V(N)%TASKS(NT)%DATA)
        ENDIF
      ENDDO
!
      DEALLOCATE(WORDS_BOUND_V_SOUTH(N)%TASKS)
      DEALLOCATE(PD_B_SOUTH_V(N)%TASKS)
      DEALLOCATE(U_B_SOUTH(N)%TASKS)
      DEALLOCATE(V_B_SOUTH(N)%TASKS)
!
!-----------------------------------------------------------------------
!
!-----------
!***  North
!-----------
!
      NMAX=UBOUND(HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV,1)
      IF(NMAX>0)THEN
        DO NT=1,NMAX
          CALL MPI_WAIT(HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
      north_h: IF(ASSOCIATED(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS))THEN
!
        NCHILD_TASKS=UBOUND(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS,1)          !<-- # of Nbndry tasks on child N that recvd H point data
!
        DO NT=1,NCHILD_TASKS
          IF(ASSOCIATED(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA))THEN
            DEALLOCATE(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,*)' DEALLOC_WORK_PARENTS N=',N,' NT=',NT,' INDX2=',INDX2 &
                       ,' Failed to deallocate CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA'
            ELSE
!             WRITE(0,*)' DEALLOC_WORK_PARENTS N=',N,' NT=',NT,' INDX2=',INDX2 &
!                      ,' Succeeded in deallocating CHILD_BOUND_H_NORTH(N,INDX2)%TASKS(NT)%DATA'
            ENDIF
          ENDIF
        ENDDO
!
        DEALLOCATE(CHILD_BOUND_H_NORTH(N,INDX2)%TASKS)
!
        IF(ASSOCIATED(HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV))THEN
          DEALLOCATE(HANDLE_H_NORTH(N,INDX2)%NTASKS_TO_RECV)
        ENDIF
!
      ENDIF north_h
!
      NCHILD_TASKS=UBOUND(PD_B_NORTH(N)%TASKS,1)                           !<-- # of Nbndry tasks on child N that recvd H point data
!
      DO NT=1,NCHILD_TASKS
        IF(ASSOCIATED(PD_B_NORTH(N)%TASKS(NT)%DATA))THEN
          DEALLOCATE(PD_B_NORTH(N)%TASKS(NT)%DATA,stat=ISTAT)
          IF(ISTAT/=0)THEN
            WRITE(0,*)' DEALLOC_WORK_PARENTS N=',N,' NT=',NT,' INDX2=',INDX2 &
                     ,' Failed to deallocate PD_B_NORTH(N)%TASKS(NT)%DATA'
          ELSE
!           WRITE(0,*)' DEALLOC_WORK_PARENTS N=',N,' NT=',NT,' INDX2=',INDX2 &
!                    ,' Succeeded in deallocating PD_B_NORTH(N)%TASKS(NT)%DATA'
          ENDIF
        ENDIF
        DEALLOCATE(FIS_CHILD_NORTH(N)%TASKS(NT)%DATA)
      ENDDO
!
      DEALLOCATE(PD_B_NORTH(N)%TASKS)
      DEALLOCATE( T_B_NORTH(N)%TASKS)
      DEALLOCATE( Q_B_NORTH(N)%TASKS)
      DEALLOCATE(CW_B_NORTH(N)%TASKS)
      DEALLOCATE(FIS_CHILD_NORTH(N)%TASKS)
      DEALLOCATE(WORDS_BOUND_H_NORTH(N)%TASKS)
!
      NMAX=UBOUND(HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV,1)
      IF(NMAX>0)THEN
        DO NT=1,NMAX
          CALL MPI_WAIT(HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV(NT)      &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
      north_v: IF(ASSOCIATED(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS))THEN
!
        NCHILD_TASKS=UBOUND(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS,1)          !<-- # of Nbndry tasks on child N that recvd V point data
        DO NT=1,NCHILD_TASKS
          IF(ASSOCIATED(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA))THEN
            DEALLOCATE(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA         &
                                                           ,stat=ISTAT)
            IF(ISTAT/=0)THEN
              WRITE(0,*)' DEALLOC_WORK_PARENTS N=',N,' NT=',NT,' INDX2=',INDX2 &
                       ,' Failed to deallocate CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA'
            ELSE
!             WRITE(0,*)' DEALLOC_WORK_PARENTS N=',N,' NT=',NT,' INDX2=',INDX2 &
!                      ,' Succeeded in deallocating CHILD_BOUND_V_NORTH(N,INDX2)%TASKS(NT)%DATA'
            ENDIF
          ENDIF
        ENDDO
!
        DEALLOCATE(CHILD_BOUND_V_NORTH(N,INDX2)%TASKS)
!
        IF(ASSOCIATED(HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV))THEN
          DEALLOCATE(HANDLE_V_NORTH(N,INDX2)%NTASKS_TO_RECV)
        ENDIF
!
      ENDIF north_v
!
      NCHILD_TASKS=UBOUND(PD_B_NORTH_V(N)%TASKS,1)                         !<-- # of Nbndry tasks on child N that recvd V point data
      DO NT=1,NCHILD_TASKS
        IF(ASSOCIATED(PD_B_NORTH_V(N)%TASKS(NT)%DATA))THEN
          DEALLOCATE(PD_B_NORTH_V(N)%TASKS(NT)%DATA,stat=ISTAT)
          IF(ISTAT/=0)THEN
            WRITE(0,*)' DEALLOC_WORK_PARENTS N=',N,' NT=',NT,' INDX2=',INDX2 &
                     ,' Failed to deallocate PD_B_NORTH_V(N)%TASKS(NT)%DATA'
          ELSE
!           WRITE(0,*)' DEALLOC_WORK_PARENTS N=',N,' NT=',NT,' INDX2=',INDX2 &
!                    ,' Succeeded in deallocating PD_B_NORTH_V(N)%TASKS(NT)%DATA'
          ENDIF
        ENDIF
      ENDDO
!
      DEALLOCATE(WORDS_BOUND_V_NORTH(N)%TASKS)
      DEALLOCATE(PD_B_NORTH_V(N)%TASKS)
      DEALLOCATE(U_B_NORTH(N)%TASKS)
      DEALLOCATE(V_B_NORTH(N)%TASKS)
!
!-----------------------------------------------------------------------
!
!----------
!***  West
!----------
!
      NMAX=UBOUND(HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV,1)
      IF(NMAX>0)THEN
        DO NT=1,NMAX
          CALL MPI_WAIT(HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV(NT)       &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
      west_h: IF(ASSOCIATED(CHILD_BOUND_H_WEST(N,INDX2)%TASKS))THEN
!
        NCHILD_TASKS=UBOUND(CHILD_BOUND_H_WEST(N,INDX2)%TASKS,1)           !<-- # of Wbndry tasks on child N that recvd H point data
        DO NT=1,NCHILD_TASKS
          IF(ASSOCIATED(CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA))THEN
            DEALLOCATE(CHILD_BOUND_H_WEST(N,INDX2)%TASKS(NT)%DATA)
          ENDIF
        ENDDO
!
        DEALLOCATE(CHILD_BOUND_H_WEST(N,INDX2)%TASKS)
!
        IF(ASSOCIATED(HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV))THEN
          DEALLOCATE(HANDLE_H_WEST(N,INDX2)%NTASKS_TO_RECV)
        ENDIF
!
      ENDIF west_h
!
      NCHILD_TASKS=UBOUND(PD_B_WEST(N)%TASKS,1)                            !<-- # of Wbndry tasks on child N that recvd H point data
      DO NT=1,NCHILD_TASKS
        IF(ASSOCIATED(PD_B_WEST(N)%TASKS(NT)%DATA))THEN
          DEALLOCATE(PD_B_WEST(N)%TASKS(NT)%DATA)
        ENDIF
        DEALLOCATE(FIS_CHILD_WEST(N)%TASKS(NT)%DATA)
      ENDDO
!
      DEALLOCATE(PD_B_WEST(N)%TASKS)
      DEALLOCATE( T_B_WEST(N)%TASKS)
      DEALLOCATE( Q_B_WEST(N)%TASKS)
      DEALLOCATE(CW_B_WEST(N)%TASKS)
      DEALLOCATE(FIS_CHILD_WEST(N)%TASKS)
      DEALLOCATE(WORDS_BOUND_H_WEST(N)%TASKS)
!
      NMAX=UBOUND(HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV,1)
      IF(NMAX>0)THEN
        DO NT=1,NMAX
          CALL MPI_WAIT(HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV(NT)       &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
      west_v: IF(ASSOCIATED(CHILD_BOUND_V_WEST(N,INDX2)%TASKS))THEN
!
        NCHILD_TASKS=UBOUND(CHILD_BOUND_V_WEST(N,INDX2)%TASKS,1)           !<-- # of Wbndry tasks on child N that recvd V point data
!
        DO NT=1,NCHILD_TASKS
          IF(ASSOCIATED(CHILD_BOUND_V_WEST(N,INDX2)%TASKS(NT)%DATA))THEN
            DEALLOCATE(CHILD_BOUND_V_WEST(N,INDX2)%TASKS(NT)%DATA)
          ENDIF
        ENDDO
!
        DEALLOCATE(CHILD_BOUND_V_WEST(N,INDX2)%TASKS)
!
        IF(ASSOCIATED(HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV))THEN
          DEALLOCATE(HANDLE_V_WEST(N,INDX2)%NTASKS_TO_RECV)
        ENDIF
!
      ENDIF west_v
!
      NCHILD_TASKS=UBOUND(PD_B_WEST_V(N)%TASKS,1)                          !<-- # of Wbndry tasks on child N that recvd V point data
      DO NT=1,NCHILD_TASKS
        IF(ASSOCIATED(PD_B_WEST_V(N)%TASKS(NT)%DATA))THEN
          DEALLOCATE(PD_B_WEST_V(N)%TASKS(NT)%DATA)
        ENDIF
      ENDDO
!
      DEALLOCATE(WORDS_BOUND_V_WEST(N)%TASKS)
      DEALLOCATE(PD_B_WEST_V(N)%TASKS)
      DEALLOCATE(U_B_WEST(N)%TASKS)
      DEALLOCATE(V_B_WEST(N)%TASKS)
!
!-----------------------------------------------------------------------
!
!----------
!***  East
!----------
!
      NMAX=UBOUND(HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV,1)
      IF(NMAX>0)THEN
        DO NT=1,NMAX
          CALL MPI_WAIT(HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV(NT)       &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
      east_h: IF(ASSOCIATED(CHILD_BOUND_H_EAST(N,INDX2)%TASKS))THEN
!
        NCHILD_TASKS=UBOUND(CHILD_BOUND_H_EAST(N,INDX2)%TASKS,1)           !<-- # of Ebndry tasks on child N that recvd H point data
!
        DO NT=1,NCHILD_TASKS
          IF(ASSOCIATED(CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA))THEN
            DEALLOCATE(CHILD_BOUND_H_EAST(N,INDX2)%TASKS(NT)%DATA)
          ENDIF
        ENDDO
!
        DEALLOCATE(CHILD_BOUND_H_EAST(N,INDX2)%TASKS)
!
        IF(ASSOCIATED(HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV))THEN
          DEALLOCATE(HANDLE_H_EAST(N,INDX2)%NTASKS_TO_RECV)
        ENDIF
!
      ENDIF east_h
!
      NCHILD_TASKS=UBOUND(PD_B_EAST(N)%TASKS,1)                            !<-- # of Ebndry tasks on child N that recvd H point data
      DO NT=1,NCHILD_TASKS
        IF(ASSOCIATED(PD_B_EAST(N)%TASKS(NT)%DATA))THEN
          DEALLOCATE(PD_B_EAST(N)%TASKS(NT)%DATA)
        ENDIF
        DEALLOCATE(FIS_CHILD_EAST(N)%TASKS(NT)%DATA)
      ENDDO
!
      DEALLOCATE(PD_B_EAST(N)%TASKS)
      DEALLOCATE( T_B_EAST(N)%TASKS)
      DEALLOCATE( Q_B_EAST(N)%TASKS)
      DEALLOCATE(CW_B_EAST(N)%TASKS)
      DEALLOCATE(FIS_CHILD_EAST(N)%TASKS)
      DEALLOCATE(WORDS_BOUND_H_EAST(N)%TASKS)
!
      NMAX=UBOUND(HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV,1)
      IF(NMAX>0)THEN
        DO NT=1,NMAX
          CALL MPI_WAIT(HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV(NT)       &  !<-- Handle for ISend from parent task to child N's task NT
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDDO
      ENDIF
!
      east_v: IF(ASSOCIATED(CHILD_BOUND_V_EAST(N,INDX2)%TASKS))THEN
!
        NCHILD_TASKS=UBOUND(CHILD_BOUND_V_EAST(N,INDX2)%TASKS,1)           !<-- # of Ebndry tasks on child N that recvd V point data
!
        DO NT=1,NCHILD_TASKS
          IF(ASSOCIATED(CHILD_BOUND_V_EAST(N,INDX2)%TASKS(NT)%DATA))THEN
            DEALLOCATE(CHILD_BOUND_V_EAST(N,INDX2)%TASKS(NT)%DATA)
          ENDIF
        ENDDO
!
        DEALLOCATE(CHILD_BOUND_V_EAST(N,INDX2)%TASKS)
!
        IF(ASSOCIATED(HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV))THEN
          DEALLOCATE(HANDLE_V_EAST(N,INDX2)%NTASKS_TO_RECV)
        ENDIF
!
      ENDIF east_v
!
      NCHILD_TASKS=UBOUND(PD_B_EAST_V(N)%TASKS,1)                          !<-- # of Ebndry tasks on child N that recvd V point data
      DO NT=1,NCHILD_TASKS
        IF(ASSOCIATED(PD_B_EAST_V(N)%TASKS(NT)%DATA))THEN
          DEALLOCATE(PD_B_EAST_V(N)%TASKS(NT)%DATA)
        ENDIF
      ENDDO
!
      DEALLOCATE(WORDS_BOUND_V_EAST(N)%TASKS)
      DEALLOCATE(PD_B_EAST_V(N)%TASKS)
      DEALLOCATE(U_B_EAST(N)%TASKS)
      DEALLOCATE(V_B_EAST(N)%TASKS)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DEALLOC_WORK_PARENTS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE DEALLOC_WORK_CHILDREN
!
!-----------------------------------------------------------------------
!***  Children deallocate all working pointers that need to be allocated 
!***  with unique dimensions at the outset of the forecast and again
!***  each time a nest moves.  These allocations took place in
!***  subroutine CHILD_RECVS_CHILD_DATA_LIMITS.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: ID_DOM,LIM_HI,N
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      ID_DOM=ID_PARENTS(MY_DOMAIN_ID)                                      !<-- Domain ID of this child's parent
      LIM_HI=FTASKS_DOMAIN(ID_DOM)                                         !<-- # of fcst tasks on this nest's parent domain
!
      DO N=1,LIM_HI
        IF(ALLOCATED(PARENT_TASK(N)%SOUTH_H%STRING))THEN
          DEALLOCATE(PARENT_TASK(N)%SOUTH_H%STRING)                        !<-- Sboundary H datastring from parent 'N'
        ENDIF
        IF(ALLOCATED(PARENT_TASK(N)%SOUTH_V%STRING))THEN
          DEALLOCATE(PARENT_TASK(N)%SOUTH_V%STRING)                        !<-- Sboundary V datastring from parent 'N'
        ENDIF
!
        IF(ALLOCATED(PARENT_TASK(N)%NORTH_H%STRING))THEN
          DEALLOCATE(PARENT_TASK(N)%NORTH_H%STRING)                        !<-- Nboundary H datastring from parent 'N'
        ENDIF
        IF(ALLOCATED(PARENT_TASK(N)%NORTH_V%STRING))THEN
          DEALLOCATE(PARENT_TASK(N)%NORTH_V%STRING)                        !<-- Nboundary V datastring from parent 'N'
        ENDIF
!
        IF(ALLOCATED(PARENT_TASK(N)%WEST_H%STRING))THEN
          DEALLOCATE(PARENT_TASK(N)%WEST_H%STRING)                         !<-- Wboundary H datastring from parent 'N'
        ENDIF
        IF(ALLOCATED(PARENT_TASK(N)%WEST_V%STRING))THEN
          DEALLOCATE(PARENT_TASK(N)%WEST_V%STRING)                         !<-- Wboundary V datastring from parent 'N'
        ENDIF
!
        IF(ALLOCATED(PARENT_TASK(N)%EAST_H%STRING))THEN
          DEALLOCATE(PARENT_TASK(N)%EAST_H%STRING)                         !<-- Eboundary H datastring from parent 'N'
        ENDIF
        IF(ALLOCATED(PARENT_TASK(N)%EAST_V%STRING))THEN
          DEALLOCATE(PARENT_TASK(N)%EAST_V%STRING)                         !<-- Eboundary V datastring from parent 'N'
        ENDIF
!
      ENDDO
!
      IF(NUM_PARENT_TASKS_SENDING_H%SOUTH>0)THEN                           !<-- Did this child task recv Sboundary H data from parent?
        DEALLOCATE(PDB_S)
        DEALLOCATE( TB_S)
        DEALLOCATE( QB_S)
        DEALLOCATE(CWB_S)
        DEALLOCATE(BOUND_1D_SOUTH_H)
      ENDIF
      IF(NUM_PARENT_TASKS_SENDING_V%SOUTH>0)THEN                           !<-- Did this child task recv Sboundary V data from parent?
        DEALLOCATE(UB_S)
        DEALLOCATE(VB_S)
        DEALLOCATE(BOUND_1D_SOUTH_V)
      ENDIF
!
      IF(NUM_PARENT_TASKS_SENDING_H%NORTH>0)THEN                           !<-- Did this child task recv Nboundary H data from parent?
        DEALLOCATE(PDB_N)
        DEALLOCATE( TB_N)
        DEALLOCATE( QB_N)
        DEALLOCATE(CWB_N)
        DEALLOCATE(BOUND_1D_NORTH_H)
      ENDIF
      IF(NUM_PARENT_TASKS_SENDING_V%NORTH>0)THEN                           !<-- Did this child task recv Nboundary V data from parent?
        DEALLOCATE(UB_N)
        DEALLOCATE(VB_N)
        DEALLOCATE(BOUND_1D_NORTH_V)
      ENDIF
!
      IF(NUM_PARENT_TASKS_SENDING_H%WEST>0)THEN                            !<-- Did this child task recv Wboundary H data from parent?
        DEALLOCATE(PDB_W)
        DEALLOCATE( TB_W)
        DEALLOCATE( QB_W)
        DEALLOCATE(CWB_W)
        DEALLOCATE(BOUND_1D_WEST_H)
      ENDIF
      IF(NUM_PARENT_TASKS_SENDING_V%WEST>0)THEN                            !<-- Did this child task recv Wboundary V data from parent?
        DEALLOCATE(UB_W)
        DEALLOCATE(VB_W)
        DEALLOCATE(BOUND_1D_WEST_V)
      ENDIF
!
      IF(NUM_PARENT_TASKS_SENDING_H%EAST>0)THEN                            !<-- Did this child task recv Eboundary H data from parent?
        DEALLOCATE(PDB_E)
        DEALLOCATE( TB_E)
        DEALLOCATE( QB_E)
        DEALLOCATE(CWB_E)
        DEALLOCATE(BOUND_1D_EAST_H)
      ENDIF
      IF(NUM_PARENT_TASKS_SENDING_V%EAST>0)THEN                            !<-- Did this child task recv Eboundary V data from parent?
        DEALLOCATE(UB_E)
        DEALLOCATE(VB_E)
        DEALLOCATE(BOUND_1D_EAST_V)
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DEALLOC_WORK_CHILDREN
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE COMPUTE_STORM_MOTION(NTIMESTEP                         &
                                     ,LAST_STEP_MOVED                   &
                                     ,DT                                &
                                     ,NUM_PES_FCST                      &
                                     ,FIS                               &
                                     ,PD                                &
                                     ,PINT                              &
                                     ,T                                 &
                                     ,Q                                 &
                                     ,CW                                &
                                     ,U                                 &
                                     ,V                                 &
                                     ,DSG2                              &
                                     ,PDSG1                             &
                                     ,DX,DY                             &
                                     ,SEA_MASK                          &
                                     ,I_SW_PARENT_CURRENT               &
                                     ,J_SW_PARENT_CURRENT               &
                                     ,I_WANT_TO_MOVE                    &
                                     ,I_SW_PARENT_NEW                   &
                                     ,J_SW_PARENT_NEW )
!
!-----------------------------------------------------------------------
!***  The nest computes the location of the center of the storm
!***  on its grid and decides if it should move.  This routine
!***  is called at the end of each timestep.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables   
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_SW_PARENT_CURRENT              &  !<-- Parent I of nest domain's SW corner; current position
                                      ,J_SW_PARENT_CURRENT              &  !<-- Parent J of nest domain's SW corner; current position
                                      ,LAST_STEP_MOVED                  &  !<-- Most recent timestep the nest moved
                                      ,NTIMESTEP                        &  !<-- The nest's current timestep
                                      ,NUM_PES_FCST                        !<-- # of forecast tasks
!
      INTEGER(kind=KINT),INTENT(OUT) :: I_SW_PARENT_NEW                 &  !<-- Parent I of nest domain's SW corner; new position
                                       ,J_SW_PARENT_NEW                    !<-- Parent J of nest domain's SW corner; new position
!
      REAL(kind=KFPT),INTENT(IN) :: DT                                  &  !<-- This domain's timestep
                                   ,DY                                     !<-- Delta Y (m) on the nest grid
!
      REAL(kind=KFPT),DIMENSION(JDS:JDE),INTENT(IN) :: DX                  !<-- Delta X (m) on the nest grid
!
      REAL(kind=KFPT),DIMENSION(1:LM),INTENT(IN) :: DSG2,PDSG1
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: FIS      &  !<-- Sfc geopotential (m2/s2)
                                                              ,PD       &  !<-- Psfc minus PTOP (Pa)
                                                              ,SEA_MASK    !<-- Sea mask (1->water)
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM),INTENT(IN) :: CW  &  !<-- Cloud condensate (kg/kg)
                                                                   ,Q   &  !<-- Specific humidity (kg/kg)
                                                                   ,T   &  !<-- Sensible temperature (K)
                                                                   ,U   &  !<-- U component of wind (m/s)
                                                                   ,V      !<-- V component of wind (m/s)
!
      REAL(kind=KFPT),DIMENSION(IMS:IME,JMS:JME,1:LM+1),INTENT(IN) ::   &
                                                                   PINT    !<-- Layer interface pressures (Pa)
!
      LOGICAL(kind=KLOG),INTENT(OUT) :: I_WANT_TO_MOVE                     !<-- Does nest want to move to a new position?
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT),SAVE :: I_EAST,I_WEST                          &  !<-- These define the nest domain search limits
                                ,J_NORTH,J_SOUTH                           !    within the central window.
!
      INTEGER(kind=KINT),SAVE :: ID_DUMMY=-999                          &  !<-- Dummy value for task ID
                                ,ITAG_PG=200                               !<-- Hardwire this tag for Sends/Recvs of Pgrad values
!
      INTEGER(kind=KINT),SAVE :: I_MAX,I_MIN                            &
                                ,J_MAX,J_MIN                            &
                                ,NPTS_NS,NPTS_WE
!
      INTEGER(kind=KINT) :: I,I_CENTER_CURRENT,I_CENTER_NEW             &
                           ,I_DIFF,ID_PE_MIN                            &
                           ,J,J_CENTER_CURRENT,J_CENTER_NEW             &
                           ,J_DIFF                                      &
                           ,L,L1,L2,L3,N,N_SEND
!
      INTEGER(kind=KINT) :: IERR,ISTAT
!
      INTEGER(kind=KINT),DIMENSION(1:4),SAVE :: I_PG,J_PG
!
      INTEGER(kind=KINT),DIMENSION(1:4) :: HANDLE_PVAL
!
      INTEGER(kind=KINT),DIMENSION(0:NUM_PES_FCST-1) :: HANDLE_PDYN     &
                                                       ,HANDLE_WIN
!
      INTEGER(kind=KINT),DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(kind=KFPT),SAVE :: DIST_GRAD=100000.                         &  !<-- Distance (m) for checking storm-scale pressure gradient
                             ,DIST_LAND=100000.                         &  !<-- Distance (m) for checking shift over land
                             ,ELAPSED_TIME_MAX=30000.                   &  !<-- Maximum time (sec) between nest shifts over land
                             ,GI=1./G                                   & 
                             ,PGRD_MIN=200.                             &  !<-- Minimum storm scale pressure gradient (Pa)
                             ,PGRD_MIN_LAND=400.                        &  !<-- Minimum storm scale pressure gradient (Pa) over land
                             ,STD_LAPSE=6.5E-3                          &  !<-- Standard atmospheric lapse rate
                             ,THIRD=1./3.                               &
                             ,Z1=2000.                                  &  !<-- In computing the dynamic pressure use the winds
                             ,Z2=1500.                                  &  !    at around 2km, 1.5km and 1km.
                             ,Z3=1000.                                     !<--
!
      REAL(kind=KFPT),SAVE :: COEF                                      &
                             ,ELAPSED_TIME_MIN                          &
                             ,RNPTS_HZ
                         
!
      REAL(kind=KFPT) :: APELP,DFDP,DIST_I,DIST_J,DZ                    &
                        ,ELAPSED_TIME,FACTOR,FRAC_SEA                   &
                        ,PARENT_DIFF,PCHECK                             &
                        ,PDYN_MIN_GBL,PMAX,PVAL                         &
                        ,SLP_MIN,SQWS,SUM_SEA,TSFC                      &
                        ,ZDIFF,ZLOW,ZSAVE1,ZSAVE2,ZSAVE3
!
      REAL(kind=KFPT),DIMENSION(:),POINTER :: PDYN_MIN_VALS                !<-- Minimum PDYN on each task subdomain
!
      REAL(kind=KFPT),DIMENSION(1:3,0:NUM_PES_FCST-1),TARGET ::         &
                                                             PDYN_MIN      !<-- Minimum PDYN and its I,J on each nest task subdomain
!
      REAL(kind=KFPT),DIMENSION(1:LM) :: ZMEAN
!
      REAL(kind=KFPT),DIMENSION(ITS:ITE,JTS:JTE) :: PDYN,SLP
!
      REAL(kind=KFPT),DIMENSION(ITS:ITE,JTS:JTE,1:LM+1) :: Z
!
      LOGICAL(kind=KLOG),SAVE :: FIRST_PASS=.TRUE.                      &
                                ,I_HOLD_CENTER_POINT
!
      LOGICAL(kind=KLOG),DIMENSION(:),ALLOCATABLE,SAVE :: IN_WINDOW
!
      LOGICAL(kind=KLOG),DIMENSION(1:4) :: NO_VALUE
!
      LOGICAL(kind=KLOG),DIMENSION(1:4),SAVE :: I_HOLD_PG_POINT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      I_WANT_TO_MOVE=.FALSE.                                               !<-- Begin with this assumption
!
!-----------------------------------------------------------------------
!***  Allow the nest to move only in one of its physics timesteps.
!-----------------------------------------------------------------------
!
      IF(MOD(NTIMESTEP,NPHS)/=0)THEN
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
!***  What fraction of the entire nest domain is now over water?
!-----------------------------------------------------------------------
!
      SUM_SEA=0.
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        SUM_SEA=SUM_SEA+SEA_MASK(I,J)
      ENDDO
      ENDDO
!
      CALL MPI_ALLREDUCE(SUM_SEA                                        &  !<-- The value on each task
                        ,FRAC_SEA                                       &  !<-- The final value on all tasks
                        ,1                                              &  !<-- Operating on single words
                        ,MPI_REAL                                       &  !<-- Datatype
                        ,MPI_SUM                                        &  !<-- Sum all the values
                        ,MPI_COMM_COMP                                  &  !<-- MPI communicator
                        ,IERR )
!
      FRAC_SEA=FRAC_SEA/(IDE*JDE)
!
!-----------------------------------------------------------------------
!***  Do not consider moving if the storm has been stationary 
!***  too long over land.
!-----------------------------------------------------------------------
!
      IF(FRAC_SEA<=0.2.AND.ELAPSED_TIME>ELAPSED_TIME_MAX)THEN
        IF(MYPE==0)THEN
          WRITE(0,*)' SKIP MOTION COMPUTATION:  Storm stationary over land.'
        ENDIF
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
!***  Which point is the current center of the moving nest domain?
!***  This is the most recent storm center point.
!-----------------------------------------------------------------------
!
      I_CENTER_CURRENT=IDS+INT(0.5*(IDE-IDS)+EPS)     
      J_CENTER_CURRENT=JDS+INT(0.5*(JDE-JDS)+EPS)    
!
!-----------------------------------------------------------------------
!***  In the first pass through this routine the nest tasks
!***  determine which of their subdomains contain points within
!***  the nest domain's central search window over the storm.  Then
!***  they exchange that information with each other.
!-----------------------------------------------------------------------
!
      prelim: IF(FIRST_PASS)THEN                                           !<-- Work needs to be done only once during the forecast
!
!-----------------------------------------------------------------------
!
        FIRST_PASS=.FALSE.
!
        I_WEST =IDE/2-IDE/3                                                !<-- These define the nest domain search limits
        I_EAST =IDE/2+IDE/3                                                !    within the central window.
        J_SOUTH=JDE/2-JDE/3                                                !
        J_NORTH=JDE/2+JDE/3                                                !<--
!
        IF(.NOT.ALLOCATED(IN_WINDOW))THEN
          ALLOCATE(IN_WINDOW(0:NUM_PES_FCST-1),stat=ISTAT)                 !<-- Must be allocated since 0:NUM_PES_FCST is passed in
          IF(ISTAT/=0)THEN
            WRITE(0,*)' COMPUTE_STORM_MOTION: ERROR'
            WRITE(0,*)' Failed to allocate IN_WINDOW stat=',ISTAT
          ENDIF
        ENDIF
!
        IN_WINDOW(MYPE)=.FALSE.
!
        IF(ITS<=I_EAST.AND.ITE>=I_WEST                                  &  !<-- Does any of nest task N's subdomain
                      .AND.                                             &  !    lie within the central window?
           JTS<=J_NORTH.AND.JTE>=J_SOUTH)THEN
!
          IN_WINDOW(MYPE)=.TRUE.                                           !<-- Yes, this task lies in the search window
!
          I_MIN=MAX(ITS,I_WEST)                                            !<-- Index limits of this task's subdomain
          I_MAX=MIN(ITE,I_EAST)                                            !    that lie inside the search window.
          J_MIN=MAX(JTS,J_SOUTH)                                           !
          J_MAX=MIN(JTE,J_NORTH)                                           !<--
!
          RNPTS_HZ=1./REAL((I_MAX-I_MIN+1)*(J_MAX-J_MIN+1))                !<-- Reciprocal of # of task's points in search window
!
        ENDIF
!
        DO N=0,NUM_PES_FCST-1                                              !<-- Nest fcst tasks send their window status to each other
          IF(N/=MYPE)THEN                                                  !<-- But not to themselves
            CALL MPI_ISEND(IN_WINDOW(MYPE)                              &  !<-- This task's central window status
                          ,1                                            &  !<-- Sending this many words
                          ,MPI_LOGICAL                                  &  !<-- Datatype
                          ,N                                            &  !<-- Sending to this local nest task ID
                          ,MYPE                                         &  !<-- Use this task's ID as the tag
                          ,MPI_COMM_COMP                                &  !<-- MPI communicator
                          ,HANDLE_WIN(N)                                &  !<-- Handle for task N's Recv
                          ,IERR )
          ENDIF
        ENDDO
!
        DO N=0,NUM_PES_FCST-1                                              !<-- Nest fcst tasks recv window status from each other
          IF(N/=MYPE)THEN                                                  !<-- But not from themselves
            CALL MPI_RECV(IN_WINDOW(N)                                  &  !<-- Nest task N's central window status
                         ,1                                             &  !<-- Receiving this many words
                         ,MPI_LOGICAL                                   &  !<-- Datatype
                         ,N                                             &  !<-- Data was sent by this nest task
                         ,N                                             &  !<-- Tag is the sender's rank
                         ,MPI_COMM_COMP                                 &  !<-- MPI communicator
                         ,JSTAT                                         &  !<-- MPI status
                         ,IERR )
          ENDIF
        ENDDO
!
        DO N=0,NUM_PES_FCST-1
          IF(N/=MYPE)THEN   
            CALL MPI_WAIT(HANDLE_WIN(N)                                 &  !<-- Proceed only after all Recvs have completed
                         ,JSTAT                                         & 
                         ,IERR )
          ENDIF
        ENDDO
!
!-----------------------------------------------------------------------
!***  We will need an approximation of the maximum pressure gradient 
!***  around the storm center.  Use SLP values roughly DIST_GRAD meters
!***  N, S, W, and E of the center.  What is the relative distance
!***  in gridpoints that this distance represents?
!-----------------------------------------------------------------------
!
        NPTS_NS=NINT(DIST_GRAD/DY)                                         !<-- # grid points to N/S of center to check SLP
        NPTS_WE=NINT(DIST_GRAD/DX(J_CENTER_CURRENT))                       !<-- # grid points to W/E of center to check SLP
!
!-----------------------------------------------------------------------
!
        COEF=-G/(R_D*STD_LAPSE)
!
!-----------------------------------------------------------------------
!***  Set the minimum time between domain shifts.  For now make it a
!***  linear relationship where the 9km nest must wait at least 15 min
!***  to move.
!-----------------------------------------------------------------------
!
        ELAPSED_TIME_MIN=45.*DT                                            !<-- Equals 900s for a 9 km nest.
!
!-----------------------------------------------------------------------
!***  Allow the minimum pressure gradient between the storm center
!***  and the cardinal points to be smaller for higher resolution.
!-----------------------------------------------------------------------
!
        IF(DY<8999.)THEN
          PGRD_MIN=PGRD_MIN*DY/9000.
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF prelim
!
!-----------------------------------------------------------------------
!***  Never let the nest move until a minimum amount of time
!***  has passed.
!-----------------------------------------------------------------------
!
      ELAPSED_TIME=(NTIMESTEP-LAST_STEP_MOVED)*DT
!
      IF(ELAPSED_TIME<ELAPSED_TIME_MIN)THEN
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
!
      DO N=0,NUM_PES_FCST-1
        PDYN_MIN(1,N)=110000.                                              !<-- Initialize the minimum value of PDYN.
        PDYN_MIN(2,N)=-1000.                                               !<-- Initialize the I index of minimum PDYN
        PDYN_MIN(3,N)=-1000.                                               !<-- Initialize the J index of minimum PDYN
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        SLP(I,J)=0.                                                        !<-- Initialize the SLP array.
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Begin the search for the new storm center location.
!***  Search inside the central window.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      window: IF(IN_WINDOW(MYPE))THEN                                      !<-- Only those tasks within the search window proceed
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Compute a 'dynamic pressure' using vertically averaged winds
!***  around Z1, Z2, and Z3 meters as well as the sea level pressure.
!-----------------------------------------------------------------------
!
        DO J=J_MIN,J_MAX
        DO I=I_MIN,I_MAX
          Z(I,J,LM+1)=FIS(I,J)*GI                                          !<-- The surface elevation (m)
        ENDDO
        ENDDO
!
        DO L=LM,1,-1
          ZMEAN(L)=0.
          DO J=J_MIN,J_MAX
          DO I=I_MIN,I_MAX
            APELP=(PINT(I,J,L)+PINT(I,J,L+1))*0.5
            DFDP=(Q(I,J,L)*P608+(1.-CW(I,J,L)))*T(I,J,L)*R_D/APELP
            DZ=GI*DFDP*(DSG2(L)*PD(I,J)+PDSG1(L))
            Z(I,J,L)=Z(I,J,L+1)+DZ
            ZMEAN(L)=ZMEAN(L)+0.5*(Z(I,J,L)+Z(I,J,L+1))
          ENDDO
          ENDDO
          ZMEAN(L)=ZMEAN(L)*RNPTS_HZ                                       !<-- The mean height (m) of midlayer L in search window
        ENDDO
!
        ZSAVE1=1.E10
        ZSAVE2=1.E10
        ZSAVE3=1.E10
!
!-----------------------------------------------------------------------
!***  Find and save the model midlayer indices nearest to the 
!***  prescribed heights of Z1, Z2, and Z3 meters above the ground.
!-----------------------------------------------------------------------
!
        DO L=1,LM
!
          ZDIFF=ABS(ZMEAN(L)-Z1)
          IF(ZDIFF<ZSAVE1)THEN
            ZSAVE1=ZDIFF
            L1=L
          ENDIF
!
          ZDIFF=ABS(ZMEAN(L)-Z2)
          IF(ZDIFF<ZSAVE2)THEN
            ZSAVE2=ZDIFF
            L2=L
          ENDIF
!
          ZDIFF=ABS(ZMEAN(L)-Z3)
          IF(ZDIFF<ZSAVE3)THEN
            ZSAVE3=ZDIFF
            L3=L
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Compute a 'dynamic pressure' using the sea level pressure
!***  and the winds nearest the heights Z1, Z2, and Z3.
!-----------------------------------------------------------------------
!
        DO J=J_MIN,J_MAX
        DO I=I_MIN,I_MAX
          ZLOW=0.5*(Z(I,J,LM)+Z(I,J,LM+1))                                 !<-- Height of lowest midlayer
          TSFC=T(I,J,LM)*(Q(I,J,LM)*P608+(1.-CW(I,J,LM)))+STD_LAPSE*ZLOW   !<-- Representative sfc T
          FACTOR=STD_LAPSE*Z(I,J,LM+1)/TSFC
          SLP(I,J)=PINT(I,J,LM+1)*(1.-FACTOR)**COEF
          SQWS=(U(I,J,L1)*U(I,J,L1)+V(I,J,L1)*V(I,J,L1)                 &
               +U(I,J,L2)*U(I,J,L2)+V(I,J,L2)*V(I,J,L2)                 &
               +U(I,J,L3)*U(I,J,L3)+V(I,J,L3)*V(I,J,L3))*THIRD
!         PDYN(I,J)=SLP(I,J)+0.55*SQWS
!
          PDYN(I,J)=SLP(I,J)                                               !<-- FOR NOW, set PDYN to the sea level pressure
!
          IF(PDYN(I,J)<PDYN_MIN(1,MYPE))THEN
            PDYN_MIN(1,MYPE)=PDYN(I,J)                                     !<-- Minimum PDYN on this task subdomain
            PDYN_MIN(2,MYPE)=REAL(I)                                       !<-- I of the local minimum of PDYN
            PDYN_MIN(3,MYPE)=REAL(J)                                       !<-- J of the local minimum of PDYN
          ENDIF
        ENDDO
        ENDDO
!
!-----------------------------------------------------------------------
!
      ENDIF window
!
!-----------------------------------------------------------------------
!***  Those nest tasks with any of their subdomains inside the search 
!***  window share their PDYN values and locations with all other
!***  nest tasks.
!-----------------------------------------------------------------------
!
      IF(IN_WINDOW(MYPE))THEN                                              !<-- Tasks only send if they are in the central window
        DO N=0,NUM_PES_FCST-1                                              !<-- Nest window tasks send their PDYN_MIN to all other tasks
          IF(N/=MYPE)THEN
            CALL MPI_ISEND(PDYN_MIN(1,MYPE)                             &  !<-- This task's minimum PDYN with its I,J
                          ,3                                            &  !<-- Sending this many words
                          ,MPI_REAL                                     &  !<-- Datatype
                          ,N                                            &  !<-- Sending to this local nest task ID
                          ,MYPE                                         &  !<-- Use this task's ID as the tag
                          ,MPI_COMM_COMP                                &  !<-- MPI communicator
                          ,HANDLE_PDYN(N)                               &  !<-- Handle for task N's Recv
                          ,IERR )
          ENDIF
        ENDDO
      ENDIF
!
      DO N=0,NUM_PES_FCST-1                                                !<-- All nest fcst tasks recv PDYN_MIN from window tasks
        IF(N/=MYPE.AND.IN_WINDOW(N))THEN
          CALL MPI_RECV(PDYN_MIN(1,N)                                   &  !<-- Nest task N's minimum PDYN
                       ,3                                               &  !<-- Receiving this many words
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,N                                               &  !<-- Data was sent by this nest task
                       ,N                                               &  !<-- Tag is the sender's rank
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
        ENDIF
      ENDDO
!
      IF(IN_WINDOW(MYPE))THEN                                              !<-- Tasks only send if they are in the central window
        DO N=0,NUM_PES_FCST-1                                              !<-- Nest window tasks send their PDYN_MIN to all other tasks
          IF(N/=MYPE)THEN
            CALL MPI_WAIT(HANDLE_PDYN(N)                                &  !<-- Proceed only after all Recvs have completed
                         ,JSTAT                                         & 
                         ,IERR )
          ENDIF
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!***  All nest tasks determine the location of the minimum value
!***  of PDYN within the central window.  This will be the new
!***  location of the storm center if the decision is made to
!***  go ahead with the move.
!-----------------------------------------------------------------------
!
      PDYN_MIN_VALS=>PDYN_MIN(1,0:NUM_PES_FCST-1)                          !<-- Select only the PDYN values from the array
!
      ID_PE_MIN=MINLOC(PDYN_MIN_VALS,1)-1                                  !<-- ID of nest task with minimum value of PDYN
!
      SLP_MIN=PDYN_MIN(1,ID_PE_MIN)                                        !<-- Minimum value of PDYN for all tasks
      I_CENTER_NEW=INT(PDYN_MIN(2,ID_PE_MIN))                              !<-- I index of minimum PDYN on the nest domain
      J_CENTER_NEW=INT(PDYN_MIN(3,ID_PE_MIN))                              !<-- J index of minimum PDYN on the nest domain
!
      I_DIFF=I_CENTER_NEW-I_CENTER_CURRENT                                 !<-- Shift in I to potential new center
      J_DIFF=J_CENTER_NEW-J_CENTER_CURRENT                                 !<-- Shift in J to potential new center 
!
!-----------------------------------------------------------------------
!***  If the nest moves then its SW corner must shift from one parent
!***  H point to another which means the I and J shifts must be in
!***  integer multiples of SPACE_RATIO_MY_PARENT.  Adjust I_DIFF and
!***  J_DIFF given this constraint.
!-----------------------------------------------------------------------
!
      IF(MOD(I_DIFF,SPACE_RATIO_MY_PARENT)/=0)THEN
        PARENT_DIFF=REAL(I_DIFF)/REAL(SPACE_RATIO_MY_PARENT)
        IF(ABS(FRACTION(PARENT_DIFF))>0.5)THEN
          I_DIFF=NINT(PARENT_DIFF)*SPACE_RATIO_MY_PARENT
        ELSE
          I_DIFF=INT(PARENT_DIFF)*SPACE_RATIO_MY_PARENT
        ENDIF
        I_CENTER_NEW=I_CENTER_CURRENT+I_DIFF
      ENDIF
!
      IF(MOD(J_DIFF,SPACE_RATIO_MY_PARENT)/=0)THEN
        PARENT_DIFF=REAL(J_DIFF)/REAL(SPACE_RATIO_MY_PARENT)
        IF(ABS(FRACTION(PARENT_DIFF))>0.5)THEN
          J_DIFF=NINT(PARENT_DIFF)*SPACE_RATIO_MY_PARENT
        ELSE
          J_DIFF=INT(PARENT_DIFF)*SPACE_RATIO_MY_PARENT
        ENDIF
        J_CENTER_NEW=J_CENTER_CURRENT+J_DIFF
      ENDIF
!
      IF(ABS(I_DIFF)==0.AND.ABS(J_DIFF)==0)THEN
        IF(MYPE==0)THEN
          WRITE(0,*)' NO MOTION: Less than one parent grid increment.'
        ENDIF
        RETURN                                                             !<-- No motion so exit.
      ENDIF
!
!-----------------------------------------------------------------------
!***  Tasks learn if they contain any of the four cardinal direction
!***  points to be used for the storm-gradient check.
!-----------------------------------------------------------------------
!
      N_SEND=4
      DO N=1,4
        NO_VALUE(N)=.FALSE.                                                !<-- Start with all cardinal directions having values inside
      ENDDO
!
      I_PG(1)=I_CENTER_NEW                                                 !<-- I coordinate of north pressure point
      J_PG(1)=J_CENTER_NEW+NPTS_NS                                         !<-- J coordinate of north pressure point
      IF(I_PG(1)<IDS.OR.I_PG(1)>IDE.OR.J_PG(1)<JDS.OR.J_PG(1)>JDE)THEN
        NO_VALUE(1)=.TRUE.                                                 !<-- North point is outside nest domain
        N_SEND=N_SEND-1
      ELSE
        CALL LOCATE_POINT_ON_TASKS(I_PG(1)                              &
                                  ,J_PG(1)                              &
                                  ,I_HOLD_PG_POINT(1))                     !<-- Does this task subdomain contain north PG point?
      ENDIF
!
      I_PG(2)=I_CENTER_NEW                                                 !<-- I coordinate of south pressure point
      J_PG(2)=J_CENTER_NEW-NPTS_NS                                         !<-- J coordinate of south pressure point
      IF(I_PG(2)<IDS.OR.I_PG(2)>IDE.OR.J_PG(2)<JDS.OR.J_PG(2)>JDE)THEN
        NO_VALUE(2)=.TRUE.                                                 !<-- South point is outside nest domain
        N_SEND=N_SEND-1
      ELSE
        CALL LOCATE_POINT_ON_TASKS(I_PG(2)                              &
                                  ,J_PG(2)                              &
                                  ,I_HOLD_PG_POINT(2))                     !<-- Does this task subdomain contain south PG point?
      ENDIF
!
      I_PG(3)=I_CENTER_NEW-NPTS_WE                                         !<-- I coordinate of west pressure point
      J_PG(3)=J_CENTER_NEW                                                 !<-- J coordinate of west pressure point
      IF(I_PG(3)<IDS.OR.I_PG(3)>IDE.OR.J_PG(3)<JDS.OR.J_PG(3)>JDE)THEN
        NO_VALUE(3)=.TRUE.                                                 !<-- West point is outside nest domain
        N_SEND=N_SEND-1
      ELSE
        CALL LOCATE_POINT_ON_TASKS(I_PG(3)                              &
                                  ,J_PG(3)                              &
                                  ,I_HOLD_PG_POINT(3))                     !<-- Does this task subdomain contain west PG point?
      ENDIF
!
      I_PG(4)=I_CENTER_NEW+NPTS_WE                                         !<-- I coordinate of east pressure point
      J_PG(4)=J_CENTER_NEW                                                 !<-- J coordinate of east pressure point
      IF(I_PG(4)<IDS.OR.I_PG(4)>IDE.OR.J_PG(4)<JDS.OR.J_PG(4)>JDE)THEN
        NO_VALUE(4)=.TRUE.                                                 !<-- East point is outside nest domain
        N_SEND=N_SEND-1
      ELSE
        CALL LOCATE_POINT_ON_TASKS(I_PG(4)                              &
                                  ,J_PG(4)                              &
                                  ,I_HOLD_PG_POINT(4))                     !<-- Does this task subdomain contain east PG point?
      ENDIF
!
!-----------------------------------------------------------------------
!***  Those tasks that hold the four points N, S, W, and E of the
!***  new storm center send their pressure values to the task
!***  whose subdomain contains the new storm center.
!-----------------------------------------------------------------------
!
      IF(N_SEND>0)THEN
        DO N=1,4                                                           !<-- The 4 cardinal direction points
          IF(NO_VALUE(N))CYCLE
          IF(I_HOLD_PG_POINT(N))THEN
            PVAL=SLP(I_PG(N),J_PG(N))
            CALL MPI_ISEND(PVAL                                         &  !<-- SLP at cardinal point N around storm center
                          ,1                                            &  !<-- It is one word
                          ,MPI_REAL                                     &  !<-- Datatype
                          ,ID_PE_MIN                                    &  !<-- ID of task holding the new storm center
                          ,ITAG_PG                                      &  !<-- Tag used for exchange of pressure data for gradient
                          ,MPI_COMM_COMP                                &  !<-- MPI communicator
                          ,HANDLE_PVAL(N)                               &  !<-- Handle for task N's Recv
                          ,IERR )
          ENDIF
        ENDDO
      ELSE
        WRITE(0,*)' ALERT: Storm has moved more than DIST_GRAD beyond'  &
                 ,' moving nest boundary!'
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
!***  The task holding the new storm center saves the maximum
!***  pressure among the four cardinal points then informs all 
!***  other tasks of that value.
!-----------------------------------------------------------------------
!
      IF(MYPE==ID_PE_MIN)THEN
        PMAX=-100000.
        DO N=1,4
          IF(NO_VALUE(N))CYCLE
          CALL MPI_RECV(PVAL                                            &  !<-- Pressure from cardinal point N
                       ,1                                               &  !<-- It is one word
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,MPI_ANY_SOURCE                                  &  !<-- Does not know ID of the sending task
                       ,ITAG_PG                                         &  !<-- Tag used for exchange of pressure data for gradient
                       ,MPI_COMM_COMP                                   &  !<-- MPI communicator
                       ,JSTAT                                           &  !<-- MPI status
                       ,IERR )
!
          IF(PVAL>PMAX)THEN
            PMAX=PVAL                                                      !<-- Save the maximum pressure.
          ENDIF
!
        ENDDO
      ENDIF
!
      DO N=1,4                                                             !<-- The 4 cardinal directions
        IF(NO_VALUE(N))CYCLE
        IF(I_HOLD_PG_POINT(N))THEN
          CALL MPI_WAIT(HANDLE_PVAL(N)                                  &  !<-- Proceed only after all Recvs have completed
                       ,JSTAT                                           & 
                       ,IERR )
        ENDIF
      ENDDO
!
      IF(N_SEND>0)THEN
        CALL MPI_BCAST(PMAX                                             &  !<-- Max cardinal pressure around storm
                      ,1                                                &  !<-- It is one word
                      ,MPI_REAL                                         &  !<-- Datatype
                      ,ID_PE_MIN                                        &  !<-- The root sender
                      ,MPI_COMM_COMP                                    &  !<-- MPI communicator
                      ,IERR )
        IF(PMAX<50000.)PMAX=1000000.
      ELSE
        PMAX=1000000.
      ENDIF
!
!-----------------------------------------------------------------------
!***  Now that a new central pressure and its location have been 
!***  identified several conditions must be met before the move
!***  can be allowed to execute.  We have already made certain
!***  that the nest domain's SW corner remains on a parent H point
!***  if the shift is executed.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  (1) Did the vortex disappear?
!-----------------------------------------------------------------------
!
      IF(PMAX<SLP_MIN)THEN
        IF(MYPE==0)THEN
          WRITE(0,*)' CANNOT MOVE:  Lost the vortex.'
        ENDIF
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
!***  (2) Is the storm-scale pressure gradient too weak? 
!-----------------------------------------------------------------------
!
      IF(PMAX-SLP_MIN<PGRD_MIN)THEN
        IF(MYPE==0)THEN
          WRITE(0,*)' DO NOT MOVE:  Pressure gradient too weak.'
        ENDIF
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
!***  (3) Has the vortex been lost over land?
!-----------------------------------------------------------------------
!
      DIST_I=I_DIFF*DX(J_CENTER_NEW)                                       !<-- Shift distance (m) in I
      DIST_J=J_DIFF*DY                                                     !<-- Shift distance (m) in J
!
      IF(FRAC_SEA<=0.2.AND.DIST_I>DIST_LAND.AND.DIST_J>DIST_LAND)THEN
        IF(MYPE==0)THEN
          WRITE(0,*)' CANNOT MOVE:  Vortex lost over land.'
        ENDIF
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
!***  (4) Is the vortex too weak over land?
!-----------------------------------------------------------------------
!
      IF(FRAC_SEA<=0.2.AND.PMAX-SLP_MIN<PGRD_MIN_LAND)THEN
        IF(MYPE==0)THEN
          WRITE(0,*)' DO NOT MOVE:  Vortex is too weak over land.'
        ENDIF
        RETURN
      ENDIF
!
!-----------------------------------------------------------------------
!***  The nest mechanics are anchored on the parent's I,J of the
!***  nest's SW corner (H point).  Given the new center point
!***  what is the parent I,J of the nest's SW corner after it moves?
!-----------------------------------------------------------------------
!
      I_WANT_TO_MOVE=.TRUE.
!
      I_SW_PARENT_NEW=I_SW_PARENT_CURRENT+I_DIFF/SPACE_RATIO_MY_PARENT
      J_SW_PARENT_NEW=J_SW_PARENT_CURRENT+J_DIFF/SPACE_RATIO_MY_PARENT
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE LOCATE_POINT_ON_TASKS(I_PG                             &
                                      ,J_PG                             &
                                      ,I_HOLD_THIS_POINT )
!
!-----------------------------------------------------------------------
!***  Find the ID of the task whose subdomain contains 
!***  point (I_PG,J_PG).
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: I_PG                             &  !<-- I index of point in question
                                      ,J_PG                                !<-- J index of point in question
!
      LOGICAL(kind=KLOG),INTENT(OUT) :: I_HOLD_THIS_POINT                  !<-- This task subdomain holds the point?
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      I_HOLD_THIS_POINT=.FALSE.
!
      IF(I_PG>=ITS.AND.I_PG<=ITE.AND.J_PG>=JTS.AND.J_PG<=JTE)THEN
        I_HOLD_THIS_POINT=.TRUE.
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE LOCATE_POINT_ON_TASKS
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE COMPUTE_STORM_MOTION
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
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: NNEW,NOLD
!
      REAL(kind=KFPT),DIMENSION(1:NOLD),INTENT(IN) :: XOLD,YOLD
      REAL(kind=KFPT),DIMENSION(1:NNEW),INTENT(IN) :: XNEW
!
      REAL(kind=KFPT),DIMENSION(1:NOLD),INTENT(INOUT) :: Y2
!
      REAL(kind=KFPT),DIMENSION(1:NNEW),INTENT(OUT) :: YNEW
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: K,K1,K2,KOLD,NOLDM1
!
      REAL(kind=KFPT) :: AK,BK,CK,DEN,DX,DXC,DXL,DXR,DYDXL,DYDXR        &
                        ,RDX,RTDXC,X,XK,XSQ,Y2K,Y2KP1
!
      REAL(kind=KFPT),DIMENSION(1:NOLD-2) :: P,Q
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
!
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
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine artificial_move(ntimestep                              &
                                ,kount_moves                            &
                                ,i_want_to_move                         &
                                ,i_sw_parent_current                    &
                                ,j_sw_parent_current                    &
                                ,i_sw_parent_new                        &
                                ,j_sw_parent_new )
!
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: i_sw_parent_current              &
                                      ,j_sw_parent_current              &
                                      ,kount_moves                      &
                                      ,ntimestep
!
      integer(kind=kint),intent(out) :: i_sw_parent_new                 &
                                       ,j_sw_parent_new
!
      logical(kind=klog),intent(out) :: i_want_to_move
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      i_want_to_move=.false.
!
      write(0,*)' enter artificial ntimestep=',ntimestep,' kount_moves=',kount_moves
      if(ntimestep>0.and.mod(ntimestep,51)==0)then
        i_want_to_move=.true.   
      write(0,*)' artificial set i_want_to_move=',i_want_to_move
        if(mod(kount_moves,16)<= 1)then                                    !<-- NW
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to NW'
        elseif(mod(kount_moves,16)<= 3)then                                !<-- N
          i_sw_parent_new=i_sw_parent_current
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to N'
        elseif(mod(kount_moves,16)<= 5)then                                !<-- NE
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to NE'
        elseif(mod(kount_moves,16)<= 7)then                                !<-- E
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current
      write(0,*)' artificial to E'
        elseif(mod(kount_moves,16)<= 9)then                                !<-- SE
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to SE'
        elseif(mod(kount_moves,16)<=11)then                                !<-- S
          i_sw_parent_new=i_sw_parent_current
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to S'
        elseif(mod(kount_moves,16)<=13)then                                !<-- SW
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to SW'
        elseif(mod(kount_moves,16)<=15)then                                !<-- W
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current
      write(0,*)' artificial to W'
        endif
      endif
      write(0,*)' exit artificial_move ntimestep=',ntimestep,' mod=',mod(ntimestep,14) &
               ,' i_want_to_move=',i_want_to_move
!
!-----------------------------------------------------------------------
      end subroutine artificial_move
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine artificial_move2(ntimestep                              &
                                ,kount_moves                            &
                                ,i_want_to_move                         &
                                ,i_sw_parent_current                    &
                                ,j_sw_parent_current                    &
                                ,i_sw_parent_new                        &
                                ,j_sw_parent_new )
!
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: i_sw_parent_current              &
                                      ,j_sw_parent_current              &
                                      ,kount_moves                      &
                                      ,ntimestep
!
      integer(kind=kint),intent(out) :: i_sw_parent_new                 &
                                       ,j_sw_parent_new
!
      logical(kind=klog),intent(out) :: i_want_to_move
!
      integer(kind=kint) :: mod_km
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      i_want_to_move=.false.
!
      write(0,*)' enter artificial ntimestep=',ntimestep,' kount_moves=',kount_moves
!!!   if(ntimestep>3.and.mod(ntimestep,14)<=2)then
      if(ntimestep>0.and.mod(ntimestep+3,51)==0)then
        i_want_to_move=.true.   
      write(0,*)' artificial set i_want_to_move=',i_want_to_move
        mod_km=mod(kount_moves,16)
        if(mod_km<=01)then                                                 !<-- E
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current
      write(0,*)' artificial to E mod_km=',mod_km
        elseif(mod_km<=03)then                                             !<-- NE
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to NE mod_km=',mod_km
        elseif(mod_km<=05)then                                             !<-- N
          i_sw_parent_new=i_sw_parent_current
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to N mod_km=',mod_km
        elseif(mod_km<=07)then                                             !<-- NW
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to NW mod_km=',mod_km
        elseif(mod_km<=09)then                                             !<-- W
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current
      write(0,*)' artificial to W mod_km=',mod_km
        elseif(mod_km<=11)then                                             !<-- SW
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to SW mod_km=',mod_km
        elseif(mod_km<=13)then                                             !<-- S
          i_sw_parent_new=i_sw_parent_current
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to S mod_km=',mod_km
        elseif(mod_km<=15)then                                             !<-- SE
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to SE mod_km=',mod_km
        endif
      endif
      write(0,*)' exit artificial_move ntimestep=',ntimestep,' mod=',mod(ntimestep+3,51) &
               ,' i_want_to_move=',i_want_to_move
!
!-----------------------------------------------------------------------
      end subroutine artificial_move2
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine artificial_move3(ntimestep                              &
                                ,kount_moves                            &
                                ,i_want_to_move                         &
                                ,i_sw_parent_current                    &
                                ,j_sw_parent_current                    &
                                ,i_sw_parent_new                        &
                                ,j_sw_parent_new )
!
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: i_sw_parent_current              &
                                      ,j_sw_parent_current              &
                                      ,kount_moves                      &
                                      ,ntimestep
!
      integer(kind=kint),intent(out) :: i_sw_parent_new                 &
                                       ,j_sw_parent_new
!
      logical(kind=klog),intent(out) :: i_want_to_move
!
      integer(kind=kint) :: mod_km
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      i_want_to_move=.false.
!
      write(0,*)' enter artificial ntimestep=',ntimestep,' kount_moves=',kount_moves
!!!   if(ntimestep>3.and.mod(ntimestep,14)<=2)then
      if(ntimestep>0.and.mod(ntimestep+3,51)==0)then
        i_want_to_move=.true.   
      write(0,*)' artificial set i_want_to_move=',i_want_to_move
        mod_km=mod(kount_moves,16)
        if(mod_km<=01)then                                                 !<-- NW
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to NW mod_km=',mod_km
        elseif(mod_km<=03)then                                             !<-- N
          i_sw_parent_new=i_sw_parent_current
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to N mod_km=',mod_km
        elseif(mod_km<=05)then                                             !<-- NE
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to NE mod_km=',mod_km
        elseif(mod_km<=07)then                                             !<-- E
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current
      write(0,*)' artificial to E mod_km=',mod_km
        elseif(mod_km<=09)then                                             !<-- SE
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to SE mod_km=',mod_km
        elseif(mod_km<=11)then                                             !<-- S
          i_sw_parent_new=i_sw_parent_current
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to S mod_km=',mod_km
        elseif(mod_km<=13)then                                             !<-- SW
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to SW mod_km=',mod_km
        elseif(mod_km<=15)then                                             !<-- W
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current
      write(0,*)' artificial to W mod_km=',mod_km
        endif
      endif
      write(0,*)' exit artificial_move ntimestep=',ntimestep,' mod=',mod(ntimestep+3,51) &
               ,' i_want_to_move=',i_want_to_move
!
!-----------------------------------------------------------------------
      end subroutine artificial_move3
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine artificial_move4(ntimestep                              &
                                ,kount_moves                            &
                                ,i_want_to_move                         &
                                ,i_sw_parent_current                    &
                                ,j_sw_parent_current                    &
                                ,i_sw_parent_new                        &
                                ,j_sw_parent_new )
!
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: i_sw_parent_current              &
                                      ,j_sw_parent_current              &
                                      ,kount_moves                      &
                                      ,ntimestep
!
      integer(kind=kint),intent(out) :: i_sw_parent_new                 &
                                       ,j_sw_parent_new
!
      logical(kind=klog),intent(out) :: i_want_to_move
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      i_want_to_move=.false.
!
      write(0,*)' enter artificial ntimestep=',ntimestep,' kount_moves=',kount_moves
      if(ntimestep>3.and.mod(ntimestep,14)<=2)then
        i_want_to_move=.true.   
      write(0,*)' artificial set i_want_to_move=',i_want_to_move
        if(mod(kount_moves,45)<=5)then                                     !<-- E
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current
      write(0,*)' artificial to NW'
        elseif(mod(kount_moves,45)<=10)then                                !<-- NE
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to N'
        elseif(mod(kount_moves,45)<=15)then                                !<-- N
          i_sw_parent_new=i_sw_parent_current  
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to NE'
        elseif(mod(kount_moves,45)<=20)then                                !<-- NW
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current+1
      write(0,*)' artificial to E'
        elseif(mod(kount_moves,45)<=25)then                                !<-- W
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current  
      write(0,*)' artificial to SE'
        elseif(mod(kount_moves,45)<=30)then                                !<-- SW
          i_sw_parent_new=i_sw_parent_current-1
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to S'
        elseif(mod(kount_moves,45)<=35)then                                !<-- S
          i_sw_parent_new=i_sw_parent_current  
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to SW'
        elseif(mod(kount_moves,45)<=40)then                                !<-- SW
          i_sw_parent_new=i_sw_parent_current+1
          j_sw_parent_new=j_sw_parent_current-1
      write(0,*)' artificial to W'
        endif
      endif
      write(0,*)' exit artificial_move ntimestep=',ntimestep,' mod=',mod(ntimestep,14) &
               ,' i_want_to_move=',i_want_to_move
!
!-----------------------------------------------------------------------
      end subroutine artificial_move4
!-----------------------------------------------------------------------
!      
      END MODULE MODULE_PARENT_CHILD_CPL_COMP
!
!-----------------------------------------------------------------------
