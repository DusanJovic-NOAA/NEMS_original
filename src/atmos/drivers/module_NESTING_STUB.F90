!-----------------------------------------------------------------------
!
      MODULE MODULE_NESTING
!
!-----------------------------------------------------------------------
!
!***  THIS MODULE CONTAINS ROUTINES THAT PERFORM VARIOUS INTERACTIONS
!***  BETWEEN PARENT DOMAINS AND THEIR CHILDREN.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!
!   2008-02-07  Black - PARENT_TO_CHILD_FILL
!   2008-03-05  Black - PARENT_CHILD_SPLIT
!   2008-03-25  Black - PARENT_TO_CHILD_INIT_NMM
!   2008-04-22  Black - Replace PARENT_CHILD_SPLIT with _COMMS
!   2008-06-18  Black - PARENT_TO_CHILD_COMPUTE
!   2008-06-18  Black - PREPARE_PARENT_TO_CHILD_INTERP
!   2008-08-14  Black - Added BOUNDARY_DATA_STATE_TO_STATE
!   2009-03-12  Black - Added Z0BASE and STDH now needed for NPS.
!
!-----------------------------------------------------------------------
!
! USAGE: 
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
!      USE MODULE_INCLUDE
!
!      USE MODULE_ERR_MSG    ,ONLY: ERR_MSG,MESSAGE_CHECK
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: BOUNDARY_DATA_STATE_TO_STATE                            &
               ,PARENT_DATA_TO_ATM                                      &
               ,PARENT_CHILD_COMMS                                      &
               ,PARENT_TO_CHILD_INIT_NMM 
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: NEAREST=0                                    &  !<-- Flag for nearest neighbor interpolation (parent to child)
                          ,BILINEAR=1                                      !<-- Flag for bilinear interpolation (parent to child)
!
      INTEGER,SAVE :: MAX_NUM_DOMAINS=99
!
      INTEGER,SAVE :: LM
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_TO_CHILD_INIT_NMM(MYPE                          &
                                         ,CF                            &
                                         ,MY_DOMAIN_ID                  &
                                         ,THIS_CHILD_ID                 &
                                         ,DYN_GRID_COMP                 &
                                         ,PHY_GRID_COMP                 & 
                                         ,COMM_MY_DOMAIN )
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN)                         :: MYPE                &  !<-- My MPI task rank
                                                   ,MY_DOMAIN_ID        &  !<-- IDs of each domain
                                                   ,THIS_CHILD_ID       &  !<-- ID of the current child
                                                   ,COMM_MY_DOMAIN         !<-- MPI intracommunicator for individual domains
!
      TYPE(ESMF_Config),DIMENSION(*),INTENT(INOUT) :: CF                   !<-- The config objects (one per domain)
!
      TYPE(ESMF_GridComp),INTENT(INOUT)          :: DYN_GRID_COMP       &  !<-- The parent's Dynamics gridded component
                                                   ,PHY_GRID_COMP          !<-- The parent's Physics gridded component
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_TO_CHILD_INIT_NMM
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_DATA_TO_ATM(EXP_STATE_DYN                       &
                                   ,EXP_STATE_PHY                       &
                                   ,EXP_STATE_ATM )
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State),INTENT(IN)  :: EXP_STATE_DYN                     &   !<-- Dynamics export state
                                     ,EXP_STATE_PHY                         !<-- Physics export state
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_ATM                       !<-- ATM export state into which fcst Arrays are transferred
!
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_DATA_TO_ATM
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE BOUNDARY_DATA_STATE_TO_STATE(CLOCK                     &
                                             ,RATIO                     &
                                             ,STATE_IN                  &
                                             ,STATE_OUT )
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Clock),INTENT(IN),OPTIONAL :: CLOCK                        !<-- ESMF Clock
!
      INTEGER,INTENT(IN),OPTIONAL          :: RATIO                        !<-- # of child timesteps per parent timestep          
!
      TYPE(ESMF_State),INTENT(INOUT)       :: STATE_IN                     !<-- Input ESMF State
!
      TYPE(ESMF_State),INTENT(INOUT)       :: STATE_OUT                    !<-- Output ESMF State
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE BOUNDARY_DATA_STATE_TO_STATE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_CHILD_COMMS(MYPE                                &
                                   ,NUM_DOMAINS                         &
                                   ,CF                                  &
                                   ,COMM_FULL_DOMAIN                    &
                                   ,MY_DOMAIN_ID                        &
                                   ,ID_DOMAINS                          &
                                   ,ID_PARENTS                          &
                                   ,NUM_CHILDREN                        &
                                   ,ID_CHILDREN                         &
                                   ,COMM_MY_DOMAIN                      &
                                   ,COMM_TO_MY_PARENT                   &
                                   ,COMM_TO_MY_CHILDREN                 &
                                   ,FTASKS_DOMAIN                       &
                                   ,NTASKS_DOMAIN                       &
                                   ,PETLIST_ATM                         &
                                         )
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT (IN) :: COMM_FULL_DOMAIN                           &   !<-- MPI intracommunicator for ALL tasks
                            ,MYPE                                       &   !<-- My task ID
                            ,NUM_DOMAINS                                    !<-- Total number of domains
!
      TYPE(ESMF_Config),DIMENSION(*),INTENT(INOUT) :: CF                    !<-- The config objects (one per domain)
!
      INTEGER,DIMENSION(:,:),POINTER,INTENT(OUT) :: ID_CHILDREN         &   !<-- Domain IDs of all domains' children
                                                   ,PETLIST_ATM             !<-- List of task IDs on each domain
!
      INTEGER                       ,INTENT(OUT) :: COMM_MY_DOMAIN      &   !<-- Communicators for each individual domain
                                                   ,COMM_TO_MY_PARENT   &   !<-- Communicators to each domain's parent
                                                   ,MY_DOMAIN_ID            !<-- ID of the domain on which this task resides
!
      INTEGER,DIMENSION(:)  ,POINTER,INTENT(OUT) :: COMM_TO_MY_CHILDREN &   !<-- Communicators to each domain's children
                                                   ,ID_DOMAINS          &   !<-- Array of the domain IDs
                                                   ,ID_PARENTS          &   !<-- Array of the domains' parent IDs
                                                   ,FTASKS_DOMAIN       &   !<-- # of forecast tasks on each domain excluding descendents
                                                   ,NTASKS_DOMAIN       &   !<-- # of tasks on each domain excluding descendents
                                                   ,NUM_CHILDREN            !<-- # of children on each domain
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_COMMS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_NESTING
!
!-----------------------------------------------------------------------
