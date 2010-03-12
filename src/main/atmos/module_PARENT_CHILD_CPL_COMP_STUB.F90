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
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_NEST_REG=ESMF_SUCCESS

!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_CPL_REGISTER
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
                                           ,ID_PARENTS                  &  !     |   
                                           ,N_CONFIGURE                 &  !     |   
                                           ,MAX_DOMAINS                 &  !   Input 
!                                                                           -----------
                                           ,IMP_STATE_CPL_NEST          &  !   Output
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
                                                           ,ID_PARENTS              !<-- IDs of parents of nested domains
!
      INTEGER(kind=KINT),DIMENSION(MAX_DOMAINS),INTENT(IN) :: N_CONFIGURE  !<-- Association of domains with configure file IDs
!
      REAL(kind=KFPT),DIMENSION(1:NUM_DOMAINS),INTENT(IN) :: DT            !<-- Timesteps for all domains (ATM Components)
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
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_COUPLER_SETUP
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_PARENT_CHILD_CPL_COMP
!
!-----------------------------------------------------------------------
