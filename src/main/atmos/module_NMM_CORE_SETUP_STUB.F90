      MODULE MODULE_NMM_CORE_SETUP
!-----------------------------------------------------------------------
!
!***  THIS MODULE HOLDS THE DYNAMICS REGISTER, INIT, RUN, AND FINALIZE
!***  ROUTINES.  THEY ARE CALLED FROM THE ATM GRIDDED COMPONENT
!***  (ATM INITIALIZE CALLS DYNAMICS INITIALIZE, ETC.)
!***  IN MODULE_MAIN_GRID_COMP.F.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
!
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
!
      PUBLIC :: NMM_SETUP   
!
!
!-----------------------------------------------------------------------
      INCLUDE '../../../inc/kind.inc'
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      SUBROUTINE NMM_SETUP(MYPE                                         &
                          ,COMM_MY_DOMAIN                               &
                          ,CF                                           &
                          ,ATM_GRID_COMP                                &
                          ,ATM_INT_STATE                                &
                          ,GRID_NMM_ATM)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE CONTAINS NMM-SPECIFIC CODE FOR THE ATM COMPONENT:
!***    (1) SETTING UP DISTRIBUTED MEMORY PARALLELISM IN THE NMM;
!***    (2) CREATING THE ESMF Grid FOR THE ATM COMPONENT;
!***    (3) SHARING LOCAL SUBDOMAIN INDEX LIMITS AMONG TASKS.
!-----------------------------------------------------------------------
!
      USE MODULE_ATM_INTERNAL_STATE
!
!
!-----------------------------------------------------------------------
!***  SET UP THE ESMF GRID FOR THE NMM AND ESTABLISH THE
!***  DISTRIBUTED MEMORY PARALLELISM.
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: MYPE                                        &
                           ,COMM_MY_DOMAIN
!
      TYPE(ESMF_Config)       ,INTENT(IN)    :: CF                         !<-- The configure object
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: ATM_GRID_COMP              !<-- The ATM gridded component
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE              !<-- The ATM Internal State
      TYPE(ESMF_Grid)         ,INTENT(OUT)   :: GRID_NMM_ATM               !<-- The ESMF GRID for the NMM integration grid
!
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_SETUP
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_NMM_CORE_SETUP
!
!-----------------------------------------------------------------------
