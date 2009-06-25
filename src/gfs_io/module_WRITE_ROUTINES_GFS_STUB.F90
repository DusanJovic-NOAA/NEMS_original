!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_ROUTINES_GFS
!
!-----------------------------------------------------------------------
!***  THIS MODULE CONTAINS ROUTINES NEEDED BY THE RUN STEP OF THE
!***  WRITE GRIDDED COMPONENT IN WHICH HISTORY OUTPUT DATA FROM
!***  THE FORECAST TASKS ARE ASSEMBLED AND WRITTEN TO HISTORY FILES
!***  BY THE WRITE TASKS.
!-----------------------------------------------------------------------
!***
!***  HISTORY   
!***
!       07 May 2009:  J. Wang - adopt write _routine from NMMB_io
!                                modified for GFS
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE MODULE_WRITE_INTERNAL_STATE_GFS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: WRITE_ASYNC_GFS,WRITE_INIT_GFS 
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_INIT_GFS(ATM_GRID_COMP,wrt_comps,imp_state_write, &
     &  exp_state_write,CLOCK_ATM,WRITE_GROUP_READY_TO_GO)
! 
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEP OF THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ESMF_GRIDComp),INTENT(INOUT)      :: wrt_comps(:)             !<-- The ATM Internal State
      TYPE(ESMF_STATE),INTENT(INOUT)         :: imp_state_write             !<-- The ATM Internal State
      TYPE(ESMF_STATE),INTENT(INOUT)         :: exp_state_write             !<-- The ATM Internal State
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
      INTEGER, INTENT(INOUT)                 :: WRITE_GROUP_READY_TO_GO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_INIT_GFS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_ASYNC_GFS(WRT_COMPS,EXP_STATE,IMP_STATE_WRITE     &
     &                      ,EXP_STATE_WRITE,CLOCK_ATM,MYPE              &
     &                      ,WRITE_GROUP_READY_TO_GO)
!
!-----------------------------------------------------------------------
!***  WRITE OUT A HISTORY FILE USING THE ASYNCHRONOUS QUILTING.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: WRT_COMPS(:)               !<-- The write_comp
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: EXP_STATE                  !<-- The export state dyn
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: IMP_STATE_WRITE            !<-- The import state of write
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: EXP_STATE_WRITE            !<-- The export state of write
      TYPE(ESMF_Clock)        ,INTENT(INOUT) :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
      INTEGER,INTENT(IN) :: MYPE
      INTEGER,INTENT(INOUT) :: WRITE_GROUP_READY_TO_GO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_ASYNC_GFS
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_WRITE_ROUTINES_GFS
!
!-----------------------------------------------------------------------
