!-----------------------------------------------------------------------
!
      MODULE MODULE_GFS_WRITE
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
!       03 Sep 2009:  W. Yang - Ensemble GEFS.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: WRITE_ASYNC_GFS, WRITE_INIT_GFS,              &
                WRITE_SETUP_GFS, WRITE_DESTROY_GFS 
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_INIT_GFS(ATM_GRID_COMP,wrt_comps,imp_state_write, &
        exp_state_write,CLOCK_ATM,WRITE_GROUP_READY_TO_GO)

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
      END SUBROUTINE WRITE_INIT_GFS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_ASYNC_GFS(WRT_COMPS,EXP_STATE,IMP_STATE_WRITE     &
                            ,EXP_STATE_WRITE,CLOCK_ATM,MYPE              &
                            ,WRITE_GROUP_READY_TO_GO)
!
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: WRT_COMPS(:)               !<-- The write_comp
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: EXP_STATE                  !<-- The export state dyn
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: IMP_STATE_WRITE            !<-- The import state of write
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: EXP_STATE_WRITE            !<-- The export state of write
      TYPE(ESMF_Clock)        ,INTENT(INOUT) :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
      INTEGER,INTENT(IN) :: MYPE
      INTEGER,INTENT(INOUT) :: WRITE_GROUP_READY_TO_GO
!
      END SUBROUTINE WRITE_ASYNC_GFS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_SETUP_GFS(ATM_GRID_COMP,WRT_COMPS           &
        ,exp_state_dyn,exp_state_phy                  &
        ,imp_state_write,exp_state_write)
! 
!-----------------------------------------------------------------------
!***  SET UP THE WRITE COMPONENTS WITH THE FORECAST TASKS AND
!***  THE GROUPS OF WRITE TASKS NEEDED FOR QUILTING THE OUTPUT
!***  AND WRITING IT TO HISTORY FILES.
!-----------------------------------------------------------------------
!
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ESMF_GridComp),INTENT(inOUT)      :: WRT_COMPS(:)              !<-- The ATM gridded component
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_DYN
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_PHY
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_state_WRITE
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_WRITE
!
      END SUBROUTINE WRITE_SETUP_GFS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_DESTROY_GFS(ATM_GRID_COMP,WRT_COMPS,             &
        IMP_STATE_WRITE,EXP_STATE_WRITE,CLOCK_ATM)
! 
!-----------------------------------------------------------------------
!***  DESTROY ALL OBJECTS RELATED TO THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ESMF_GridComp),DIMENSION(:),INTENT(INOUT)      ::WRT_COMPS
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_STATE_WRITE
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_STATE_WRITE
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
!
!
      END SUBROUTINE WRITE_DESTROY_GFS
!
!-----------------------------------------------------------------------
      END MODULE MODULE_GFS_WRITE
!
!-----------------------------------------------------------------------
