!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_GRID_COMP_GFS
!
!-----------------------------------------------------------------------
!***  THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!***  DATA WAS PUT INTO THIS COMPONENT'S IMPORT STATE DESTINED FOR
!***  HISTORY OUTPUT.  THIS COMPONENT EXTRACTS THAT INFORMATION
!***  FROM THE IMPORT STATE WHOSE CONTENTS ARE SEEN ONLY BY THE
!***  FORECAST TASKS AND TRANSFERS 2D DATA TO GROUPS OF WRITE TASKS
!***  WHERE IT IS PARTIALLY REASSEMBLED.  THE WRITE TASKS THEN
!***  TRANSFER THEIR SUBSECTIONS TO THE LEAD WRITE TASK WHICH
!***  ASSEBMLES THE 2D DATA ONTO THE FULL DOMAIN AND WRITES OUT
!***  ALL SCALAR/1D/2D DATA TO A HISTORY FILE.
!-----------------------------------------------------------------------
!***
!***  HISTORY   
!***
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions in CPL_REGISTER  
!                                and CPL_INITIALIZE
!       14 Aug 2007:  T. Black - Revised CPL_RUN for general output
!                                selection and added documentation
!                                for users.
!       12 Sep 2007:  T. Black - Replaced the write component and the
!                                write gridded component with only
!                                a gridded component that contains
!                                quilting.
!          Mar 2008:  R. Vasic - Convert from ESMF 3.0.1 to 3.1.0
!       15 Aug 2008:  J. Wang  - Revised for addition of NEMS-IO
!       16 Sep 2008:  J. Wang  - Output array reverts from 3-D to 2-D
!       14 Oct 2008:  R. Vasic - Add restart capability
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_WRITE_INTERNAL_STATE_GFS
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: WRITE_REGISTER_GFS
!
      PUBLIC :: WRITE_SETUP_GFS,WRITE_DESTROY_GFS
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_REGISTER_GFS(WRT_COMP,RC_WRT)
! 
!-----------------------------------------------------------------------
!***  REGISTER THE WRITE COMPONENT'S 
!***  INITIALIZE, RUN, AND FINALIZE SUBROUTINE NAMES.
!-----------------------------------------------------------------------
!
!***  HISTORY   
!       xx Feb 2007:  W. Yang  - Originator
!       30 Jun 2007:  T. Black - Modified to share same traits as
!                                rest of code.
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: WRT_COMP                     ! The write component
      INTEGER,INTENT(OUT)               :: RC_WRT                       ! Final return code
!     
!----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_REGISTER_GFS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRT_INITIALIZE_GFS(WRT_COMP                            &
                               ,IMP_STATE_WRITE                         &
                               ,EXP_STATE_WRITE                         &
                               ,CLOCK                                   &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  INITIALIZE THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
!***  HISTORY   
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions
!       29 Jun 2007:  T. Black - Generalize output types; add comment
!                                descriptions.
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE_WRITE  
      TYPE(ESMF_GridComp),INTENT(INOUT) :: WRT_COMP
      TYPE(ESMF_State)   ,INTENT(INOUT) :: EXP_STATE_WRITE  
!
      TYPE(ESMF_Clock)   ,INTENT(IN)    :: CLOCK
!
      INTEGER,INTENT(OUT)               :: RC_INIT
!
!----------------------------------------------------------------------- 
!
      END SUBROUTINE WRT_INITIALIZE_GFS
!
!----------------------------------------------------------------------- 
!####################################################################### 
!----------------------------------------------------------------------- 
!
      SUBROUTINE WRT_RUN_GFS(WRT_COMP                                   &
                        ,IMP_STATE_WRITE                                &
                        ,EXP_STATE_WRITE                                &
                        ,CLOCK                                          &
                        ,RC_RUN)
!
!----------------------------------------------------------------------- 
!***  THE RUN STEP FOR THE WRITE GRIDDED COMPONENT.  
!***  MOVE DATA INTENDED FOR HISTORY OUTPUT FROM THE IMPORT STATE
!***  TO THE WRITE TASKS.
!----------------------------------------------------------------------- 
!-----------------------------------------------------------------------
!
!***  HISTORY   
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions
!       14 Aug 2007:  T. Black - Major revisions for generalized
!                                selectable history output with
!                                quilting.  Add descriptive comments.
!
!-----------------------------------------------------------------------
!
      USE ESMF_FieldGetMOD
!
      TYPE(ESMF_GridComp),INTENT(IN) :: WRT_COMP
      TYPE(ESMF_Clock)   ,INTENT(IN) :: CLOCK
! 
      TYPE(ESMF_State),INTENT(INOUT) :: IMP_STATE_WRITE  
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_WRITE                  !<-- The write component export state.
                                                                         !    Although it is loaded up only as output from
                                                                         !    this subroutine, its INTENT needs to be INOUT
                                                                         !    to function properly.
      INTEGER,INTENT(OUT)            :: RC_RUN 
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRT_RUN_GFS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRT_FINALIZE_GFS(WRT_COMP                                  &
                             ,IMP_STATE_WRITE                           &
                             ,EXP_STATE_WRITE                           &
                             ,CLOCK                                     &
                             ,RCFINAL)
!
!-----------------------------------------------------------------------
!***  FINALIZE THE WRITE GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
!***  HISTORY
!       xx Feb 2007:  W. Yang - Originator
!       13 Jun 2007:  T. Black - Name revisions
!
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT):: WRT_COMP
      TYPE(ESMF_State)  ,INTENT(INOUT) :: IMP_STATE_WRITE  
      TYPE(ESMF_State),  INTENT(INOUT) :: EXP_STATE_WRITE  
      TYPE(ESMF_Clock)  ,INTENT(INOUT) :: CLOCK
!
      INTEGER,INTENT(OUT)              :: RCFINAL
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRT_FINALIZE_GFS
!
!-----------------------------------------------------------------------
!#######################################################################
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
      TYPE(ESMF_GridComp),INTENT(inOUT)      :: WRT_COMPS(:)  !<-- The ATM gridded component
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_DYN
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_PHY
      TYPE(ESMF_State),INTENT(INOUT)         :: IMP_state_WRITE
      TYPE(ESMF_State),INTENT(INOUT)         :: EXP_state_WRITE
!
!
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_DESTROY_GFS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      END MODULE MODULE_WRITE_GRID_COMP_GFS
!
!-----------------------------------------------------------------------
