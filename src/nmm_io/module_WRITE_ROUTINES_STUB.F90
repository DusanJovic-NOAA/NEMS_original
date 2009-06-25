!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_ROUTINES
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
!       06 Sep 2007:  T. Black - Created module.
!                                Moved the FIRST block from the old
!                                write component here for clarity.
!       15 Aug 2008:  J. Wang  - Revised for NEMS-IO
!       20 Aug 2008:  J. Wang  - Output start date first instead of
!                                forecast date.
!       16 Sep 2008:  J. Wang  - WRITE_NEMSIO_RUNHISTORY_OPEN only
!                                opens file and writes metadata.
!       30 Sep 2008:  E. Colon - Generalize counts for nemsio
!       14 Oct 2008:  R. Vasic - Added restart capability
!       05 Jan 2009:  J. Wang  - Added 10-m wind factor to NMMB
!                                runhistory and restart files.
!       06 Jan 2009:  T. Black - Replace Max # of words recv'd by
!                                Write tasks with actual # of words.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE MODULE_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE
!
!jw      USE MODULE_INCLUDE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: WRITE_ASYNC                                             &
               ,WRITE_INIT                                              
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_INIT(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
! 
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEP OF THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE             !<-- The ATM Internal State
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
!
      END SUBROUTINE WRITE_INIT
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_ASYNC(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM,MYPE &
                            ,CWRT)
!
!-----------------------------------------------------------------------
!***  WRITE OUT A HISTORY FILE USING THE ASYNCHRONOUS QUILTING.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE             !<-- The ATM Internal State
      TYPE(ESMF_Clock)        ,INTENT(INOUT) :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
      INTEGER,INTENT(IN) :: MYPE
      CHARACTER(ESMF_MAXSTR) :: CWRT                                      !<-- Restart/History label
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_ASYNC
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_WRITE_ROUTINES
!
!-----------------------------------------------------------------------
