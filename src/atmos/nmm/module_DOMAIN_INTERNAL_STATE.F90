!---------------------------------------------------------------------------
!
      MODULE MODULE_DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!***  Define all quantities that lie within the DOMAIN component's
!***  internal state.
!---------------------------------------------------------------------------
!
      USE ESMF_Mod
!
!---------------------------------------------------------------------------
!
      IMPLICIT NONE
!
!---------------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: DOMAIN_INTERNAL_STATE                                       &
               ,WRAP_DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
      TYPE DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
        TYPE(ESMF_GridComp),ALLOCATABLE,DIMENSION(:) :: DOMAIN_CHILD_COMP      !<-- DOMAIN components of child domains
!
        TYPE(ESMF_GridComp) :: DYN_GRID_COMP                                   !<-- The Dynamics gridded component
        TYPE(ESMF_GridComp) :: PHY_GRID_COMP                                   !<-- The Physics gridded component
        TYPE(ESMF_CplComp)  :: COUPLER_DYN_PHY_COMP                            !<-- The Dynamics-Physics coupler component
!
        TYPE(ESMF_State) :: IMP_STATE_DYN                                      !<-- The import state of the Dynamics component
        TYPE(ESMF_State) :: IMP_STATE_PHY                                      !<-- The import state of the Physics component
        TYPE(ESMF_State) :: IMP_STATE_WRITE                                    !<-- The import state of the write components
!
        TYPE(ESMF_State) :: EXP_STATE_DYN                                      !<-- The export state of the Dynamics component
        TYPE(ESMF_State) :: EXP_STATE_PHY                                      !<-- The export state of the Physics component
        TYPE(ESMF_State) :: EXP_STATE_WRITE                                    !<-- The export state of the write components
!
        INTEGER :: LEAD_TASK_DOMAIN                                            !<-- The first task on a given domain
        INTEGER :: NUM_PES_FCST                                                !<-- The number of forecast tasks
!
!---------------------------------------------------------------------------
!***  THE FOLLOWING ARE SPECIFIC TO ASYNCHRONOUS QUILTING/WRITING
!---------------------------------------------------------------------------
!
        LOGICAL :: QUILTING                                                    !<-- Is the user selecting asynchronous quilting/writing?
!
        TYPE(ESMF_GridComp),DIMENSION(:),POINTER :: WRITE_COMPS                !<-- The array of Write gridded components
!
        INTEGER :: WRITE_GROUPS                                                !<-- The number of write groups
        INTEGER :: WRITE_GROUP_READY_TO_GO                                     !<-- The active group of write tasks
        INTEGER :: WRITE_TASKS_PER_GROUP                                       !<-- The number of write tasks in each write group
!
        INTEGER,DIMENSION(:)  ,POINTER :: PETLIST_FCST                         !<-- Task ID list of fcst tasks (for Dyn and Phy components)
        INTEGER,DIMENSION(:,:),POINTER :: PETLIST_WRITE                        !<-- Task ID list of fcst tasks w/ write tasks by group
!
        INTEGER,DIMENSION(:),POINTER :: LOCAL_ISTART,LOCAL_IEND             &  !<-- The local I,J limits of the forecast tasks
                                       ,LOCAL_JSTART,LOCAL_JEND
!
!---------------------------------------------------------------------------
!
      END TYPE DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
!---------------------------------------------------------------------------
!***  This state is supported by C pointers but not by F90 pointers
!***  therefore we use this "WRAP".
!---------------------------------------------------------------------------
!
      TYPE WRAP_DOMAIN_INTERNAL_STATE
        TYPE(DOMAIN_INTERNAL_STATE),POINTER :: DOMAIN_INT_STATE
      END TYPE WRAP_DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
!
      END MODULE MODULE_DOMAIN_INTERNAL_STATE
!
!---------------------------------------------------------------------------
