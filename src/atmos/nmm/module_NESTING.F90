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
!   2009-10-12  Black - Fix for generalized of parent-child space ratios.
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
      USE MODULE_LS_NOAHLSM ,ONLY: NUM_SOIL_LAYERS
!
      USE MODULE_ERR_MSG    ,ONLY: ERR_MSG,MESSAGE_CHECK
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
               ,PARENT_CHILD_COMMS                                      &
               ,PARENT_DATA_TO_DOMAIN                                   &
               ,PARENT_TO_CHILD_INIT_NMM                                &
               ,SET_NEST_GRIDS
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER :: NEAREST=0                         &  !<-- Flag for nearest neighbor interpolation (parent to child)
                                     ,BILINEAR=1                           !<-- Flag for bilinear interpolation (parent to child)
!
      INTEGER(kind=KINT),SAVE :: MAX_NUM_DOMAINS=99
!
      INTEGER(kind=KINT),SAVE :: LM
!
!-----------------------------------------------------------------------
!
      CONTAINS
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
                                   ,RANK_TO_DOMAIN_ID                   &
                                         )
!
!-----------------------------------------------------------------------
!***  Create MPI intercommunicators between the tasks of a parent domain
!***  and those of all its 1st generation nests (children).
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument Variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: COMM_FULL_DOMAIN                 &   !<-- MPI intracommunicator for ALL tasks
                                      ,MYPE                             &   !<-- My task ID
                                      ,NUM_DOMAINS                          !<-- Total number of domains
!
      INTEGER(kind=KINT),DIMENSION(*),INTENT(IN) :: RANK_TO_DOMAIN_ID       !<-- Domain ID for each configure file
!
      TYPE(ESMF_Config),DIMENSION(*),INTENT(INOUT) :: CF                    !<-- The config objects (one per domain)
!
      INTEGER(kind=KINT),DIMENSION(:,:),POINTER,INTENT(OUT) :: ID_CHILDREN &   !<-- Domain IDs of all domains' children
                                                              ,PETLIST_ATM     !<-- List of task IDs on each domain
!
      INTEGER(kind=KINT),INTENT(OUT) :: COMM_MY_DOMAIN                  &   !<-- Communicators for each individual domain
                                       ,COMM_TO_MY_PARENT               &   !<-- Communicators to each domain's parent
                                       ,MY_DOMAIN_ID                        !<-- ID of the domain on which this task resides
!
      INTEGER(kind=KINT),DIMENSION(:),POINTER,INTENT(OUT) :: COMM_TO_MY_CHILDREN &   !<-- Communicators to each domain's children
                                                            ,ID_DOMAINS          &   !<-- Array of the domain IDs
                                                            ,ID_PARENTS          &   !<-- Array of the domains' parent IDs
                                                            ,FTASKS_DOMAIN       &   !<-- # of forecast tasks on each domain
                                                            ,NTASKS_DOMAIN       &   !<-- # of tasks on each domain excluding descendents
                                                            ,NUM_CHILDREN            !<-- # of children on each domain
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER(kind=KINT) :: I,ID_DOM,IERR,ISTAT,N,N1,N2,NN,RC
!
      INTEGER(kind=KINT) :: COMM_ALL                                    &
                           ,COMM_DUP                                    &
                           ,COMM_INTRA                                  &
                           ,COMM_PARENT_CHILD                           &
                           ,GROUP_ALL                                   &
                           ,GROUP_NEW                                   &
                           ,I_COLOR                                     &
                           ,ID_CHILD                                    &
                           ,ID_DOMX                                     &
                           ,INPES                                       &
                           ,JNPES                                       &
                           ,KOUNT_EXCLUDE                               &
                           ,KOUNT_TASKS                                 &
                           ,LEAD_REMOTE                                 &
                           ,N_CHILDREN                                  &
                           ,NTASKS_PARENT                               &
                           ,NUM_EXCLUDE                                 &
                           ,RC_COMMS                                    &
                           ,TASK_EXCLUDE                                &
                           ,TASK_X                                      &
                           ,TOTAL_TASKS                                 &
                           ,WRITE_GROUPS                                &
                           ,WRITE_TASKS_PER_GROUP 
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: EXCLUDE_TASKS      &  !<-- Tasks excluded from each parent-child couplet
                                                    ,LAST_TASK          &  !<-- ID of last task on each domain
                                                    ,LEAD_TASK          &  !<-- ID od first task on each domain
                                                    ,PARENT_OF_DOMAIN      !<-- IDs of parents of each domain
!
      CHARACTER(2)      :: NUM_DOMAIN
      CHARACTER(6),SAVE :: FMT='(I2.2)'
!
      LOGICAL(kind=KLOG) :: I_AM_PARENT                                 &
                           ,INCLUDE_MYPE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_COMMS=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!
      ALLOCATE(ID_DOMAINS      (1:NUM_DOMAINS))       
      ALLOCATE(ID_PARENTS      (1:NUM_DOMAINS))       
      ALLOCATE(LEAD_TASK       (1:NUM_DOMAINS))       
      ALLOCATE(LAST_TASK       (1:NUM_DOMAINS))       
      ALLOCATE(FTASKS_DOMAIN   (1:NUM_DOMAINS))       
      ALLOCATE(NTASKS_DOMAIN   (1:NUM_DOMAINS))       
      ALLOCATE(NUM_CHILDREN    (1:NUM_DOMAINS))                          
!
      TOTAL_TASKS=0
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Incoming tasks extract relevant information from all config files.
!
!***  This loop is general thus the domain IDs do not need to correspond
!***  to the number in the configure file name.  The user may assign
!***  IDs monotonic starting with 1 to the domains in any order
!***  desired except that the uppermost parent must have an ID of 1.
!***  However the rank/element of each domain in the CF array is equal
!***  to the given domain's ID.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      read_configs: DO N=1,NUM_DOMAINS                                     !<-- Loop through all configure objects
!
!-----------------------------------------------------------------------
!***  Save the domain IDs.
!***  These are simply integers each domain will use to keep track
!***  of itself with respect to others.
!-----------------------------------------------------------------------
!
        ID_DOMAINS(N)=RANK_TO_DOMAIN_ID(N)
        ID_DOM=ID_DOMAINS(N)
!
!-----------------------------------------------------------------------
!***  Who is the parent of each domain?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract ID of Parent of this Domain"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =ID_PARENTS(ID_DOM)          &  !<-- The ID of the parent of this domain
                                    ,label ='my_parent_id:'             &  !<-- Take values from this config label 
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!***  How many children does each domain have?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract # of Children of this Domain"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =NUM_CHILDREN(ID_DOM)        &  !<-- # of children on this domain
                                    ,label ='n_children:'               &  !<-- Take value from this config label 
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  How many forecast/write tasks will be active on each domain?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent_Child_Comms: Extract INPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =INPES                       &  !<-- The domain's fcst tasks in I
                                    ,label ='inpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent_Child_Comms: Extract JNPES From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =JNPES                       &  !<-- The domain's fcst tasks in J
                                    ,label ='jnpes:'                    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent_Child_Comms: Extract Write_Groups From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =WRITE_GROUPS                &  !<-- The number of Write groups on this domain
                                    ,label ='write_groups:'             &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Parent_Child_Comms: Extract Write_Task_Per_Group From Config File"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOM)                  &  !<-- The config object
                                    ,value =WRITE_TASKS_PER_GROUP       &  !<-- The number of tasks per Write group 
                                    ,label ='write_tasks_per_group:'    &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_COMMS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        FTASKS_DOMAIN(ID_DOM)=INPES*JNPES                                  !<-- # of forecast tasks on each domain
!
        NTASKS_DOMAIN(ID_DOM)=INPES*JNPES+                              &  !<-- Total # of tasks on each domain
                              WRITE_GROUPS*WRITE_TASKS_PER_GROUP
!
        TOTAL_TASKS=TOTAL_TASKS+NTASKS_DOMAIN(ID_DOM)                      !<-- Sum the total task count among all domains
!
!-----------------------------------------------------------------------
!
      ENDDO read_configs
!
      ALLOCATE(PETLIST_ATM(1:TOTAL_TASKS,1:NUM_DOMAINS))
!
!-----------------------------------------------------------------------
!***  Assign tasks to all domains in the order of their configure files.
!-----------------------------------------------------------------------
!
      DO N=1,NUM_DOMAINS
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
!-----------------------------------------------------------------------
!***  We must assign a "COLOR" to each task.
!***  Simply set the color equal to the domain ID.
!***  At the same time we are determining the global IDs
!***  of the lead and last tasks on each domain which
!***  in turn lets us fill the PETLIST for each domain.
!-----------------------------------------------------------------------
!
        IF(N==1)THEN
          LEAD_TASK(N)=0                                                   !<-- Task 0 is first in line
        ELSE
          LEAD_TASK(ID_DOM)=LAST_TASK(ID_DOMAINS(N-1))+1                   !<-- Lead task on domain follows last task on previous domain
        ENDIF
!
        LAST_TASK(ID_DOM)=LEAD_TASK(ID_DOM)+NTASKS_DOMAIN(ID_DOM)-1        !<-- The last task on each domain
!
        IF(MYPE>=LEAD_TASK(ID_DOM).AND.MYPE<=LAST_TASK(ID_DOM))THEN        !<-- Associate tasks with each domain              
          I_COLOR=ID_DOM                                                   !<-- Set color to domain ID
          MY_DOMAIN_ID=ID_DOM                                              !<-- Tell this task the ID of the domain it is on 
        ENDIF
!
        KOUNT_TASKS=0
        DO N2=LEAD_TASK(ID_DOM),LAST_TASK(ID_DOM)     
          KOUNT_TASKS=KOUNT_TASKS+1
          PETLIST_ATM(KOUNT_TASKS,ID_DOM)=N2
        ENDDO
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  Now split the input intracommunicator shared by all tasks
!***  into local intracommunicators called COMM_MY_DOMAIN on 
!***  each domain.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_DUP  (COMM_FULL_DOMAIN,COMM_DUP,IERR)                  !<-- Duplicate incoming intracommunicator to play it safe
      CALL MPI_COMM_SPLIT(COMM_DUP,I_COLOR,MYPE,COMM_MY_DOMAIN,IERR)       !<-- COMM_FULL_DOMAIN has been split (but still exists)
!
!-----------------------------------------------------------------------
!***  All tasks know the task counts and IDs of each domain as well as
!***  the parents of each domain.
!
!***  Loop through all domains in order to associate all parents
!***  with their children through intercommunicators.
!-----------------------------------------------------------------------
!
      ALLOCATE(ID_CHILDREN(1:NUM_DOMAINS,1:NUM_DOMAINS))                   !<-- Array to hold all domains' children's IDs
!
      DO N1=1,NUM_DOMAINS
      DO N2=1,NUM_DOMAINS
        ID_CHILDREN(N1,N2)=0                                               !<-- All genuine Domain IDs are >0
      ENDDO
      ENDDO
!
      COMM_TO_MY_PARENT=-999                                               !<-- Initialize communicator to parent to nonsense
!
!-----------------------------------------------------------------------
!
      main_loop: DO N=1,NUM_DOMAINS
!
        ID_DOM=RANK_TO_DOMAIN_ID(N)
!
!-----------------------------------------------------------------------
!
        N_CHILDREN=NUM_CHILDREN(ID_DOM)                                    !<-- The # of children on this domain
!
        IF(N_CHILDREN==0)CYCLE main_loop                                   !<-- If this domain has no children, move on
!
        NTASKS_PARENT=NTASKS_DOMAIN(ID_DOM)                                !<-- Total # of fcst and write tasks on this domain
!
!-----------------------------------------------------------------------
!***  All domain IDs will be searched in the configure files to find 
!***  matches between the current domain's ID and the parent IDs of 
!***  the other domains.
!***  Matches will identify parent-child couplets.
!-----------------------------------------------------------------------
!
        NN=0                                                               !<-- Index of children of the parent domain
!
        DO N2=1,NUM_DOMAINS                                                !<-- Search for children who have this parent
          ID_CHILD=ID_DOMAINS(N2)
!
          IF(ID_PARENTS(ID_CHILD)==ID_DOM)THEN                             !<-- If yes then we found a nest that is this domain's child
            NN=NN+1
            ID_CHILDREN(NN,ID_DOM)=ID_CHILD                                !<-- IDs of this parent's (domain N's) children
          ENDIF
!
          IF(NN==N_CHILDREN)EXIT                                           !<-- We have found all of this domain's children 
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Be sure some essentials are present before we proceed.
!-----------------------------------------------------------------------
!
        IF(NTASKS_PARENT<1)THEN
          WRITE(0,*)' ERROR: THIS PARENT DOMAIN HAS NO TASKS'
          RETURN
        ENDIF
!
!-----------------------------------------------------------------------
!***  Parent tasks create their intercommunicator array.
!-----------------------------------------------------------------------
!
        I_AM_PARENT=.FALSE.
!
        IF(MYPE>=LEAD_TASK(ID_DOM).AND.MYPE<=LAST_TASK(ID_DOM))THEN        !<-- Select parent tasks              
          ALLOCATE(COMM_TO_MY_CHILDREN(1:N_CHILDREN))                      !<-- Parent allocates parent-to-child intercommunicators
          I_AM_PARENT=.TRUE.
!
          DO N2=1,N_CHILDREN
            COMM_TO_MY_CHILDREN(N2)=-999                                   !<-- Initialize intercommunicators to nonsense
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  Create intercommunicators between the parent and its
!***  children.  For each child exclude all tasks that are not
!***  associated with that child or its parent.
!
!***  First set the local group ID of the remote task as seen by
!***  the parent and by its children.
!-----------------------------------------------------------------------
!
        IF(I_AM_PARENT)THEN
          LEAD_REMOTE=NTASKS_PARENT                                        !<-- The lead child task as seen by its parent is always
                                                                           !      relative to the parent/child groupsize (starting 
                                                                           !      with the parent) and NOT to the global task ranks.
        ELSE
          LEAD_REMOTE=0                                                    !<-- The lead parent task as seen from its children;
                                                                           !      it is relative to groupsize and thus is always 0.
        ENDIF
!
!-----------------------------------------------------------------------
!
        inter_comm: DO N2=1,N_CHILDREN
!
          INCLUDE_MYPE =.TRUE.
!
!-----------------------------------------------------------------------
!***  Only if there are more than two domains (a parent and a single
!***  child) do we need to exclude tasks outside of that union of tasks.
!-----------------------------------------------------------------------
!
          exclude: IF(NUM_DOMAINS>2)THEN
!
            ID_CHILD=ID_CHILDREN(N2,ID_DOM)
            NUM_EXCLUDE=TOTAL_TASKS-NTASKS_PARENT-NTASKS_DOMAIN(ID_CHILD)  !<-- Task count of non-parent non-child tasks
            ALLOCATE(EXCLUDE_TASKS(1:NUM_EXCLUDE))                         !<-- This will hold excluded tasks
!
            KOUNT_EXCLUDE=0
            TASK_X       =0                                                !<-- First task that might be excluded from parent/child couplet
!
!-----------------------------------------------------------------------
!
            DO NN=1,NUM_DOMAINS                                            !<-- Loop through all domains
              ID_DOMX=ID_DOMAINS(NN)
!
              IF(ID_DOMX/=ID_DOM.AND.ID_DOMX/=ID_CHILD)THEN                !<-- Select tasks of all domains who are not parent/child
                TASK_EXCLUDE=TASK_X
!
                DO I=1,NTASKS_DOMAIN(ID_DOMX)                              !<-- Loop through tasks of domain that is not child "N" or parent
                  KOUNT_EXCLUDE=KOUNT_EXCLUDE+1                            !<-- Add up all tasks not in parent/child couplet
                  EXCLUDE_TASKS(KOUNT_EXCLUDE)=TASK_EXCLUDE                !<-- This task is not on child "N" or its parent thus is excluded
                  IF(MYPE==TASK_EXCLUDE)INCLUDE_MYPE=.FALSE.               !<-- Excluded tasks set their Include flag to false
                  TASK_EXCLUDE=TASK_EXCLUDE+1
                ENDDO
!
              ENDIF
!
              TASK_X=TASK_X+NTASKS_DOMAIN(ID_DOMX)
!
            ENDDO
!
            IF(NUM_EXCLUDE/=KOUNT_EXCLUDE)THEN
              WRITE(0,*)' ERROR: Count of excluded tasks went wrong. '  &
                       ,' NUM_EXCLUDE=',NUM_EXCLUDE                     &
                       ,' KOUNT_EXCLUDE=',KOUNT_EXCLUDE
            ENDIF
!
          ENDIF exclude
!
!-----------------------------------------------------------------------
!***  We now have array EXCLUDE_TASKS holding the IDs OF ALL tasks
!***  that are not part of the union of child "N" and its parent.
!***  The number of these excluded tasks is KOUNT_EXCLUDE.
!
!***  Create a new MPI group with only the tasks of the parent
!***  and child "N" and then create an intracommunicator that
!***  contains only that union of tasks.
!-----------------------------------------------------------------------
!
          COMM_ALL=COMM_FULL_DOMAIN                                        !<-- Communicator for all tasks entering this routine
          CALL MPI_COMM_GROUP(COMM_ALL,GROUP_ALL,IERR)                     !<-- GROUP_ALL contains all tasks entering this routine
!
          IF(NUM_DOMAINS>2)THEN                                            !<-- Exclusion of tasks possible only if NUM_DOMAINS>2
            CALL MPI_GROUP_EXCL(GROUP_ALL,KOUNT_EXCLUDE,EXCLUDE_TASKS   &
                               ,GROUP_NEW,IERR)                            !<-- GROUP_NEW contains tasks on parent and child "N"
            CALL MPI_COMM_CREATE(COMM_ALL,GROUP_NEW,COMM_INTRA,IERR)       !<-- COMM_INTRA is intracommunicator for tasks in GROUP_NEW
            CALL MPI_GROUP_FREE(GROUP_NEW,IERR)                            !<-- Clear (set to null) GROUP_NEW
          ELSE
            CALL MPI_COMM_CREATE(COMM_ALL,GROUP_ALL,COMM_INTRA,IERR)       !<-- COMM_INTRA is intracommunicator for tasks in GROUP_NEW
          ENDIF
!
          CALL MPI_GROUP_FREE(GROUP_ALL,IERR)                              !<-- Clear (set to null) GROUP_ALL
!
!-----------------------------------------------------------------------
!***  The intracommunicator called COMM_INTRA excludes all tasks
!***  outside the union of the current parent and its child "N".
!***  Create an intercommunicator between this parent and its child "N"
!***  Through a process in which only those tasks can participate.
!***  Save the intercommunicator in an array holding all of the
!***  intercommunicators between the current parent and its children.
!-----------------------------------------------------------------------
!
          IF(INCLUDE_MYPE)THEN
            CALL MPI_INTERCOMM_CREATE(COMM_MY_DOMAIN                    &  !<-- Each task's local intracommunicator
                                     ,0                                 &  !<-- Rank of lead task in each local intracommunicator
                                     ,COMM_INTRA                        &  !<-- Intracommunicator between tasks of current domain
                                                                           !      and child "N"
                                     ,LEAD_REMOTE                       &  !<-- Rank of leader in the remote group in COMM_INTRA
                                     ,0                                 &  !<-- A tag
                                     ,COMM_PARENT_CHILD                 &  !<-- The new intercommunicator between tasks of current
                                                                           !      domain and those of its child "N"
                                     ,IERR )
!
            IF(I_AM_PARENT)THEN                                  
              COMM_TO_MY_CHILDREN(N2)=COMM_PARENT_CHILD                    !<-- Parent: The intercommunicator is stored as parent-to-child
            ELSE
              COMM_TO_MY_PARENT      =COMM_PARENT_CHILD                    !<-- Child: The intercommunicator is stored as child-to-parent
            ENDIF
!
          ENDIF
!
          CALL MPI_BARRIER(COMM_MY_DOMAIN,IERR)
!
          IF(NUM_DOMAINS>2)DEALLOCATE(EXCLUDE_TASKS)
!
!-----------------------------------------------------------------------
!
        ENDDO inter_comm
!
!-----------------------------------------------------------------------
!
      ENDDO main_loop
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(LEAD_TASK)
      DEALLOCATE(LAST_TASK)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_CHILD_COMMS
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
!***  Initialize data in a child's domain with data from its parent.
!***  Rows and columns on the child's grid lie parallel to rows and
!***  colums on the parent's grid.
!***  Only parent tasks are needed for this.
!-----------------------------------------------------------------------
!
      USE MODULE_DYNAMICS_INTERNAL_STATE,ONLY: DYNAMICS_INTERNAL_STATE  &
                                              ,WRAP_DYN_INT_STATE
!
      USE MODULE_PHYSICS_INTERNAL_STATE ,ONLY: PHYSICS_INTERNAL_STATE   &
                                              ,WRAP_PHY_INT_STATE
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
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,L,M,N,NCHILD,NFCST,RC,RC_CHILD,LNSH=1    ! Tom, change lnsh value later
      INTEGER :: INPES,JNPES
      INTEGER :: IDS,IDE,JDS,JDE
      INTEGER :: IMS,IME,JMS,JME
      INTEGER :: ITS,ITE,JTS,JTE
      INTEGER :: NUM_PES_PARENT
      INTEGER :: IM_CHILD,JM_CHILD
      INTEGER :: IM_PARENT,JM_PARENT
      INTEGER :: I_PARENT_START,J_PARENT_START
      INTEGER :: IHREND,NTSD
      INTEGER :: PARENT_CHILD_SPACE_RATIO
!
      INTEGER,DIMENSION(:),POINTER      :: LOCAL_ISTART                 &
                                          ,LOCAL_IEND                   &
                                          ,LOCAL_JSTART                 &
                                          ,LOCAL_JEND
!
      INTEGER,ALLOCATABLE,DIMENSION(:,:)   :: IDUMMY_2D
!
      REAL                              :: CHILD_PARENT_SPACE_RATIO     &
                                          ,COL_0                        &
                                          ,DLMD                         &
                                          ,DLMD_CHILD                   &
                                          ,DLMD_PARENT                  &
                                          ,DPHD                         &
                                          ,DPHD_CHILD                   &
                                          ,DPHD_PARENT                  &
                                          ,ROW_0                        &
                                          ,SBD                          &
                                          ,SBD_CHILD                    &
                                          ,SBD_PARENT                   &
                                          ,SW_LATD_CHILD                &
                                          ,SW_LOND_CHILD                &
                                          ,TLM0D_CHILD                  &
                                          ,TLM0D_PARENT                 &
                                          ,TPH0D_CHILD                  &
                                          ,TPH0D_PARENT                 &
                                          ,WBD                          &
                                          ,WBD_CHILD                    &
                                          ,WBD_PARENT
!
      REAL,ALLOCATABLE,DIMENSION(:)     :: DUMMY_SOIL
!
      REAL,ALLOCATABLE,DIMENSION(:,:)   :: SEA_MASK,SEA_ICE
!
      REAL,ALLOCATABLE,DIMENSION(:,:,:) :: DUMMY_2D_IN                  &
                                          ,DUMMY_2D_OUT                 &
                                          ,DUMMY_3D                     &
                                          ,DUMMY_3DS                    &
                                          ,PD_BILINEAR                  &
                                          ,PD_NEAREST                   &
                                          ,TEMPSOIL
!
      CHARACTER(2)  :: INT_TO_CHAR
      CHARACTER(6)  :: FMT
      CHARACTER(50) :: GLOBAL_FLAG                                      &
                      ,OUTFILE    
!
      LOGICAL :: GLOBAL,OPENED
      LOGICAL,ALLOCATABLE,DIMENSION(:,:) :: LOWER_TOPO
!
      TYPE(WRAP_DYN_INT_STATE)  :: WRAP_DYN
      TYPE(WRAP_PHY_INT_STATE)  :: WRAP_PHY
      TYPE(DYNAMICS_INTERNAL_STATE),POINTER :: DYN_INT_STATE
      TYPE(PHYSICS_INTERNAL_STATE) ,POINTER :: PHY_INT_STATE
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE PROVIDES DATA TO A CHILD DOMAIN WHEN NO NORMAL
!***  PRE-PROCESSED INPUT HAS BEEN PREPARED.  THIS IS DONE BY SIMPLE
!***  BILINEAR INTERPOLATION FROM THE PARENT DOMAIN'S DATA.
!
!***  THE RECORD ORDER OF THE PARENT'S INPUT IS DUPLICATED AND THOSE 
!***  RECORDS ARE WRITTEN TO A FILE SO THAT THE CHILD CAN READ IN AND 
!***  DISTRIBUTE THE DATA JUST AS IF A NORMAL PRE-PROCESSED FILE WERE 
!***  BEING USED.
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC      =ESMF_SUCCESS
      RC_CHILD=ESMF_SUCCESS
!
      NULLIFY(DYN_INT_STATE)
      NULLIFY(PHY_INT_STATE)
!
!-----------------------------------------------------------------------
!***  ONLY FORECAST TASKS ARE RELEVANT AND HAVE CORRECT DATA
!***  FOR THIS WORK.  FIND THE NUMBER OF FORECAST TASKS FROM
!***  CONFIGURE FILE (ALL FORECAST TASKS ALREADY HAVE THIS
!***  IN THEIR INTERNAL STATE BUT THE WRITE TASKS DO NOT).
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract INPES From Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
                                  ,value =INPES                         &  !<-- The variable filled (parent fcst tasks in I)
                                  ,label ='inpes:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract JNPES From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
                                  ,value =JNPES                         &  !<-- The variable filled (parent fcst tasks in J)
                                  ,label ='jnpes:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      NUM_PES_PARENT=INPES*JNPES                                           !<-- # of parent fcst tasks
!
      IF(MYPE>=NUM_PES_PARENT)RETURN                                       !<-- Parent's quilt/write tasks may leave
!
!-----------------------------------------------------------------------
!***  WE NEED THE SPATIAL RESOLUTION OF THE PARENT GRID SO EXTRACT
!***  ITS DIMENSIONS AND ITS BOUNDS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract IM From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
                                  ,value =IM_PARENT                     &  !<-- The variable filled (I dimension of parent grid)
                                  ,label ='im:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract JM From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
                                  ,value =JM_PARENT                     &  !<-- The variable filled (J dimension of parent grid)
                                  ,label ='jm:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract SBD From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!  CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
!!!!                              ,value =SBD_PARENT                    &  !<-- The variable filled (South boundary of parent grid)
!!!!                              ,label ='sbd:'                        &  !<-- Give this label's value to the previous variable
!!!!                              ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract WBD From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!  CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
!!!!                              ,value =WBD_PARENT                    &  !<-- The variable filled (West boundary of parent grid)
!!!!                              ,label ='wbd:'                        &  !<-- Give this label's value to the previous variable
!!!!                              ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract TPH0D From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!  CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
!!!!                              ,value =TPH0D_PARENT                  &  !<-- The variable filled (Central lat of parent grid)
!!!!                              ,label ='tph0d:'                      &  !<-- Give this label's value to the previous variable
!!!!                              ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract TLM0D From Config File"
!!!   CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!!!!  CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)              &  !<-- The parent's config object
!!!!                              ,value =TLM0D_PARENT                  &  !<-- The variable filled (Central lon of parent grid)
!!!!                              ,label ='tlm0d:'                      &  !<-- Give this label's value to the previous variable
!!!!                              ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Global Flag for Parent Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(MY_DOMAIN_ID)             &  !<-- The configure object of my parent
                                  ,value =GLOBAL_FLAG                  &  !<-- The variable filled 
                                  ,label ='global:'                    &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD) 
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE PARENT GRID RESOLUTION.
!-----------------------------------------------------------------------
!
      IF(TRIM(GLOBAL_FLAG)=='true')THEN                                    !<-- Parent is global 
        GLOBAL=.TRUE.
        IDE=IM_PARENT+2
        JDE=JM_PARENT+2
!!!!    DPHD_PARENT=-SBD_PARENT*2./REAL(JDE-3)
!!!!    DLMD_PARENT=-WBD_PARENT*2./REAL(IDE-3)
      ELSE                                                                 !<-- Parent is regional
        GLOBAL=.FALSE.
        IDE=IM_PARENT
        JDE=JM_PARENT
!!!!    DPHD_PARENT=-SBD_PARENT*2./REAL(JDE-1)
!!!!    DLMD_PARENT=-WBD_PARENT*2./REAL(IDE-1)
      ENDIF
!
      ROW_0=0.5*(JDE+1)
      COL_0=0.5*(IDE+1)
!
!-----------------------------------------------------------------------
!***  EXTRACT THE DYNAMICS AND PHYSICS INTERNAL STATES OF THE PARENT
!***  SO WE CAN USE THEIR DATA FOR THE NESTS.
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGetInternalState(DYN_GRID_COMP                  &
                                        ,WRAP_DYN                       &
                                        ,RC )
!
      CALL ESMF_GridCompGetInternalState(PHY_GRID_COMP                  &
                                        ,WRAP_PHY                       &
                                        ,RC )
!
!-----------------------------------------------------------------------
!
      DYN_INT_STATE=>wrap_dyn%INT_STATE
      PHY_INT_STATE=>wrap_phy%INT_STATE
!
      IMS=dyn_int_state%IMS                                                !<-- Horizontal memory limits on parent tasks
      IME=dyn_int_state%IME                                                !
      JMS=dyn_int_state%JMS                                                !
      JME=dyn_int_state%JME                                                !<--
!
      ITS=dyn_int_state%ITS                                                !<-- Horizontal integration limits on parent tasks
      ITE=dyn_int_state%ITE                                                !
      JTS=dyn_int_state%JTS                                                !
      JTE=dyn_int_state%JTE                                                !<--
!
      LM=dyn_int_state%LM                                                  !<-- Number of atmospheric layers
!
      LOCAL_ISTART=>dyn_int_state%LOCAL_ISTART                             !<-- Local integration limits for all parent tasks
      LOCAL_IEND  =>dyn_int_state%LOCAL_IEND                               !
      LOCAL_JSTART=>dyn_int_state%LOCAL_JSTART                             !
      LOCAL_JEND  =>dyn_int_state%LOCAL_JEND                               !<--
!
!-----------------------------------------------------------------------
!***  DPHD/DLMD AND SBD/WBD ARE USED ONLY FOR STAND-ALONE, INDEPENDENT
!***  ROTATED PARENT/NEST GRIDS (I.E., NOT GRID-ASSOCIATED NESTS).
!-----------------------------------------------------------------------
!
      DPHD_PARENT=dyn_int_state%DPHD
      DLMD_PARENT=dyn_int_state%DLMD
      SBD_PARENT=dyn_int_state%SBD
      WBD_PARENT=dyn_int_state%WBD
      TPH0D_PARENT=dyn_int_state%TPH0D
      TLM0D_PARENT=dyn_int_state%TLM0D
!     write(0,*)' _INIT_NMM dphd_parent=',dphd_parent,' dlmd_parent=',dlmd_parent
!     write(0,*)' _INIT_NMM sbd_parent=',sbd_parent,' wbd_parent=',wbd_parent &
!              ,' tph0d_parent=',tph0d_parent,' tlm0d_parent=',tlm0d_parent
!
!-----------------------------------------------------------------------
!***  EXTRACT RELEVANT INFORMATION FROM THIS CHILD'S CONFIGURE FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract I of SW Point on Parent"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =I_PARENT_START                &  !<-- The variable filled (parent I of child's SW corner)
                                  ,label ='i_parent_start:'             &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract J of SW Point on Parent"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =J_PARENT_START                &  !<-- The variable filled (parent J of child's SW corner)
                                  ,label ='j_parent_start:'             &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract Child/Parent Grid Ratio"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =PARENT_CHILD_SPACE_RATIO      &  !<-- The variable filled (child grid increment / parent's)
                                  ,label ='parent_child_space_ratio:'   &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CHILD_PARENT_SPACE_RATIO=1./REAL(PARENT_CHILD_SPACE_RATIO)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract Global IM of Child"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =IM_CHILD                      &  !<-- The variable filled (IM of child domain)
                                  ,label ='im:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Child_Init: Extract Global JM of Child"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(THIS_CHILD_ID)             &  !<-- The child's config object
                                  ,value =JM_CHILD                      &  !<-- The variable filled (JM of child domain)
                                  ,label ='jm:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CHILD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ONLY FOR FREE-STANDING NESTS:
!
!***  WHAT IS THE PARENT LAT/LON OF THE SW CORNER H POINT OF THE
!***  CHILD GRID?
!***  FIND THE RESOLUTION, BOUNDS, AND CENTER OF THE CHILD GRID.
!-----------------------------------------------------------------------
!
!     CALL CONVERT_IJ_TO_LATLON  (I_PARENT_START                        &
!                                ,J_PARENT_START                        &
!                                ,IM_PARENT                             &
!                                ,JM_PARENT                             &
!                                ,TPH0D_PARENT                          &
!                                ,TLM0D_PARENT                          &
!                                ,DPHD_PARENT                           &
!                                ,DLMD_PARENT                           &
!                                ,SW_LATD_CHILD                         &
!                                ,SW_LOND_CHILD )
!     write(0,*)' I_PARENT_START=',I_PARENT_START,' J_PARENT_START=',J_PARENT_START,' IM_PARENT=',IM_PARENT,' JM_PARENT=',JM_PARENT
!     write(0,*)' TPH0D_PARENT=',TPH0D_PARENT,' TLM0D_PARENT=',TLM0D_PARENT,' DPHD_PARENT=',DPHD_PARENT,' DLMD_PARENT=',DLMD_PARENT &
!              ,' global=',global
!
!!!   DPHD_CHILD=DPHD_PARENT*CHILD_PARENT_SPACE_RATIO
!!!   DLMD_CHILD=DLMD_PARENT*CHILD_PARENT_SPACE_RATIO
!
!!!   SBD_CHILD=-0.5*(JM_CHILD-1)*DPHD_CHILD
!!!   WBD_CHILD=-0.5*(IM_CHILD-1)*DLMD_CHILD
!
!!!   CALL CENTER_NEST(SBD_CHILD                                        &
!!!                   ,WBD_CHILD                                        &
!!!                   ,SW_LATD_CHILD                                    &
!!!                   ,SW_LOND_CHILD                                    &
!!!                   ,TPH0D_CHILD                                      &
!!!                   ,TLM0D_CHILD )
!
!-----------------------------------------------------------------------
!***  ALLOCATE 2-D AND 3-D DUMMY ARRAYS FOR CHILD QUANTITIES.
!-----------------------------------------------------------------------
!
      ALLOCATE(SEA_MASK(1:IM_CHILD,1:JM_CHILD))
      ALLOCATE(SEA_ICE(1:IM_CHILD,1:JM_CHILD))
!
      ALLOCATE(IDUMMY_2D(1:IM_CHILD,1:JM_CHILD))
      ALLOCATE(DUMMY_2D_IN (IMS:IME,JMS:JME,1:1))
      ALLOCATE(DUMMY_2D_OUT(1:IM_CHILD,1:JM_CHILD,1:1))
      ALLOCATE(DUMMY_3D (1:IM_CHILD,1:JM_CHILD,1:LM))
      ALLOCATE(DUMMY_3DS(1:IM_CHILD,1:JM_CHILD,1:NUM_SOIL_LAYERS))
      ALLOCATE(DUMMY_SOIL(1:NUM_SOIL_LAYERS))
      ALLOCATE(TEMPSOIL (1:NUM_SOIL_LAYERS,1:IM_CHILD,1:JM_CHILD))
!
      ALLOCATE(PD_NEAREST (1:IM_CHILD,1:JM_CHILD,1:1))
      ALLOCATE(PD_BILINEAR(1:IM_CHILD,1:JM_CHILD,1:1))
!
      ALLOCATE(LOWER_TOPO(1:IM_CHILD,1:JM_CHILD))
      DO J=1,JM_CHILD
      DO I=1,IM_CHILD
        LOWER_TOPO(I,J)=.FALSE.
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  PARENT TASK 0 OPENS A FILE FOR WRITING OUT THE CHILD'S INPUT.
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
!
        select_unit: DO N=51,59
          INQUIRE(N,OPENED=OPENED)
          IF(.NOT.OPENED)THEN
            NFCST=N
            EXIT select_unit
          ENDIF
        ENDDO select_unit
!
        FMT='(I2.2)'
        WRITE(INT_TO_CHAR,FMT)THIS_CHILD_ID
        OUTFILE='input_domain_'//INT_TO_CHAR
!
        OPEN(unit=NFCST,file=OUTFILE,status='replace',form='unformatted')
        write(0,*)' PARENT_CHILD_INIT opened file ',trim(outfile),' unit ',nfcst
!
!-----------------------------------------------------------------------
!***  THE FOLLOWING VARIABLES ARE FOR THE VERTICAL GRID STRUCTURE
!***  AND ARE SHARED BY THE PARENT AND ITS CHILDREN.
!-----------------------------------------------------------------------
!
        IHREND=0                                                           !<-- Not used 
        NTSD  =0                                                           !<-- Not used
!
        WRITE(NFCST)dyn_int_state%RUN                                   &
                   ,dyn_int_state%IDAT                                  &
                   ,dyn_int_state%IHRST                                 &
!                  ,dyn_int_state%IHREND                                &
!                  ,dyn_int_state%NTSD
                   ,IHREND                                              &
                   ,NTSD
!
        WRITE(NFCST)dyn_int_state%PT                                    &
                   ,dyn_int_state%PDTOP                                 &
                   ,dyn_int_state%LPT2                                  &
                   ,dyn_int_state%SGM                                   &
                   ,dyn_int_state%SG1                                   &
                   ,dyn_int_state%DSG1                                  &
                   ,dyn_int_state%SGML1                                 &
                   ,dyn_int_state%SG2                                   &
                   ,dyn_int_state%DSG2                                  &
                   ,dyn_int_state%SGML2
!
        WRITE(NFCST)I_PARENT_START,J_PARENT_START
!
        DLMD=DLMD_PARENT*CHILD_PARENT_SPACE_RATIO
        DPHD=DPHD_PARENT*CHILD_PARENT_SPACE_RATIO
!
        IF(GLOBAL)THEN
          SBD=SBD_PARENT+(J_PARENT_START-2)*DPHD_PARENT
          WBD=WBD_PARENT+(I_PARENT_START-2)*DLMD_PARENT
        ELSE 
          SBD=SBD_PARENT+(J_PARENT_START-1)*DPHD_PARENT
          WBD=WBD_PARENT+(I_PARENT_START-1)*DLMD_PARENT
        ENDIF
!
        WRITE(NFCST)DLMD,DPHD                                           &
                   ,WBD,SBD                                             &
                   ,TLM0D_PARENT,TPH0D_PARENT
!
        WRITE(NFCST)IM_CHILD,JM_CHILD,LM,LNSH
!
!-----------------------------------------------------------------------
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Sea Mask
!-----------------------------------------------------------------------
!***  The Sea Mask is needed for the Sfc Geopotential so compute it now.
!***  If there are adjacent water points with different elevations
!***  after Sfc Geopotential is computed then the WATERFALL routine
!***  will level them by changing the sfc elevations.  At such points
!***  the atmospheric column will need adjusting so save the locations
!***  of those points along with the preliminary values of the nest's
!***  PD, T, Q, CW, U, and V which will then be modified.  
!***  The Sea Mask will be written out in its proper place following
!***  the Stnd Deviation of Sfc Height.  
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=dyn_int_state%SM(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'SeaMask'                               &
                               ,DUMMY_2D_OUT                            &
                               ,NEAREST)
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          SEA_MASK(I,J)=REAL(NINT(DUMMY_2D_OUT(I,J,1)))
        ENDDO
        ENDDO
      ENDIF
!
!-----------------------------------------------------------------------
!***  Sfc Geopotential
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=dyn_int_state%FIS(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'FIS'                                   &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL WATERFALLS(DUMMY_2D_OUT                                    &  !<-- Level adjacent water points with different elevations
                       ,SEA_MASK                                        &
                       ,LOWER_TOPO                                      &
                       ,1,IM_CHILD,1,JM_CHILD)
!
        WRITE(NFCST)DUMMY_2D_OUT
      ENDIF
!
!     write(0,*)' after Sfc Geo'
!
!-----------------------------------------------------------------------
!***  Stnd Deviation of Sfc Height
!-----------------------------------------------------------------------

      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%STDH(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'STDH'                                  &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after STDH'
!
!-----------------------------------------------------------------------
!***  Sea Mask
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        WRITE(NFCST)SEA_MASK
      ENDIF
!     write(0,*)' after Sea Mask'
!
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=dyn_int_state%PD(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'PD'                                    &
                               ,PD_BILINEAR                             &
                               ,BILINEAR)
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &  !<-- Save nearest neighbors for topo adjustment
                               ,'PD'                                    &
                               ,PD_NEAREST                              &
                               ,NEAREST) 
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(LOWER_TOPO(I,J))THEN
            DUMMY_2D_OUT(I,J,1)=PD_NEAREST(I,J,1)
          ELSE
            DUMMY_2D_OUT(I,J,1)=PD_BILINEAR(I,J,1)
          ENDIF
        ENDDO
        ENDDO
      ENDIF
!
!     write(0,*)' after PD'
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!
!-----------------------------------------------------------------------
!***  U
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(dyn_int_state%U, LM                     &
                               ,'Uwind'                                 &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,dyn_int_state%PT                            &
                           ,dyn_int_state%PDTOP                         &
                           ,dyn_int_state%SG1                           &
                           ,dyn_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after U'
!
!-----------------------------------------------------------------------
!***  V
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(dyn_int_state%V, LM                     &
                               ,'Vwind'                                 &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,dyn_int_state%PT                            &
                           ,dyn_int_state%PDTOP                         &
                           ,dyn_int_state%SG1                           &
                           ,dyn_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after V'
!
!-----------------------------------------------------------------------
!***  T
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(dyn_int_state%T, LM                     &
                               ,'Temperature'                           &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,dyn_int_state%PT                            &
                           ,dyn_int_state%PDTOP                         &
                           ,dyn_int_state%SG1                           &
                           ,dyn_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after T'
!
!-----------------------------------------------------------------------
!***  Q
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(dyn_int_state%Q, LM                     &
                               ,'SpecHum'                               &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,dyn_int_state%PT                            &
                           ,dyn_int_state%PDTOP                         &
                           ,dyn_int_state%SG1                           &
                           ,dyn_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after Q'
!
!-----------------------------------------------------------------------
!***  CW
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(dyn_int_state%CW, LM                    &
                               ,'CW'                                    &
                               ,DUMMY_3D                                &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        CALL ADJUST_COLUMNS(PD_NEAREST                                  &
                           ,PD_BILINEAR                                 &
                           ,LOWER_TOPO                                  &
                           ,DUMMY_3D                                    &
                           ,dyn_int_state%PT                            &
                           ,dyn_int_state%PDTOP                         &
                           ,dyn_int_state%SG1                           &
                           ,dyn_int_state%SG2                           &
                           ,IM_CHILD,JM_CHILD)
      ENDIF
!
      DO L=1,LM
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after CW'
!
!-----------------------------------------------------------------------
!***  O3
!-----------------------------------------------------------------------
!
!     CALL PARENT_TO_CHILD_FILL(dyn_int_state%O3, LM                    &
!                              ,'O3'                                    &
!                              ,DUMMY_3D                                &
!                              ,BILINEAR)
!
!     IF(MYPE==0)THEN
!       CALL ADJUST_COLUMNS(PD_NEAREST                                  &
!                          ,PD_BILINEAR                                 &
!                          ,LOWER_TOPO                                  &
!                          ,DUMMY_3D                                    &
!                          ,dyn_int_state%PT                            &
!                          ,dyn_int_state%PDTOP                         &
!                          ,dyn_int_state%SG1                           &
!                          ,dyn_int_state%SG2                           &
!                          ,IM_CHILD,JM_CHILD)
!     ENDIF
!
      DO L=1,LM
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          DUMMY_3D(I,J,L)=0.    ! for now keep O3 = 0.
        ENDDO
        ENDDO
        IF(MYPE==0)WRITE(NFCST)DUMMY_3D(:,:,L)
      ENDDO
!     write(0,*)' after O3'
!
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%ALBEDO(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                           &
                               ,'ALBEDO'                                 &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after Albedo'
!
!-----------------------------------------------------------------------
!***  EPSR
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%EPSR(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                           &
                               ,'EPSR'                                   &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after EPSR'
!
!-----------------------------------------------------------------------
!***  MXSNAL
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%MXSNAL(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'MXSNAL'                                &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after MXSNAL'
!
!-----------------------------------------------------------------------
!***  TSKIN  
!-----------------------------------------------------------------------
!
!     write(0,*)' before TSKIN'
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%TSKIN(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'TSKIN'                                 &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(DUMMY_2D_OUT(I,J,1)<150.)THEN
            SEA_MASK(I,J)=1.0
            DUMMY_2D_OUT(I,J,1)=0.
          ENDIF
          if(dummy_2d_out(i,j,1)<173..and.sea_mask(i,j)<0.5)then
            write(0,*)' Very cold TSKIN=',dummy_2d_out(i,j,1) &
                     ,' at (',i,',',j,')'
          endif
        ENDDO
        ENDDO
      ENDIF
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after TSKIN'
!
!-----------------------------------------------------------------------
!***  SST
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%SST(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'SST'                                   &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SST'
!
!-----------------------------------------------------------------------
!***  SNO
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%SNO(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                           &
                               ,'SNO'                                    &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SNO'
!
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%SI(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                           &
                               ,'SI'                                     &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SI'
!
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%SICE(I,J)
      ENDDO
      ENDDO
!
!     write(0,*)' PARENT_TO_CHILD_INIT SICE max=',maxval(DUMMY_2D_IN(IMS:IME,JMS:JME,1)) &
!              ,' min=',minval(DUMMY_2D_IN(IMS:IME,JMS:JME,1))
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                           &
                               ,'SICE'                                   &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          SEA_ICE(I,J)=DUMMY_2D_OUT(I,J,1)
        ENDDO
        ENDDO
      ENDIF
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SICE'
!
!-----------------------------------------------------------------------
!***  TG  
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%TG(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'TG'                                    &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after TG'
!
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%CMC(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'CMC'                                   &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after CMC'
!
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%SR(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'SR'                                    &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after SR'
!
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%USTAR(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'USTAR'                                 &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after USTAR'
!
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%Z0(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'Z0'                                    &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after Z0'
!
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%Z0BASE(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'Z0BASE'                                &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after Z0BASE'
!
!-----------------------------------------------------------------------
!***  STC
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(phy_int_state%STC, NUM_SOIL_LAYERS      &
                               ,'STC'                                   &
                               ,DUMMY_3DS                               &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO L=1,NUM_SOIL_LAYERS
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          TEMPSOIL(L,I,J)=DUMMY_3DS(I,J,L)
        ENDDO
        ENDDO
        ENDDO
!
        WRITE(NFCST)TEMPSOIL
      ENDIF
!     write(0,*)' after STC'
!
!-----------------------------------------------------------------------
!***  SMC
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(phy_int_state%SMC, NUM_SOIL_LAYERS      &
                               ,'SMC'                                   &
                               ,DUMMY_3DS                               &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO L=1,NUM_SOIL_LAYERS
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          TEMPSOIL(L,I,J)=DUMMY_3DS(I,J,L)
        ENDDO
        ENDDO
        ENDDO
!
        WRITE(NFCST)TEMPSOIL
      ENDIF
!     write(0,*)' after SMC'
!
!-----------------------------------------------------------------------
!***  SH2O
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_FILL(phy_int_state%SH2O, NUM_SOIL_LAYERS     &
                               ,'SH2O'                                  &
                               ,DUMMY_3DS                               &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO L=1,NUM_SOIL_LAYERS
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          TEMPSOIL(L,I,J)=DUMMY_3DS(I,J,L)
        ENDDO
        ENDDO
        ENDDO
!
        WRITE(NFCST)TEMPSOIL
      ENDIF
!     write(0,*)' after SH2O'
!
!-----------------------------------------------------------------------
!***  ISLTYP
!-----------------------------------------------------------------------
!
      CALL PARENT_TO_CHILD_IFILL(phy_int_state%ISLTYP                   &
                                ,'ISLTYP'                               &
                                ,IDUMMY_2D )
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(IDUMMY_2D(I,J)<1.AND.SEA_MASK(I,J)<0.5)THEN
            IDUMMY_2D(I,J)=1     !<--------- Bandaid for interpolated soil value=0 while interpolated seamask=0 (i.e, a land point)
!           if(abs(IDUMMY_2D(I,J))>50)write(0,*)' write ISLTYP i=',i,' j=',j,' sea_mask=',SEA_MASK(I,J),' isltyp=',IDUMMY_2D(I,J)
          ENDIF
        ENDDO
        ENDDO
!
        WRITE(NFCST)IDUMMY_2D
      ENDIF
!     write(0,*)' after ISLTYP'
!
!-----------------------------------------------------------------------
!***  IVGTYP
!-----------------------------------------------------------------------
!
!     write(0,*)' PARENT_TO_CHILD_INIT IVGTYP max=',maxval(phy_int_state%IVGTYP) &
!              ,' min=',minval(phy_int_state%IVGTYP),' maxloc=',maxloc(phy_int_state%IVGTYP) &
!              ,' minloc=',minloc(phy_int_state%IVGTYP)
      CALL PARENT_TO_CHILD_IFILL(phy_int_state%IVGTYP                   &
                                ,'IVGTYP'                               &
                                ,IDUMMY_2D )
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(IDUMMY_2D(I,J)<1.AND.SEA_MASK(I,J)<0.5)THEN
            IDUMMY_2D(I,J)=1     !<--------- Bandaid for interpolated vegetation value=0 while interpolated seamask=0 (i.e, a land point)
          ENDIF
        ENDDO
        ENDDO
!
        WRITE(NFCST)IDUMMY_2D
      ENDIF
!     write(0,*)' after IVGTYP'
!
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        DUMMY_2D_IN(I,J,1)=phy_int_state%VEGFRC(I,J)
      ENDDO
      ENDDO
!
      CALL PARENT_TO_CHILD_FILL(DUMMY_2D_IN, 1                          &
                               ,'VEGFRC'                                &
                               ,DUMMY_2D_OUT                            &
                               ,BILINEAR)
!
      IF(MYPE==0)THEN
        DO J=1,JM_CHILD
        DO I=1,IM_CHILD
          IF(DUMMY_2D_OUT(I,J,1)>0..AND.(SEA_MASK(I,J)>0.5.OR.SEA_ICE(I,J)>0.))THEN
            DUMMY_2D_OUT(I,J,1)=0.    !<--------- Bandaid for interpolated veg frac value >0 while interpolated seamask or sice >0
          ENDIF
        ENDDO
        ENDDO
      ENDIF
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_2D_OUT
!     write(0,*)' after VEGFRC'
!
!-----------------------------------------------------------------------
!
      IF(MYPE==0)WRITE(NFCST)DUMMY_SOIL
      IF(MYPE==0)WRITE(NFCST)DUMMY_SOIL
!
!-----------------------------------------------------------------------
!
      IF(MYPE==0)CLOSE(NFCST)
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(IDUMMY_2D)
      DEALLOCATE(DUMMY_2D_IN)
      DEALLOCATE(DUMMY_2D_OUT)
      DEALLOCATE(DUMMY_3D)
      DEALLOCATE(DUMMY_3DS)
      DEALLOCATE(DUMMY_SOIL)
      DEALLOCATE(TEMPSOIL)
      DEALLOCATE(SEA_MASK)
      DEALLOCATE(PD_BILINEAR)
      DEALLOCATE(PD_NEAREST)
      DEALLOCATE(LOWER_TOPO)
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!!!   SUBROUTINE PARENT_TO_CHILD_FILL_ASSOC(PARENT_ARRAY                &
      SUBROUTINE PARENT_TO_CHILD_FILL      (PARENT_ARRAY                &
                                           ,NLEV                        &
                                           ,VBL_NAME                    &
                                           ,CHILD_ARRAY                 &
                                           ,METHOD)
!
!-----------------------------------------------------------------------
!***  Rows and columns of the child's grid lie directly on top of
!***  rows and colums of the parent (thus 'ASSOCIATED').
!
!***  Fill a child's domain with data from the parent.  Only the parent
!***  tasks are needed in this routine.
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: NLEV                                           !<-- Vertical dimension of the data array 
!
!!!   REAL,DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(INOUT) :: DATA_ARRAY   !<-- The parent array that will initialize the child array
      REAL,DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(IN) :: PARENT_ARRAY    !<-- The parent array that will initialize the child array
!
      CHARACTER(*),INTENT(IN) :: VBL_NAME
!
      REAL,DIMENSION(1:IM_CHILD,1:JM_CHILD,1:NLEV),INTENT(OUT) ::       &  !<-- Data from parent tasks interpolated to child grid
                                                           CHILD_ARRAY     !      but still on parent task 0
                                                                  
!
      INTEGER,INTENT(IN) :: METHOD                                         !<-- Interpolaton method (bilinear or nearest neighbor)
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,IERR,II,IPE,IPE_LOCAL,ISTAT,J,JJ,L,N,NN
      INTEGER :: I_COPY,I_END,I_END_COPY,I_EXTENT,I_PARENT_END          &
                ,I_START_COPY
      INTEGER :: J_COPY,J_END,J_END_COPY,J_EXTENT,J_PARENT_END          &
                ,J_START_COPY
      INTEGER :: INDX_EAST,INDX_NORTH,INDX_SOUTH,INDX_WEST
      INTEGER :: NWORDS_RECV,NWORDS_SEND
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL :: DELTA_I_EAST,DELTA_I_WEST,DELTA_J_NORTH,DELTA_J_SOUTH
      REAL :: RATIO,REAL_INDX_I_PARENT,REAL_INDX_J_PARENT
      REAL :: WEIGHT_EAST,WEIGHT_NORTH,WEIGHT_SOUTH,WEIGHT_WEST
      REAL :: WEIGHT_NE,WEIGHT_NW,WEIGHT_SE,WEIGHT_SW
      REAL :: WEIGHT_MAX,WEIGHT_SUM
!
      REAL,DIMENSION(:)    ,ALLOCATABLE :: DATA_BUFFER
      REAL,DIMENSION(:,:,:),ALLOCATABLE :: ARRAY_STAGE_PARENT    
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  To simplify matters somewhat, isolate the minimum subset of
!***  points on the parent domain that underlie the child's domain.
!
!***  The southwest corner of the child always lies directly on a
!***  point in the parent domain.  We already know the I,J of that
!***  parent point since it was specified in the configure file.
!***  The number of parent points that are covered by the child is
!***  determined by the child-to-parent grid ratio and the lateral
!***  dimensions of the child's domain.
!-----------------------------------------------------------------------
!
      I_PARENT_END=I_PARENT_START                                       &  !<-- Easternmost I on parent domain surrounding child domain
                     +INT((IM_CHILD-1)*CHILD_PARENT_SPACE_RATIO)+1
!
      I_EXTENT=I_PARENT_END-I_PARENT_START+1
!
      J_PARENT_END=J_PARENT_START                                       &  !<-- Northernmost J on parent domain surrounding child domain
                     +INT((JM_CHILD-1)*CHILD_PARENT_SPACE_RATIO)+1
!
      J_EXTENT=J_PARENT_END-J_PARENT_START+1
!
!-----------------------------------------------------------------------
!***  Create a staging array on parent task 0 that will hold the entire
!***  subset of the parent domain underlying the child.
!***  Then all parent tasks with points in the intersecting region
!***  send their data to parent task 0.
!-----------------------------------------------------------------------
!
      parent_stage: IF(MYPE==0)THEN                                        !<-- Parent task 0
! 
!-----------------------------------------------------------------------
!
        ALLOCATE(ARRAY_STAGE_PARENT(1:I_EXTENT,1:J_EXTENT,1:NLEV))         !<-- Array holding all parent points in staging region
                                                                           !    Note that this array begins at (1,1,1), i.e.,
                                                                           !      its indices are relative to the nest.
!
!-----------------------------------------------------------------------
!***  If parent task 0 holds some of the staging region, copy it to
!***  the staging array.
!-----------------------------------------------------------------------
!
        IF(I_PARENT_START<=ITE.AND.J_PARENT_START<=JTE)THEN
          I_END=MIN(ITE,I_PARENT_END)
          J_END=MIN(JTE,J_PARENT_END)
!
          DO L=1,NLEV
            JJ=0
            DO J=J_PARENT_START,J_END
              JJ=JJ+1
!
              II=0
              DO I=I_PARENT_START,I_END
                II=II+1     
                ARRAY_STAGE_PARENT(II,JJ,L)=PARENT_ARRAY(I,J,L)
              ENDDO
            ENDDO
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  If there are points in the staging region outside of parent task 0
!***  then task 0 receives those points from the other parent tasks that
!***  contain those points.
!-----------------------------------------------------------------------
!
        parent_search: DO IPE=1,NUM_PES_PARENT-1                           !<-- Parent task 0 checks other parent fcst tasks for points
!
          remote_stage: IF(I_PARENT_START<=LOCAL_IEND  (IPE).AND.       &  !<-- Does remote parent task IPE contain any staging region?
                           I_PARENT_END  >=LOCAL_ISTART(IPE)            &  !
                            .AND.                                       &  !
                           J_PARENT_START<=LOCAL_JEND  (IPE).AND.       &  !
                           J_PARENT_END  >=LOCAL_JSTART(IPE))THEN          !<--
! 
            I_START_COPY=MAX(I_PARENT_START,LOCAL_ISTART(IPE))             !<-- I index of first point in staging region on remote parent task
            I_END_COPY  =MIN(I_PARENT_END  ,LOCAL_IEND  (IPE))             !<-- I index of last point in staging region on remote parent task
            I_COPY      =I_END_COPY-I_START_COPY+1                         !<-- I range of points to receive
!
            J_START_COPY=MAX(J_PARENT_START,LOCAL_JSTART(IPE))             !<-- J index of first point in staging region on remote parent task
            J_END_COPY  =MIN(J_PARENT_END  ,LOCAL_JEND  (IPE))             !<-- J index of last point in staging region on remote parent task
            J_COPY      =J_END_COPY-J_START_COPY+1                         !<-- J range of points to receive
!
            NWORDS_RECV=I_COPY*J_COPY*NLEV                                 !<-- Total # of words from remote parent task in staging region
!
            ALLOCATE(DATA_BUFFER(1:NWORDS_RECV))                           !<-- Allocate buffer array to hold remote task's staging data
            CALL MPI_RECV(DATA_BUFFER                                   &  !<-- The staging region data from remote parent task IPE
                         ,NWORDS_RECV                                   &  !<-- Total words received
                         ,MPI_REAL                                      &  !<-- Datatype
                         ,IPE                                           &  !<-- Receive from this parent task
                         ,IPE                                           &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- The MPI communicator
                         ,JSTAT                                         &  !<-- MPI status object
                         ,IERR )
!
            NN=0                                                           !<-- Counter for received words
!
            DO L=1,NLEV
!
              JJ=J_START_COPY-J_PARENT_START
              DO J=1,J_COPY
                JJ=JJ+1
!
                II=I_START_COPY-I_PARENT_START
                DO I=1,I_COPY
                  II=II+1
                  NN=NN+1
                  ARRAY_STAGE_PARENT(II,JJ,L)=DATA_BUFFER(NN)              !<-- Fill in array with staging region data from parent task IPE
                ENDDO
              ENDDO
            ENDDO
!
            DEALLOCATE(DATA_BUFFER)
!
          ENDIF remote_stage
!
        ENDDO parent_search
!
!-----------------------------------------------------------------------
!***  Now the remaining parent tasks check to see if they contain
!***  any points in the staging region.  If they do, gather them
!***  and send them to parent task 0.
!-----------------------------------------------------------------------
!
      ELSEIF(MYPE>0.AND.MYPE<=NUM_PES_PARENT-1)THEN  parent_stage          !<-- All parent forecast tasks other than 0
!
!-----------------------------------------------------------------------
        IF(I_PARENT_START<=ITE.AND.I_PARENT_END>=ITS                    &  !<-- Does this parent task contain any staging region?
           .AND.                                                        &  !
           J_PARENT_START<=JTE.AND.J_PARENT_END>=JTS)THEN                  !<--
! 
          I_START_COPY=MAX(I_PARENT_START,ITS)                             !<-- I index of first point in staging region on this parent task
          I_END_COPY  =MIN(I_PARENT_END  ,ITE)                             !<-- I index of last point in staging region on this parent task
          I_COPY=I_END_COPY-I_START_COPY+1                                 !<-- I range of points to send to parent task 0
!
          J_START_COPY=MAX(J_PARENT_START,JTS)                             !<-- J index of first point in staging region on remote parent task
          J_END_COPY  =MIN(J_PARENT_END  ,JTE)                             !<-- J index of last point in staging region on remote parent task
          J_COPY=J_END_COPY-J_START_COPY+1                                 !<-- J range of copied points
!
          NWORDS_SEND=I_COPY*J_COPY*NLEV                                   !<-- Total number of words from this parent task in staging region
          ALLOCATE(DATA_BUFFER(1:NWORDS_SEND),stat=ISTAT)                  !<-- Allocate the buffer array to hold this task's staging data
!
          NN=0
!
          DO L=1,NLEV
          DO J=J_START_COPY,J_END_COPY
          DO I=I_START_COPY,I_END_COPY
            NN=NN+1
            DATA_BUFFER(NN)=PARENT_ARRAY(I,J,L)
          ENDDO
          ENDDO
          ENDDO
!
          CALL MPI_SEND(DATA_BUFFER                                     &  !<-- The staging region data from this parent task to parent task 0
                       ,NWORDS_SEND                                     &  !<-- Total words sent
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
          DEALLOCATE(DATA_BUFFER)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF parent_stage
!
!-----------------------------------------------------------------------
!***  The subset of the input array on the parent domain that lies 
!***  under the child's domain has been mirrored onto parent task 0.
!***  Parent task 0 will fill out the array to match the child
!***  domain's horizontal grid increments and then parcel out the
!***  appropriate pieces to the corresponding tasks of the child.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First fill in the southern and western sides of the array.
!***  If bilinear interpolation is specified then only linear 
!***  interpolation needs to be used.
!-----------------------------------------------------------------------
!
      parent_task_0: IF(MYPE==0)THEN                                      !<-- Parent task 0
!
!-----------------------------------------------------------------------
!
        RATIO=CHILD_PARENT_SPACE_RATIO
!
        DO L=1,NLEV
!
          CHILD_ARRAY(1,1,L)=ARRAY_STAGE_PARENT(1,1,L)                    !<-- SW corner of child's array coincides with a parent point
!
          DO I=2,IM_CHILD                                                 !<-- Move along southern boundary of child's domain
            REAL_INDX_I_PARENT=1+(I-1)*RATIO                              !<-- Exact I index of child point on parent 
            INDX_WEST=INT(REAL_INDX_I_PARENT)                             !<-- The parent point's I index west of the child's point
            INDX_EAST=INDX_WEST+1                                         !<-- The parent point's I index east of the child's point
            WEIGHT_WEST=INDX_EAST-REAL_INDX_I_PARENT                      !<-- Interpolation weight given parent's point to the west
            WEIGHT_EAST=1.-WEIGHT_WEST                                    !<-- Interpolation weight given parent's point to the east
!
            IF(METHOD==NEAREST)THEN                                       !<-- Assign points using nearest neighbors
              WEIGHT_MAX=MAX(WEIGHT_WEST,WEIGHT_EAST)
              IF(WEIGHT_WEST==WEIGHT_MAX)THEN
                CHILD_ARRAY(I,1,L)=ARRAY_STAGE_PARENT(INDX_WEST,1,L)
              ELSEIF(WEIGHT_EAST==WEIGHT_MAX)THEN
                CHILD_ARRAY(I,1,L)=ARRAY_STAGE_PARENT(INDX_EAST,1,L)
              ENDIF
!
            ELSEIF(METHOD==BILINEAR)THEN                                         !<-- Assign points using (bi)linear interpolation
              CHILD_ARRAY(I,1,L)=WEIGHT_WEST*ARRAY_STAGE_PARENT(INDX_WEST,1,L) & !<-- Value at points along child's southern boundary
                                +WEIGHT_EAST*ARRAY_STAGE_PARENT(INDX_EAST,1,L)
            ELSE
              WRITE(0,*)" Attempting to use unknown interpolation method: ",METHOD
              CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
!
            ENDIF
!
          ENDDO
!
          DO J=2,JM_CHILD                                                 !<-- Move along western boundary of child's domain
            REAL_INDX_J_PARENT=1+(J-1)*RATIO                              !<-- Exact J index of child point on parent 
            INDX_SOUTH=INT(REAL_INDX_J_PARENT)                            !<-- The parent point's J index south of the child's point
            INDX_NORTH=INDX_SOUTH+1                                       !<-- The parent point's J index north of the child's point
            WEIGHT_SOUTH=INDX_NORTH-REAL_INDX_J_PARENT                    !<-- Interpolation weight of parent's point to the south
            WEIGHT_NORTH=1.-WEIGHT_SOUTH                                  !<-- Interpolation weight of parent's point to the north
!
            IF(METHOD==NEAREST)THEN                                       !<-- Assign points using nearest neighbors
              WEIGHT_MAX=MAX(WEIGHT_SOUTH,WEIGHT_NORTH)
              IF(WEIGHT_SOUTH==WEIGHT_MAX)THEN
                CHILD_ARRAY(1,J,L)=ARRAY_STAGE_PARENT(1,INDX_SOUTH,L)
              ELSE
                CHILD_ARRAY(1,J,L)=ARRAY_STAGE_PARENT(1,INDX_NORTH,L)
              ENDIF
!
            ELSEIF(METHOD==BILINEAR)THEN                                           !<-- Assign points using (bi)linear interpolation
              CHILD_ARRAY(1,J,L)=WEIGHT_SOUTH*ARRAY_STAGE_PARENT(1,INDX_SOUTH,L) & !<-- Value at points along child's western boundary
                              +WEIGHT_NORTH*ARRAY_STAGE_PARENT(1,INDX_NORTH,L)
            ELSE
              WRITE(0,*)" Attempting to use unknown interpolation method: ",METHOD
              CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
!
            ENDIF
!
          ENDDO
!
!-----------------------------------------------------------------------
!***  Fill in the interior of the staging array.
!-----------------------------------------------------------------------
!
          DO J=2,JM_CHILD
            REAL_INDX_J_PARENT=1+(J-1)*RATIO                              !<-- Exact J index of child point in parent staging region 
            INDX_SOUTH=INT(REAL_INDX_J_PARENT)                            !<-- The parent point's J index south of the child's point
            INDX_NORTH=INDX_SOUTH+1                                       !<-- The parent point's J index north of the child's point
!
            DELTA_J_NORTH=INDX_NORTH-REAL_INDX_J_PARENT                   !<-- Parent grid increment from child point to parent point north
            DELTA_J_SOUTH=REAL_INDX_J_PARENT-INDX_SOUTH                   !<-- Parent grid increment from child point to parent point south
!
            DO I=2,IM_CHILD
              REAL_INDX_I_PARENT=1+(I-1)*RATIO                            !<-- Exact I index of child point in parent staging region
              INDX_WEST=INT(REAL_INDX_I_PARENT)                           !<-- The parent point's I index west of the child's point
              INDX_EAST=INDX_WEST+1                                       !<-- The parent point's I index east of the child's point
!
              DELTA_I_EAST=INDX_EAST-REAL_INDX_I_PARENT
              DELTA_I_WEST=REAL_INDX_I_PARENT-INDX_WEST
!
              WEIGHT_SW=DELTA_I_EAST*DELTA_J_NORTH                        !<-- Interpolation weight of parent's point to SW 
              WEIGHT_SE=DELTA_I_WEST*DELTA_J_NORTH                        !<-- Interpolation weight of parent's point to SE
              WEIGHT_NW=DELTA_I_EAST*DELTA_J_SOUTH                        !<-- Interpolation weight of parent's point to NW
              WEIGHT_NE=DELTA_I_WEST*DELTA_J_SOUTH                        !<-- Interpolation weight of parent's point to NE
!
!-----------------------------------------------------------------------
!
              assign: IF(METHOD==NEAREST)THEN                             !<-- Assign points using nearest neighbors      
                WEIGHT_MAX=MAX(WEIGHT_SW,WEIGHT_SE                     &
                              ,WEIGHT_NW,WEIGHT_NE)
                IF(WEIGHT_SW==WEIGHT_MAX)THEN
                  CHILD_ARRAY(I,J,L)=ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L)
                ELSEIF(WEIGHT_SE==WEIGHT_MAX)THEN
                  CHILD_ARRAY(I,J,L)=ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L)
                ELSEIF(WEIGHT_NW==WEIGHT_MAX)THEN
                  CHILD_ARRAY(I,J,L)=ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L)
                ELSEIF(WEIGHT_NE==WEIGHT_MAX)THEN
                  CHILD_ARRAY(I,J,L)=ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L)
                ENDIF
!
              ELSEIF(METHOD==BILINEAR)THEN                                                 !<-- Assign points using bilinear interpolation
                IF(VBL_NAME/='FIS')THEN
                  IF(ABS(ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L))<1.E-12)WEIGHT_SW=0.
                  IF(ABS(ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L))<1.E-12)WEIGHT_SE=0.
                  IF(ABS(ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L))<1.E-12)WEIGHT_NW=0.
                  IF(ABS(ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L))<1.E-12)WEIGHT_NE=0.
                ENDIF
!
                CHILD_ARRAY(I,J,L)=WEIGHT_SW*ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L)  & !<-- Value at points in child's interior
                                  +WEIGHT_SE*ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L)  &
                                  +WEIGHT_NW*ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L)  &
                                  +WEIGHT_NE*ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L)
!
                WEIGHT_SUM=WEIGHT_SW+WEIGHT_SE+WEIGHT_NW+WEIGHT_NE
                IF(WEIGHT_SUM<0.99.AND.WEIGHT_SUM>0.01)THEN
                  CHILD_ARRAY(I,J,L)=CHILD_ARRAY(I,J,L)/WEIGHT_SUM        !<-- Normalize if some weights are zero (e.g., coastal land Temp)
                ENDIF
!
                IF(VBL_NAME=='SST')THEN                                   !<-- Include only realistic SST temperatures
                  WEIGHT_SUM=0.
                  CHILD_ARRAY(I,J,L)=0.
                  IF(ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L)>200.)THEN
                    WEIGHT_SUM=WEIGHT_SUM+WEIGHT_SW
                    CHILD_ARRAY(I,J,L)=WEIGHT_SW*ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH,L) &
                                       +CHILD_ARRAY(I,J,L)
                  ENDIF
                  IF(ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L)>200.)THEN
                    WEIGHT_SUM=WEIGHT_SUM+WEIGHT_SE
                    CHILD_ARRAY(I,J,L)=WEIGHT_SE*ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH,L) &
                                       +CHILD_ARRAY(I,J,L)
                  ENDIF
                  IF(ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L)>200.)THEN
                    WEIGHT_SUM=WEIGHT_SUM+WEIGHT_NW
                    CHILD_ARRAY(I,J,L)=WEIGHT_NW*ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH,L) &
                                       +CHILD_ARRAY(I,J,L)
                  ENDIF
                  IF(ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L)>200.)THEN
                    WEIGHT_SUM=WEIGHT_SUM+WEIGHT_NE
                    CHILD_ARRAY(I,J,L)=WEIGHT_NE*ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH,L) &
                                       +CHILD_ARRAY(I,J,L)
                  ENDIF
                  IF(WEIGHT_SUM<0.99.AND.WEIGHT_SUM>0.01)THEN
                    CHILD_ARRAY(I,J,L)=CHILD_ARRAY(I,J,L)/WEIGHT_SUM
                  ENDIF
                ENDIF
!
              ELSE
                WRITE(0,*)" Attempting to use unknown interpolation method: ",METHOD
                CALL ESMF_Finalize(terminationflag=ESMF_ABORT)
!
              ENDIF assign
!
!-----------------------------------------------------------------------
!
            ENDDO
          ENDDO
!
        ENDDO
!
        DEALLOCATE(ARRAY_STAGE_PARENT)
!
!-----------------------------------------------------------------------
!
      ENDIF parent_task_0
!
!-----------------------------------------------------------------------
!
!!!   END SUBROUTINE PARENT_TO_CHILD_FILL_ASSOC
      END SUBROUTINE PARENT_TO_CHILD_FILL
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!!!   SUBROUTINE PARENT_TO_CHILD_IFILL_ASSOC(PARENT_ARRAY               &
      SUBROUTINE PARENT_TO_CHILD_IFILL      (PARENT_ARRAY               &
                                            ,VBL_NAME                   &
                                            ,CHILD_ARRAY)
!
!-----------------------------------------------------------------------
!***  Rows and columns of the child's grid lie directly on top of
!***  rows and colums of the parent (thus 'ASSOCIATED').
!***  Fill a child's domain with data from the parent.  Only the parent
!***  tasks are needed in this routine.
!-----------------------------------------------------------------------
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PARENT_ARRAY        !<-- The parent array that will initialize the child array
!
      CHARACTER(*)                      ,INTENT(IN) :: VBL_NAME            !<-- The variable's name
!
      INTEGER,DIMENSION(1:IM_CHILD,1:JM_CHILD),INTENT(OUT) ::           &  !<-- Data from parent tasks interpolated to child grid
                                                          CHILD_ARRAY      !      but still on parent task 0
                                                                  
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,IERR,II,IPE,IPE_LOCAL,J,JJ,N,NN
      INTEGER :: I_COPY,I_END,I_END_COPY,I_EXTENT,I_PARENT_END          &
                ,I_START_COPY
      INTEGER :: J_COPY,J_END,J_END_COPY,J_EXTENT,J_PARENT_END          &
                ,J_START_COPY
      INTEGER :: INDX_EAST,INDX_NORTH,INDX_SOUTH,INDX_WEST
      INTEGER :: NWORDS_RECV,NWORDS_SEND
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL :: DELTA_I_EAST,DELTA_I_WEST,DELTA_J_NORTH,DELTA_J_SOUTH
      REAL :: RATIO,REAL_INDX_I_PARENT,REAL_INDX_J_PARENT
      REAL :: WEIGHT_EAST,WEIGHT_NORTH,WEIGHT_SOUTH,WEIGHT_WEST
      REAL :: WEIGHT_NE,WEIGHT_NW,WEIGHT_SE,WEIGHT_SW
      REAL :: WEIGHT_MAX
!
      INTEGER,DIMENSION(:)  ,ALLOCATABLE :: DATA_BUFFER
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ARRAY_STAGE_PARENT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  To simplify matters somewhat, isolate the minimum subset of
!***  points on the parent domain that underlie the child's domain.
!
!***  The southwest corner of the child always lies directly on a
!***  point in the parent domain.  We already know the I,J of that
!***  parent point since it was specified in the configure file.
!***  The number of parent points that are covered by the child is
!***  determined by the child-to-parent grid ratio and the lateral
!***  dimensions of the child's domain.
!-----------------------------------------------------------------------
!
      I_PARENT_END=I_PARENT_START                                       &  !<-- Easternmost I on parent domain surrounding child domain
                     +INT((IM_CHILD-1)*CHILD_PARENT_SPACE_RATIO)+1
!
      I_EXTENT=I_PARENT_END-I_PARENT_START+1
!
      J_PARENT_END=J_PARENT_START                                       &  !<-- Northernmost J on parent domain surrounding child domain
                     +INT((JM_CHILD-1)*CHILD_PARENT_SPACE_RATIO)+1
!
      J_EXTENT=J_PARENT_END-J_PARENT_START+1
!
!-----------------------------------------------------------------------
!***  Create a staging array on parent task 0 that will hold the entire
!***  subset of the parent domain underlying the child.
!***  Then all parent tasks with points in the intersecting region
!***  send their data to parent task 0.
!-----------------------------------------------------------------------
!
      parent_stage: IF(MYPE==0)THEN                                        !<-- Parent task 0
! 
!-----------------------------------------------------------------------
!
        ALLOCATE(ARRAY_STAGE_PARENT(1:I_EXTENT,1:J_EXTENT))                !<-- Array holding all parent points in staging region
                                                                           !    Note that this array begins at (1,1,1), i.e.,
                                                                           !      its indices are relative to the nest.
!
!-----------------------------------------------------------------------
!***  If parent task 0 holds some of the staging region, copy it to
!***  the staging array.
!-----------------------------------------------------------------------
!
        IF(I_PARENT_START<=ITE.AND.J_PARENT_START<=JTE)THEN
          I_END=MIN(ITE,I_PARENT_END)
          J_END=MIN(JTE,J_PARENT_END)
!
          JJ=0
          DO J=J_PARENT_START,J_END
            JJ=JJ+1
!
            II=0
            DO I=I_PARENT_START,I_END
              II=II+1     
              ARRAY_STAGE_PARENT(II,JJ)=PARENT_ARRAY(I,J)
            ENDDO
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!***  If there are points in the staging region outside of parent task 0
!***  then task 0 receives those points from the other parent tasks that
!***  contain those points.
!-----------------------------------------------------------------------
!
        parent_search: DO IPE=1,NUM_PES_PARENT-1                           !<-- Parent task 0 checks other parent tasks for points
!
          remote_stage: IF(I_PARENT_START<=LOCAL_IEND  (IPE).AND.       &  !<-- Does remote parent task IPE contain any staging region?
                           I_PARENT_END  >=LOCAL_ISTART(IPE)            &  !
                            .AND.                                       &  !
                           J_PARENT_START<=LOCAL_JEND  (IPE).AND.       &  !
                           J_PARENT_END  >=LOCAL_JSTART(IPE))THEN          !<--
! 
            I_START_COPY=MAX(I_PARENT_START,LOCAL_ISTART(IPE))             !<-- I index of first point in staging region on remote parent task
            I_END_COPY  =MIN(I_PARENT_END  ,LOCAL_IEND  (IPE))             !<-- I index of last point in staging region on remote parent task
            I_COPY      =I_END_COPY-I_START_COPY+1                         !<-- I range of points to receive
!
            J_START_COPY=MAX(J_PARENT_START,LOCAL_JSTART(IPE))             !<-- J index of first point in staging region on remote parent task
            J_END_COPY  =MIN(J_PARENT_END  ,LOCAL_JEND  (IPE))             !<-- J index of last point in staging region on remote parent task
            J_COPY      =J_END_COPY-J_START_COPY+1                         !<-- J range of points to receive
!
            NWORDS_RECV=I_COPY*J_COPY                                      !<-- Total # of words from remote parent task in staging region
!
            ALLOCATE(DATA_BUFFER(1:NWORDS_RECV))                           !<-- Allocate buffer array to hold remote task's staging data
            CALL MPI_RECV(DATA_BUFFER                                   &  !<-- The staging region data from remote parent task IPE
                         ,NWORDS_RECV                                   &  !<-- Total words received
                         ,MPI_INTEGER                                   &  !<-- Datatype
                         ,IPE                                           &  !<-- Receive from this parent task
                         ,IPE                                           &  !<-- MPI tag
                         ,COMM_MY_DOMAIN                                &  !<-- The MPI communicator
                         ,JSTAT                                         &  !<-- MPI status object
                         ,IERR )
!
            NN=0                                                           !<-- Counter for received words
!
              JJ=J_START_COPY-J_PARENT_START
              DO J=1,J_COPY
                JJ=JJ+1
!
                II=I_START_COPY-I_PARENT_START
                DO I=1,I_COPY
                  II=II+1
                  NN=NN+1
                  ARRAY_STAGE_PARENT(II,JJ)=DATA_BUFFER(NN)                !<-- Fill in array with staging region data from parent task IPE
                ENDDO
              ENDDO
!
            DEALLOCATE(DATA_BUFFER)
!
          ENDIF remote_stage
!
        ENDDO parent_search
!
!-----------------------------------------------------------------------
!***  Now the remaining parent tasks check to see if they contain
!***  any points in the staging region.  If they do, gather them
!***  and send them to parent task 0.
!-----------------------------------------------------------------------
!
      ELSEIF(MYPE>0.AND.MYPE<=NUM_PES_PARENT-1)THEN  parent_stage          !<-- All parent tasks other than 0
!
!-----------------------------------------------------------------------
        IF(I_PARENT_START<=ITE.AND.I_PARENT_END>=ITS                    &  !<-- Does this parent task contain any staging region?
           .AND.                                                        &  !
           J_PARENT_START<=JTE.AND.J_PARENT_END>=JTS)THEN                  !<--
! 
          I_START_COPY=MAX(I_PARENT_START,ITS)                             !<-- I index of first point in staging region on this parent task
          I_END_COPY  =MIN(I_PARENT_END  ,ITE)                             !<-- I index of last point in staging region on this parent task
          I_COPY=I_END_COPY-I_START_COPY+1                                 !<-- I range of points to send to parent task 0
!
          J_START_COPY=MAX(J_PARENT_START,JTS)                             !<-- J index of first point in staging region on remote parent task
          J_END_COPY  =MIN(J_PARENT_END  ,JTE)                             !<-- J index of last point in staging region on remote parent task
          J_COPY=J_END_COPY-J_START_COPY+1                                 !<-- J range of copied points
!
          NWORDS_SEND=I_COPY*J_COPY                                        !<-- Total number of words from this parent task in staging region
          ALLOCATE(DATA_BUFFER(1:NWORDS_SEND))                             !<-- Allocate the buffer array to hold this task's staging data
!
          NN=0
!
          DO J=J_START_COPY,J_END_COPY
          DO I=I_START_COPY,I_END_COPY
            NN=NN+1
            DATA_BUFFER(NN)=PARENT_ARRAY(I,J)
          ENDDO
          ENDDO
!
          CALL MPI_SEND(DATA_BUFFER                                     &  !<-- The staging region data from this parent task to parent task 0
                       ,NWORDS_SEND                                     &  !<-- Total words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
          DEALLOCATE(DATA_BUFFER)
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF parent_stage
!
!-----------------------------------------------------------------------
!***  The subset of the input array on the parent domain that lies 
!***  under the child's domain has been mirrored onto parent task 0.
!***  Parent task 0 will fill out the array to match the child
!***  domain's horizontal grid increments and then parcel out the
!***  appropriate pieces to the corresponding tasks of the child.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  First fill in the southern and western sides of the array
!***  choosing the nearest parent points.
!-----------------------------------------------------------------------
!
      parent_task_0: IF(MYPE==0)THEN                                      !<-- Parent task 0
!
!-----------------------------------------------------------------------
!
        RATIO=CHILD_PARENT_SPACE_RATIO
!
        CHILD_ARRAY(1,1)=ARRAY_STAGE_PARENT(1,1)                          !<-- SW corner of child's array coincides with a parent point
!
!***  Choose nearest parent point along child's southern boundary  
!
        DO I=2,IM_CHILD                                                   !<-- Move along southern boundary of child's domain
          REAL_INDX_I_PARENT=1+(I-1)*RATIO                                !<-- Exact I index of child point on parent 
          INDX_WEST=INT(REAL_INDX_I_PARENT)                               !<-- The parent point's I index west of the child's point
          INDX_EAST=INDX_WEST+1                                           !<-- The parent point's I index east of the child's point
          WEIGHT_WEST=INDX_EAST-REAL_INDX_I_PARENT                        !<-- Interpolation weight given parent's point to the west
          WEIGHT_EAST=1.-WEIGHT_WEST                                      !<-- Interpolation weight given parent's point to the east
!
          WEIGHT_MAX=MAX(WEIGHT_WEST,WEIGHT_EAST)
          IF(WEIGHT_WEST==WEIGHT_MAX)THEN
            CHILD_ARRAY(I,1)=ARRAY_STAGE_PARENT(INDX_WEST,1)
          ELSEIF(WEIGHT_EAST==WEIGHT_MAX)THEN
            CHILD_ARRAY(I,1)=ARRAY_STAGE_PARENT(INDX_EAST,1)
          ENDIF
        ENDDO
!
!***  Choose nearest parent point along child's western boundary
!
        DO J=2,JM_CHILD                                                   !<-- Move along western boundary of child's domain
          REAL_INDX_J_PARENT=1+(J-1)*RATIO                                !<-- Exact J index of child point on parent 
          INDX_SOUTH=INT(REAL_INDX_J_PARENT)                              !<-- The parent point's J index south of the child's point
          INDX_NORTH=INDX_SOUTH+1                                         !<-- The parent point's J index north of the child's point
          WEIGHT_SOUTH=INDX_NORTH-REAL_INDX_J_PARENT                      !<-- Interpolation weight of parent's point to the south
          WEIGHT_NORTH=1.-WEIGHT_SOUTH                                       !<-- Interpolation weight of parent's point to the north
          WEIGHT_MAX=MAX(WEIGHT_SOUTH,WEIGHT_NORTH)
!
          IF(WEIGHT_SOUTH==WEIGHT_MAX)THEN
            CHILD_ARRAY(1,J)=ARRAY_STAGE_PARENT(1,INDX_SOUTH)
          ELSE
            CHILD_ARRAY(1,J)=ARRAY_STAGE_PARENT(1,INDX_NORTH)
          ENDIF
!
        ENDDO
!
!-----------------------------------------------------------------------
!***  Fill in the interior of the staging array choosing the
!***  nearest parent point.
!-----------------------------------------------------------------------
!
        DO J=2,JM_CHILD
          REAL_INDX_J_PARENT=1+(J-1)*RATIO                                !<-- Exact J index of child point in parent staging region 
          INDX_SOUTH=INT(REAL_INDX_J_PARENT)                              !<-- The parent point's J index south of the child's point
          INDX_NORTH=INDX_SOUTH+1                                         !<-- The parent point's J index north of the child's point
!
          DELTA_J_NORTH=INDX_NORTH-REAL_INDX_J_PARENT
          DELTA_J_SOUTH=REAL_INDX_J_PARENT-INDX_SOUTH
!
          DO I=2,IM_CHILD
            REAL_INDX_I_PARENT=1+(I-1)*RATIO                              !<-- Exact I index of child point in parent staging region
            INDX_WEST=INT(REAL_INDX_I_PARENT)                             !<-- The parent point's I index west of the child's point
            INDX_EAST=INDX_WEST+1                                         !<-- The parent point's I index east of the child's point
!
            DELTA_I_EAST=INDX_EAST-REAL_INDX_I_PARENT
            DELTA_I_WEST=REAL_INDX_I_PARENT-INDX_WEST
!
            WEIGHT_SW=DELTA_I_EAST*DELTA_J_NORTH                          !<-- Interpolation weight of parent's point to SW 
            WEIGHT_SE=DELTA_I_WEST*DELTA_J_NORTH                          !<-- Interpolation weight of parent's point to SE
            WEIGHT_NW=DELTA_I_EAST*DELTA_J_SOUTH                          !<-- Interpolation weight of parent's point to NW
            WEIGHT_NE=DELTA_I_WEST*DELTA_J_SOUTH                          !<-- Interpolation weight of parent's point to NE
!
            WEIGHT_MAX=MAX(WEIGHT_SW,WEIGHT_SE                         & 
                          ,WEIGHT_NW,WEIGHT_NE)
!
            IF(WEIGHT_SW==WEIGHT_MAX)THEN
              CHILD_ARRAY(I,J)=ARRAY_STAGE_PARENT(INDX_WEST,INDX_SOUTH)
            ELSEIF(WEIGHT_SE==WEIGHT_MAX)THEN
              CHILD_ARRAY(I,J)=ARRAY_STAGE_PARENT(INDX_EAST,INDX_SOUTH)
            ELSEIF(WEIGHT_NW==WEIGHT_MAX)THEN
              CHILD_ARRAY(I,J)=ARRAY_STAGE_PARENT(INDX_WEST,INDX_NORTH)
            ELSEIF(WEIGHT_NE==WEIGHT_MAX)THEN
              CHILD_ARRAY(I,J)=ARRAY_STAGE_PARENT(INDX_EAST,INDX_NORTH)
            ENDIF
          ENDDO
        ENDDO
!
        DEALLOCATE(ARRAY_STAGE_PARENT)
!
!-----------------------------------------------------------------------
!
      ENDIF parent_task_0
!
!-----------------------------------------------------------------------
!
!!!   END SUBROUTINE PARENT_TO_CHILD_IFILL_ASSOC
      END SUBROUTINE PARENT_TO_CHILD_IFILL
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_TO_CHILD_FILL_GENERAL(PARENT_ARRAY              &
!!!   SUBROUTINE PARENT_TO_CHILD_FILL        (PARENT_ARRAY              &
                                             ,NLEV                      &
                                             ,VBL_NAME                  &
                                             ,CHILD_ARRAY)
!
!-----------------------------------------------------------------------
!***  Parent tasks interpolate their data to the locations of their
!***  children's gridpoints.  The child grids are unique rotated
!***  lat/lon grids with their own centers.  The southwest H point
!***  of the child grid lies directly on an H point of the parent. 
!
!***  Only parent tasks participate in this work.
!-----------------------------------------------------------------------
!
!---------------
!***  Arguments
!---------------
!
      INTEGER,INTENT(IN) :: NLEV                                           !<-- Vertical dimension of the data array 
!
      REAL,DIMENSION(IMS:IME,JMS:JME,1:NLEV),INTENT(IN) :: PARENT_ARRAY    !<-- The parent array that will initialize the child array
!
      CHARACTER(*),INTENT(IN) :: VBL_NAME                                  !<-- The variable's name 
!
      REAL,DIMENSION(1:IM_CHILD,1:JM_CHILD,1:NLEV),INTENT(OUT) :: CHILD_ARRAY  !<-- Data from parent tasks interpolated to child grid
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(KIND=KINT) :: I,I_END,ISTART,J,J_END,JSTART               &
                           ,KOUNT,L,NIJ,NN,NTOT                         &
                           ,NUM_DATA                                    &
                           ,NUM_CHILD_POINTS                            &
                           ,NUM_IJ                                      &
                           ,NUM_POINTS_REMOTE
!
      INTEGER(KIND=KINT) :: I_PARENT_SW,I_PARENT_SE                     &
                           ,I_PARENT_NW,I_PARENT_NE                     &
                           ,J_PARENT_SW,J_PARENT_SE                     &
                           ,J_PARENT_NW,J_PARENT_NE
!
      INTEGER(KIND=KINT) :: IERR,IPE,ISTAT
!
      INTEGER(KIND=KINT),DIMENSION(:),ALLOCATABLE :: CHILD_POINT_INDICES &
                                                    ,IJ_REMOTE
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(KIND=KFPT) :: CHILD_LATD_ON_PARENT                           &
                        ,CHILD_LOND_ON_PARENT                           &
                        ,DIST                                           &
                        ,R_DLMD,R_DPHD                                  &
                        ,REAL_I_PARENT                                  &
                        ,REAL_J_PARENT                                  &
                        ,RLATD_SW,RLOND_SW                              &
                        ,RLATD_SE,RLOND_SE                              &
                        ,RLATD_NW,RLOND_NW                              &
                        ,RLATD_NE,RLOND_NE                              &
                        ,SUM,SUM_RECIP                                  &
                        ,WEIGHT_SW,WEIGHT_SE                            &
                        ,WEIGHT_NW,WEIGHT_NE                            &
                        ,WEIGHT_SUM,WEIGHT_SUM_RECIP
!
      REAL(KIND=KFPT),DIMENSION(4) :: RLATD,RLOND,WGT
!
      REAL(KIND=KFPT),DIMENSION(:),ALLOCATABLE :: CHILD_STRING          &
                                                 ,DATA_REMOTE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!     write(0,*)' enter PARENT_TO_CHILD_FILL_GENERAL'
!
      R_DPHD=1./DPHD_PARENT
      R_DLMD=1./DLMD_PARENT
!
      NUM_CHILD_POINTS=0
!
      ISTART=1
      JSTART=1
      IF(GLOBAL)THEN
        ISTART=2
        JSTART=2
      ENDIF
!
!-----------------------------------------------------------------------
!***  Each parent task is responsible for searching the parent domain 
!***  extending from ITS and JTS to ITE+1 and JTE+1.  Those latter +1's
!***  are needed in order to reach the next gridpoint in each direction.
!***  We cannot go outside the full domain of course plus the wind 
!***  points have no values at IDE and JDE.
!-----------------------------------------------------------------------
!
      IF(VBL_NAME=='Uwind'.OR.VBL_NAME=='Vwind')THEN
        I_END=MIN(ITE+1,IDE-1)
        J_END=MIN(JTE+1,JDE-1)
      ELSE
        I_END=MIN(ITE+1,IDE)
        J_END=MIN(JTE+1,JDE)
      ENDIF
!
!-----------------------------------------------------------------------
!
      NTOT=2*IM_CHILD*JM_CHILD
      ALLOCATE(CHILD_POINT_INDICES(1:NTOT),stat=ISTAT)
      IF(ISTAT/=0)WRITE(0,*)' Failed to allocate CHILD_POINT_INDICES in PARENT_TO_CHILD_FILL_GENERAL'
!
      NTOT=IM_CHILD*JM_CHILD*NLEV
      ALLOCATE(CHILD_STRING(1:NTOT),stat=ISTAT)
      IF(ISTAT/=0)WRITE(0,*)' Failed to allocate CHILD_STRING in PARENT_TO_CHILD_FILL_GENERAL'
!
!-----------------------------------------------------------------------
!***  Compute the parent's lat/lon of each child gridpoint in order to
!***  determine if that gridpoint lies on a given parent task. 
!***  Save the I,J of each gridpoint found since ultimately parent
!***  task 0 will need that to properly place all interpolated child 
!***  data onto the full child grid for writing out.
!-----------------------------------------------------------------------
!
      NN=0
!
      DO J=1,JM_CHILD
      DO I=1,IM_CHILD
!
        CALL CONVERT_IJ_TO_LATLON(I,J                                   &  !<-- A point on the child grid
                                 ,IM_CHILD,JM_CHILD                     &  !<-- Dimensions of child grid
                                 ,TPH0D_CHILD,TLM0D_CHILD               &  !<-- Parent lat/lon (deg) of child grid central point
                                 ,DPHD_CHILD,DLMD_CHILD                 &  !<-- Angular grid increments (deg) on child grid
                                 ,CHILD_LATD_ON_PARENT                  &  !<-- Parent latitude of child point
                                 ,CHILD_LOND_ON_PARENT )                   !<-- Parent longitude of child point
!
        REAL_I_PARENT=(CHILD_LOND_ON_PARENT-WBD_PARENT)*R_DLMD+ISTART      !<-- REAL I index of child point on parent grid
        REAL_J_PARENT=(CHILD_LATD_ON_PARENT-SBD_PARENT)*R_DPHD+JSTART      !<-- REAL J index of child point on parent grid
!
!-----------------------------------------------------------------------
!
        IF(REAL(ITS)<=REAL_I_PARENT.AND.REAL(I_END)>REAL_I_PARENT.AND.  &  !<-- Is child gridpoint on this parent task?
           REAL(JTS)<=REAL_J_PARENT.AND.REAL(J_END)>REAL_J_PARENT)THEN     !<--
!
          NUM_CHILD_POINTS=NUM_CHILD_POINTS+1                              !<-- Add up number of child points on this parent task
!
          CHILD_POINT_INDICES(2*NUM_CHILD_POINTS-1)=I                      !<-- Save I index of this child
          CHILD_POINT_INDICES(2*NUM_CHILD_POINTS  )=J                      !<-- Save J index of this child point
!
!-----------------------------------------------------------------------
!***  Compute the distance from the child point location to each of
!***  the four surrounding parent points and generate the bilinear
!***  interpolation weights.
!***  The indices 1-->4 indicate the parent points to the SW, SE,
!***  NW, and NE in that order.
!-----------------------------------------------------------------------
!
          I_PARENT_SW=INT(REAL_I_PARENT)
          J_PARENT_SW=INT(REAL_J_PARENT)
          RLATD(1)=(J_PARENT_SW-ROW_0)*DPHD_PARENT                         !<-- Parent latitude (deg) of parent point SW of child point
          RLOND(1)=(I_PARENT_SW-COL_0)*DLMD_PARENT                         !<-- Parent longitude (deg) of parent point SW of child point
!
          I_PARENT_SE=I_PARENT_SW+1         
          J_PARENT_SE=J_PARENT_SW         
          RLATD(2)=RLATD(1)                                                !<-- SE and SW on same line of parent latitude
          RLOND(2)=RLOND(1)+DLMD_PARENT                                    !<-- SE is one gridpoint east of SW parent point
!
          I_PARENT_NW=I_PARENT_SW           
          J_PARENT_NW=J_PARENT_SW+1       
          RLATD(3)=RLATD(1)+DPHD_PARENT                                    !<-- NW is one gridpoint north of SW parent point
          RLOND(3)=RLOND(1)                                                !<-- NW and SW on same line of parent longitude
!
          I_PARENT_NE=I_PARENT_SE           
          J_PARENT_NE=J_PARENT_NW         
          RLATD(4)=RLATD(3)                                                !<-- NE and NW on same line of parent latitude
          RLOND(4)=RLOND(2)                                                !<-- NE and SE on same line of parent longitude
!
          SUM=0.
!
          DO N=1,4                                                         !<-- Loop over SW, SE, NW, and NE parent points
!
            CALL DISTANCE_ON_SPHERE(CHILD_LATD_ON_PARENT               &   !<-- Parent latitiude (deg) of child gridpoint
                                   ,CHILD_LOND_ON_PARENT               &   !<-- Parent latitiude (deg) of child gridpoint
                                   ,RLATD(N)                           &   !<-- Latitude (deg) of surrounding parent point N
                                   ,RLOND(N)                           &   !<-- Longitude (deg) of surrounding parent point N
                                   ,DIST )                                 !<-- Distance (radians) from child point to parent point N
!
            WGT(N)=1./DIST
            SUM=SUM+WGT(N)
!
          ENDDO
!
          SUM_RECIP=1./SUM
!
!-----------------------------------------------------------------------
!***  The bilinear interpolation weights of the four parent points
!***  surrounding the child point.
!-----------------------------------------------------------------------
!
          WEIGHT_SW=WGT(1)*SUM_RECIP
          WEIGHT_SE=WGT(2)*SUM_RECIP
          WEIGHT_NW=WGT(3)*SUM_RECIP
          WEIGHT_NE=WGT(4)*SUM_RECIP
!
          IF(ABS(PARENT_ARRAY(I_PARENT_SW,J_PARENT_SW,1))<1.E-12)WEIGHT_SW=0.
          IF(ABS(PARENT_ARRAY(I_PARENT_SE,J_PARENT_SE,1))<1.E-12)WEIGHT_SE=0.
          IF(ABS(PARENT_ARRAY(I_PARENT_NW,J_PARENT_NW,1))<1.E-12)WEIGHT_NW=0.
          IF(ABS(PARENT_ARRAY(I_PARENT_NE,J_PARENT_NE,1))<1.E-12)WEIGHT_NE=0.
          WEIGHT_SUM=WEIGHT_SW+WEIGHT_SE+WEIGHT_NW+WEIGHT_NE
          WEIGHT_SUM_RECIP=1./WEIGHT_SUM
!
          DO L=1,NLEV
            NN=NN+1
!
            CHILD_STRING(NN)=WEIGHT_SW*PARENT_ARRAY(I_PARENT_SW,J_PARENT_SW,L)  & !<-- Value at points on child's grid
                            +WEIGHT_SE*PARENT_ARRAY(I_PARENT_SE,J_PARENT_SE,L)  &
                            +WEIGHT_NW*PARENT_ARRAY(I_PARENT_NW,J_PARENT_NW,L)  &
                            +WEIGHT_NE*PARENT_ARRAY(I_PARENT_NE,J_PARENT_NE,L)
!
            IF(WEIGHT_SUM<0.99.AND.WEIGHT_SUM>0.01)THEN
              CHILD_STRING(NN)=CHILD_STRING(NN)*WEIGHT_SUM_RECIP           !<-- Normalize if some weights are zero (e.g., coastal land Temp)
            ENDIF
!
            IF(VBL_NAME=='SeaMask')THEN
              IF(CHILD_STRING(NN)>=0.5)CHILD_STRING(NN)=1.0
              IF(CHILD_STRING(NN)< 0.5)CHILD_STRING(NN)=0.0
            ENDIF
!
          ENDDO
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Each parent task that contains a child grid point has now done
!***  the horizontal interpolation of the parent variable to the child.
!***  Now parent task 0 receives all the interpolated data from the
!***  other parent tasks.
!-----------------------------------------------------------------------
!
      data_fill: IF(MYPE==0)THEN
!
!-----------------------------------------------------------------------
!
        remote_tasks: DO IPE=1,NUM_PES_PARENT-1
!        
          CALL MPI_RECV(NUM_POINTS_REMOTE                               &  !<-- # of child points on parent task IPE
                       ,1                                               &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(NUM_POINTS_REMOTE==0)CYCLE remote_tasks
!
          NUM_IJ  =2*NUM_POINTS_REMOTE
          NUM_DATA=NUM_POINTS_REMOTE*NLEV
!
          ALLOCATE(DATA_REMOTE(1:NUM_DATA))
          ALLOCATE(IJ_REMOTE  (1:NUM_IJ  ))
!
          CALL MPI_RECV(DATA_REMOTE                                     &  !<-- Interpolated data on child grid from parent task IPE
                       ,NUM_DATA                                        &  !<-- Total words received
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(IJ_REMOTE                                       &  !<-- Interpolated data on child grid from parent task IPE
                       ,NUM_IJ                                          &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
!-----------------------------------------------------------------------
!***  Parent task 0 fills in section of child data from remote 
!***  parent task IPE.
!-----------------------------------------------------------------------
!
          KOUNT=0
!
          DO L=1,NLEV
            DO NIJ=1,NUM_POINTS_REMOTE
              KOUNT=KOUNT+1
              I=IJ_REMOTE(2*NIJ-1)
              J=IJ_REMOTE(2*NIJ  )
              CHILD_ARRAY(I,J,L)=DATA_REMOTE(KOUNT)
            ENDDO
          ENDDO
!
!-----------------------------------------------------------------------
!
          DEALLOCATE(DATA_REMOTE)
          DEALLOCATE(IJ_REMOTE)
!
!-----------------------------------------------------------------------
!
        ENDDO remote_tasks
!
!-----------------------------------------------------------------------
!***  Finally parent task 0 fills in its own section of the child array.
!-----------------------------------------------------------------------
!
        IF(NUM_CHILD_POINTS>0)THEN
!
          KOUNT=0          
          DO L=1,NLEV
            DO N=1,NUM_CHILD_POINTS
              KOUNT=KOUNT+1
              I=CHILD_POINT_INDICES(2*N-1)
              J=CHILD_POINT_INDICES(2*N  )
              CHILD_ARRAY(I,J,L)=CHILD_STRING(KOUNT)
            ENDDO
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ELSE data_fill
!
!-----------------------------------------------------------------------
!***  Remote parent tasks send their sections of interpolated child
!***  data to parent task 0.
!-----------------------------------------------------------------------
!
        CALL MPI_SEND(NUM_CHILD_POINTS                                  &  !<-- # of child points on this parent task
                     ,1                                                 &  
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Send to parent task 0
                     ,MYPE                                              &  !<-- MPI tag
                     ,COMM_MY_DOMAIN                                    &  !<-- The MPI communicator
                     ,IERR )
!
        IF(NUM_CHILD_POINTS>0)THEN
          NUM_IJ  =2*NUM_CHILD_POINTS
          NUM_DATA=NUM_CHILD_POINTS*NLEV
!
          CALL MPI_SEND(CHILD_STRING                                    &  !<-- Interpolated data on child grid for this parent task
                       ,NUM_DATA                                        &  !<-- Total words sent
                       ,MPI_REAL                                        &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
          CALL MPI_SEND(CHILD_POINT_INDICES                             &  !<-- Indices of child points for this parent task
                       ,NUM_IJ                                          &  !<-- Total words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF data_fill
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(CHILD_POINT_INDICES)
      DEALLOCATE(CHILD_STRING)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_TO_CHILD_FILL_GENERAL
!!!   END SUBROUTINE PARENT_TO_CHILD_FILL
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_TO_CHILD_IFILL_GENERAL(PARENT_ARRAY             &
!!!   SUBROUTINE PARENT_TO_CHILD_IFILL        (PARENT_ARRAY             &
                                              ,VBL_NAME                 &
                                              ,CHILD_ARRAY)
!
!-----------------------------------------------------------------------
!***  Parent tasks interpolate their data to the locations of their
!***  children's gridpoints.  The child grids are unique rotated
!***  lat/lon grids with their own centers.  The southwest H point
!***  of the child grid lies directly on an H point of the parent. 
!
!***  Only parent tasks participate in this work.
!-----------------------------------------------------------------------
!
!---------------
!***  Arguments
!---------------
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: PARENT_ARRAY        !<-- The parent array that will initialize the child array
!
      CHARACTER(*),INTENT(IN) :: VBL_NAME                                  !<-- The variable's name 
!
      INTEGER,DIMENSION(1:IM_CHILD,1:JM_CHILD),INTENT(OUT) :: CHILD_ARRAY  !<-- Data from parent tasks interpolated to child grid
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(KIND=KINT) :: I,I_END,ISTART,J,J_END,JSTART               &
                           ,KOUNT,NIJ,NN,NTOT                           &
                           ,NUM_DATA                                    &
                           ,NUM_CHILD_POINTS                            &
                           ,NUM_IJ                                      &
                           ,NUM_POINTS_REMOTE
!
      INTEGER(KIND=KINT) :: I_PARENT_SW,I_PARENT_SE                     &
                           ,I_PARENT_NW,I_PARENT_NE                     &
                           ,J_PARENT_SW,J_PARENT_SE                     &
                           ,J_PARENT_NW,J_PARENT_NE
!
      INTEGER(KIND=KINT) :: IERR,IPE,ISTAT
!
      INTEGER(KIND=KINT),DIMENSION(:),ALLOCATABLE :: CHILD_POINT_INDICES &
                                                    ,CHILD_STRING        &
                                                    ,DATA_REMOTE         &
                                                    ,IJ_REMOTE  
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      REAL(KIND=KFPT) :: CHILD_LATD_ON_PARENT                           &
                        ,CHILD_LOND_ON_PARENT                           &
                        ,DIST                                           &
                        ,R_DLMD,R_DPHD                                  &
                        ,REAL_I_PARENT                                  &
                        ,REAL_J_PARENT                                  &
                        ,RLATD_SW,RLOND_SW                              &
                        ,RLATD_SE,RLOND_SE                              &
                        ,RLATD_NW,RLOND_NW                              &
                        ,RLATD_NE,RLOND_NE                              &
                        ,SUM,SUM_RECIP                                  &
                        ,WEIGHT_SW,WEIGHT_SE                            &
                        ,WEIGHT_NW,WEIGHT_NE                            &
                        ,WEIGHT_MAX
!
      REAL(KIND=KFPT),DIMENSION(4) :: RLATD,RLOND,WGT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      R_DPHD=1./DPHD_PARENT
      R_DLMD=1./DLMD_PARENT
!
      NUM_CHILD_POINTS=0
!
      ISTART=1
      JSTART=1
      IF(GLOBAL)THEN
        ISTART=2
        JSTART=2
      ENDIF
!
!-----------------------------------------------------------------------
!***  Each parent task is responsible for searching the parent domain 
!***  extending from ITS and JTS to ITE+1 and JTE+1.  Those latter +1's
!***  are needed in order to reach the next gridpoint in each direction.
!***  We cannot go outside the full domain of course.
!-----------------------------------------------------------------------
!
      I_END=MIN(ITE+1,IDE)
      J_END=MIN(JTE+1,JDE)
      if(vbl_name=='ISLTYP')then
        write(0,*)' computed j_end=',j_end,' jte+1=',jte+1,' jde=',jde
      endif
!
!-----------------------------------------------------------------------
!
      NTOT=2*IM_CHILD*JM_CHILD
      ALLOCATE(CHILD_POINT_INDICES(1:NTOT),stat=ISTAT)
      IF(ISTAT/=0)WRITE(0,*)' Failed to allocate CHILD_POINT_INDICES in PARENT_TO_CHILD_IFILL_GENERAL'
!
      NTOT=IM_CHILD*JM_CHILD
      ALLOCATE(CHILD_STRING(1:NTOT),stat=ISTAT)
      IF(ISTAT/=0)WRITE(0,*)' Failed to allocate CHILD_STRING in PARENT_TO_CHILD_IFILL_GENERAL'
!
!-----------------------------------------------------------------------
!***  Compute the parent's lat/lon of each child gridpoint in order to
!***  determine if that gridpoint lies on a given parent task. 
!***  Save the I,J of each gridpoint found since ultimately parent
!***  task 0 will need that to properly place all interpolated child 
!***  data onto the full child grid for writing out.
!-----------------------------------------------------------------------
!
      NN=0
!
      DO J=1,JM_CHILD
      DO I=1,IM_CHILD
!
        CALL CONVERT_IJ_TO_LATLON(I,J                                   &  !<-- A point on the child grid
                                 ,IM_CHILD,JM_CHILD                     &  !<-- Dimensions of child grid
                                 ,TPH0D_CHILD,TLM0D_CHILD               &  !<-- Parent lat/lon (deg) of child grid central point
                                 ,DPHD_CHILD,DLMD_CHILD                 &  !<-- Angular grid increments (deg) on child grid
                                 ,CHILD_LATD_ON_PARENT                  &  !<-- Parent latitude of child point
                                 ,CHILD_LOND_ON_PARENT )                   !<-- Parent longitude of child point
!
        REAL_I_PARENT=(CHILD_LOND_ON_PARENT-WBD_PARENT)*R_DLMD+ISTART      !<-- REAL I index of child point on parent grid
        REAL_J_PARENT=(CHILD_LATD_ON_PARENT-SBD_PARENT)*R_DPHD+JSTART      !<-- REAL J index of child point on parent grid
      if(i==01.and.j==01.and.vbl_name=='ISLTYP')then
      write(0,*)' ISLTYP i=',i,' j=',j,' CHILD_LATD_ON_PARENT=',CHILD_LATD_ON_PARENT,' CHILD_LOND_ON_PARENT=',CHILD_LOND_ON_PARENT &
               ,' REAL_I_PARENT=',REAL_I_PARENT,' REAL_J_PARENT=',REAL_J_PARENT
      write(0,*)' its=',its,' real_i_parent=',real_i_parent,' i_end=',i_end &
               ,' jts=',jts,' real_j_parent=',real_j_parent,' j_end=',j_end
      endif
!
!-----------------------------------------------------------------------
!
        IF(REAL(ITS)<=REAL_I_PARENT.AND.REAL(I_END)>REAL_I_PARENT.AND.  &  !<-- Is child gridpoint on this parent task?
           REAL(JTS)<=REAL_J_PARENT.AND.REAL(J_END)>REAL_J_PARENT)THEN     !<--
!
          NUM_CHILD_POINTS=NUM_CHILD_POINTS+1                              !<-- Add up number of child points on this parent task
!
          CHILD_POINT_INDICES(2*NUM_CHILD_POINTS-1)=I                      !<-- Save I index of this child
          CHILD_POINT_INDICES(2*NUM_CHILD_POINTS  )=J                      !<-- Save J index of this child point
      if(i==01.and.j==01.and.vbl_name=='ISLTYP')then
      write(0,*)' found child point ISLTYP i=',i,' j=',j,' num_child_points=',num_child_points
      write(0,*)' CHILD_LATD_ON_PARENT=',CHILD_LATD_ON_PARENT,' SBD_PARENT=',SBD_PARENT,' JSTART=',JSTART,' REAL_J_PARENT=',REAL_J_PARENT
      write(0,*)' CHILD_LOND_ON_PARENT=',CHILD_LOND_ON_PARENT,' WBD_PARENT=',WBD_PARENT,' ISTART=',ISTART,' REAL_I_PARENT=',REAL_I_PARENT
      endif
!
!-----------------------------------------------------------------------
!***  Compute the distance from the child point location to each of
!***  the four surrounding parent points and generate the bilinear
!***  interpolation weights.
!***  The indices 1-->4 indicate the parent points to the SW, SE,
!***  NW, and NE in that order.
!-----------------------------------------------------------------------
!
          I_PARENT_SW=INT(REAL_I_PARENT)
          J_PARENT_SW=INT(REAL_J_PARENT)
          RLATD(1)=(J_PARENT_SW-ROW_0)*DPHD_PARENT                         !<-- Parent latitude (deg) of parent point SW of child point
          RLOND(1)=(I_PARENT_SW-COL_0)*DLMD_PARENT                         !<-- Parent longitude (deg) of parent point SW of child point
!
          I_PARENT_SE=I_PARENT_SW+1         
          J_PARENT_SE=J_PARENT_SW         
          RLATD(2)=RLATD(1)                                                !<-- SE and SW on same line of parent latitude
          RLOND(2)=RLOND(1)+DLMD_PARENT                                    !<-- SE is one gridpoint east of SW parent point
!
          I_PARENT_NW=I_PARENT_SW           
          J_PARENT_NW=J_PARENT_SW+1       
          RLATD(3)=RLATD(1)+DPHD_PARENT                                    !<-- NW is one gridpoint north of SW parent point
          RLOND(3)=RLOND(1)                                                !<-- NW and SW on same line of parent longitude
!
          I_PARENT_NE=I_PARENT_SE           
          J_PARENT_NE=J_PARENT_NW         
          RLATD(4)=RLATD(3)                                                !<-- NE and NW on same line of parent latitude
          RLOND(4)=RLOND(2)                                                !<-- NE and SE on same line of parent longitude
!
          SUM=0.
!
          DO N=1,4                                                         !<-- Loop over SW, SE, NW, and NE parent points
!
            CALL DISTANCE_ON_SPHERE(CHILD_LATD_ON_PARENT               &   !<-- Parent latitiude (deg) of child gridpoint
                                   ,CHILD_LOND_ON_PARENT               &   !<-- Parent latitiude (deg) of child gridpoint
                                   ,RLATD(N)                           &   !<-- Latitude (deg) of surrounding parent point N
                                   ,RLOND(N)                           &   !<-- Longitude (deg) of surrounding parent point N
                                   ,DIST )                                 !<-- Distance (radians) from child point to parent point N
!
            WGT(N)=1./DIST
            SUM=SUM+WGT(N)
!
          ENDDO
!
          SUM_RECIP=1./SUM
!
!-----------------------------------------------------------------------
!***  The bilinear interpolation weights of the four parent points
!***  surrounding the child point.
!-----------------------------------------------------------------------
!
          WEIGHT_SW=WGT(1)*SUM_RECIP
          WEIGHT_SE=WGT(2)*SUM_RECIP
          WEIGHT_NW=WGT(3)*SUM_RECIP
          WEIGHT_NE=WGT(4)*SUM_RECIP
          WEIGHT_MAX=MAX(WEIGHT_SW,WEIGHT_SE,WEIGHT_NW,WEIGHT_NE)
!
!-----------------------------------------------------------------------
!***  Using the bilinear interpolation weights, assign the value of
!***  the nearest parent point to the child point.
!-----------------------------------------------------------------------
!
          NN=NN+1
!
          IF(WEIGHT_SW==WEIGHT_MAX)THEN
            CHILD_STRING(NN)=PARENT_ARRAY(I_PARENT_SW,J_PARENT_SW)
!     write(0,*)' SW parent=',PARENT_ARRAY(I_PARENT_SW,J_PARENT_SW),' nn=',nn,' i=',i,' j=',j
          ELSEIF(WEIGHT_SE==WEIGHT_MAX)THEN
            CHILD_STRING(NN)=PARENT_ARRAY(I_PARENT_SE,J_PARENT_SE)
!     write(0,*)' SE parent=',PARENT_ARRAY(I_PARENT_SE,J_PARENT_SE),' nn=',nn,' i=',i,' j=',j
          ELSEIF(WEIGHT_NW==WEIGHT_MAX)THEN
            CHILD_STRING(NN)=PARENT_ARRAY(I_PARENT_NW,J_PARENT_NW)
!     write(0,*)' NW parent=',PARENT_ARRAY(I_PARENT_NW,J_PARENT_NW),' nn=',nn,' i=',i,' j=',j
          ELSEIF(WEIGHT_NE==WEIGHT_MAX)THEN
            CHILD_STRING(NN)=PARENT_ARRAY(I_PARENT_NE,J_PARENT_NE)
!     write(0,*)' NE parent=',PARENT_ARRAY(I_PARENT_NE,J_PARENT_NE),' nn=',nn,' i=',i,' j=',j
          ENDIF
      if(i==01.and.j==01.and.vbl_name=='ISLTYP')then
        write(0,*)' parent interp value to ISLTYP is ',CHILD_STRING(NN),' nn=',nn
      endif
!     if(vbl_name=='ISLTYP'.and.child_string(nn)<1.and.sea_mask(i,j)<0.5)then
!       write(0,*)' Parent creating bad value of ISLTYP=',CHILD_STRING(NN),' nn=',nn,' i=',i,' j=',j &
!                ,' SeaMask=',sea_mask(i,j)
!     endif
        
!
!-----------------------------------------------------------------------
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  Each parent task that contains a child grid point has now done
!***  the horizontal interpolation of the parent variable to the child.
!***  Now parent task 0 receives all the interpolated data from the
!***  other parent tasks.
!-----------------------------------------------------------------------
!
      data_fill: IF(MYPE==0)THEN
!
!-----------------------------------------------------------------------
!
        remote_tasks: DO IPE=1,NUM_PES_PARENT-1
!        
          CALL MPI_RECV(NUM_POINTS_REMOTE                               &  !<-- # of child points on parent task IPE
                       ,1                                               &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          IF(NUM_POINTS_REMOTE==0)CYCLE remote_tasks
!
          NUM_IJ  =2*NUM_POINTS_REMOTE
          NUM_DATA=NUM_POINTS_REMOTE
!
          ALLOCATE(DATA_REMOTE(1:NUM_DATA))
          ALLOCATE(IJ_REMOTE  (1:NUM_IJ  ))
!
          CALL MPI_RECV(DATA_REMOTE                                     &  !<-- Interpolated data on child grid from parent task IPE
                       ,NUM_DATA                                        &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
          CALL MPI_RECV(IJ_REMOTE                                       &  !<-- Interpolated data on child grid from parent task IPE
                       ,NUM_IJ                                          &  !<-- Total words received
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,IPE                                             &  !<-- Receive from this parent task
                       ,IPE                                             &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,JSTAT                                           &  !<-- MPI status object
                       ,IERR )
!
!-----------------------------------------------------------------------
!***  Parent task 0 fills in section of child data from remote 
!***  parent task IPE.
!-----------------------------------------------------------------------
!
          KOUNT=0
!
          DO NIJ=1,NUM_POINTS_REMOTE
            KOUNT=KOUNT+1
            I=IJ_REMOTE(2*NIJ-1)
            J=IJ_REMOTE(2*NIJ  )
            CHILD_ARRAY(I,J)=DATA_REMOTE(KOUNT)
!     write(0,*)' new child i=',i,' j=',j,' kount=',kount,' data=',data_remote(kount)
      if(vbl_name=='ISLTYP'.and.i==01.and.j==01)then
        write(0,*)' parent fills ISLTYP with ',DATA_REMOTE(KOUNT),' kount=',kount
      endif
          ENDDO
!
!-----------------------------------------------------------------------
!
          DEALLOCATE(DATA_REMOTE)
          DEALLOCATE(IJ_REMOTE)
!
!-----------------------------------------------------------------------
!
        ENDDO remote_tasks
!
!-----------------------------------------------------------------------
!***  Finally parent task 0 fills in its own section of the child array.
!-----------------------------------------------------------------------
!
        IF(NUM_CHILD_POINTS>0)THEN
!
          KOUNT=0          
          DO N=1,NUM_CHILD_POINTS
            KOUNT=KOUNT+1
            I=CHILD_POINT_INDICES(2*N-1)
            J=CHILD_POINT_INDICES(2*N  )
            CHILD_ARRAY(I,J)=CHILD_STRING(KOUNT)
          ENDDO
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ELSE data_fill
!
!-----------------------------------------------------------------------
!***  Remote parent tasks send their sections of interpolated child
!***  data to parent task 0.
!-----------------------------------------------------------------------
!
        CALL MPI_SEND(NUM_CHILD_POINTS                                  &  !<-- # of child points on this parent task
                     ,1                                                 &  
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Send to parent task 0
                     ,MYPE                                              &  !<-- MPI tag
                     ,COMM_MY_DOMAIN                                    &  !<-- The MPI communicator
                     ,IERR )
!
        IF(NUM_CHILD_POINTS>0)THEN
          NUM_IJ  =2*NUM_CHILD_POINTS
          NUM_DATA=NUM_CHILD_POINTS
!
          CALL MPI_SEND(CHILD_STRING                                    &  !<-- Interpolated data on child grid for this parent task
                       ,NUM_DATA                                        &  !<-- Total words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
          CALL MPI_SEND(CHILD_POINT_INDICES                             &  !<-- Indices of child points for this parent task
                       ,NUM_IJ                                          &  !<-- Total words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,0                                               &  !<-- Send to parent task 0
                       ,MYPE                                            &  !<-- MPI tag
                       ,COMM_MY_DOMAIN                                  &  !<-- The MPI communicator
                       ,IERR )
!
        ENDIF
!
!-----------------------------------------------------------------------
!
      ENDIF data_fill
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(CHILD_POINT_INDICES)
      DEALLOCATE(CHILD_STRING)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_TO_CHILD_IFILL_GENERAL
!!!   END SUBROUTINE PARENT_TO_CHILD_IFILL
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_TO_CHILD_INIT_NMM
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE CONVERT_IJ_TO_LATLON  (I_INDEX                         &
                                       ,J_INDEX                         &
                                       ,IM                              &
                                       ,JM                              &
                                       ,TPH0D                           &
                                       ,TLM0D                           &
                                       ,DPHD                            &
                                       ,DLMD                            &
                                       ,RLATD                           &
                                       ,RLOND)
!
!-----------------------------------------------------------------------
!
!***  GIVEN THE (I,J) OF MASS POINTS ON AN ARAKAWA B-GRID,
!***  COMPUTE THE LATITUDES AND LONGITUDES BEFORE ROTATION.
!
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),INTENT(IN) :: I_INDEX                          &  !<-- I value on the grid
                                      ,J_INDEX                          &  !<-- J value on the grid
                                      ,IM                               &  !<-- Full I dimension
                                      ,JM                                  !<-- Full J dimension
!
      REAL(KIND=KFPT),INTENT(IN) :: DPHD                                &  !<-- Latitude grid increment (degrees)
                                   ,DLMD                                &  !<-- Longitude grid increment (degrees)
                                   ,TPH0D                               &  ! Central latitude (deg, positive north), unrotated system
                                   ,TLM0D                                  ! Central longitude (deg, positive east), unrotated system
!
      REAL,INTENT(OUT) :: RLATD                                         &  !<-- Latitude (deg, positive north) of point, unrotated system
                         ,RLOND                                            !<-- Longitude (deg, positive east) of point, unrotated system
!
!-----------------------------------------------------------------------
!
      INTEGER :: I,IEND,ISTART,J,JEND,JSTART
!
      REAL(KIND=DOUBLE) :: ARG1,ARG2,COL_MID,D2R,FCTR,GLATR,GLATD,GLOND &
                          ,HALF,ONE,PI,R2D,ROW_MID,TLATD,TLOND          &
                          ,TLATR,TLONR,TLM0,TPH0
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***
!***  CONVERT FROM TRANSFORMED GRID LOCATION (I,J) 
!***  TO GEOGRAPHIC COORDINATES (DEGREES).
!***
!-----------------------------------------------------------------------
!
      ONE=1.0
      HALF=1./2.
      PI=DACOS(-ONE)
      D2R=PI/180.
      R2D=1./D2R
      TPH0=TPH0D*D2R
      TLM0=TLM0D*D2R
!
      ROW_MID=(JM+ONE)*HALF
      COL_MID=(IM+ONE)*HALF
!
!-----------------------------------------------------------------------
!
      J=J_INDEX
      I=I_INDEX
!
!-----------------------------------------------------------------------
!***  FIND THE ROTATED LATITUDE (POSITIVE NORTH) AND 
!***  LONGITUDE (POSITIVE EAST).
!-----------------------------------------------------------------------
!
      TLATD=(J-ROW_MID)*DPHD
      TLOND=(I-COL_MID)*DLMD
!     write(0,*)' row_mid=',row_mid,' dphd=',dphd,' tlatd=',tlatd
!     write(0,*)' col_mid=',col_mid,' dlmd=',dlmd,' tlond=',tlond
!
!     WRITE(0,50)I,J,TLATD,TLOND
   50 FORMAT(' I=',I4,' J=',I4,' ROTATED LATITUDE IS',F8.3              &
                              ,4X,'LONGITUDE IS',F8.3)
!
!-----------------------------------------------------------------------
!***  NOW CONVERT TO GEOGRAPHIC LATITUDE (POSITIVE NORTH) AND
!***  LONGITUDE (POSITIVE WEST) IN DEGREES.
!-----------------------------------------------------------------------
!
      TLATR=TLATD*D2R
      TLONR=TLOND*D2R
      ARG1=SIN(TLATR)*COS(TPH0)+COS(TLATR)*SIN(TPH0)*COS(TLONR)
      GLATR=ASIN(ARG1)
!
      GLATD=GLATR*R2D
!
      ARG2=DCOS(TLATR)*DCOS(TLONR)/(DCOS(GLATR)*DCOS(TPH0))-            &
           DTAN(GLATR)*DTAN(TPH0)
      IF(ABS(ARG2)>1.)ARG2=ABS(ARG2)/ARG2
      FCTR=1.
      IF(TLOND>0.)FCTR=1.
      IF(TLOND>180.)FCTR=-1.
!
      GLOND=-TLM0D+FCTR*DACOS(ARG2)*R2D
!
!     WRITE(6,100)I,J,GLATD,GLOND
  100 FORMAT(' I=',I3,' J=',I3                                          &
            ,'  PARENT LATITUDE=',F9.5,'  LONGITUDE=',F10.5)
!-----------------------------------------------------------------------
!
      RLATD=GLATD
      RLOND=-GLOND
!     write(0,*)' exit CONVERT_ rlatd=',rlatd,' rlond=',rlond
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CONVERT_IJ_TO_LATLON
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE DISTANCE_ON_SPHERE(RLAT_1,RLON_1                      &
                                   ,RLAT_2,RLON_2                      &
                                   ,DISTANCE )                  
!
!-----------------------------------------------------------------------
!***  COMPUTE THE GREAT CIRCLE DISTANCE BETWEEN TWO POINTS ON THE EARTH.
!-----------------------------------------------------------------------
!
!---------------
!***  Arguments
!---------------
!
      REAL(KIND=KFPT),INTENT(IN) :: RLAT_1,RLON_1                       &  !<-- Lat/lon (deg, +east) of point 1
                        ,RLAT_2,RLON_2                                     !<-- Lat/lon (deg, +east) of point 2
!
      REAL(KIND=KFPT),INTENT(OUT) :: DISTANCE                              !<-- Distance (radians) between points 1 and 2
!
!-----------------------------------------------------------------------
!
!--------------------
!*** Local Variables
!--------------------
!
      REAL(KIND=KDBL) :: ALPHA,ARG,BETA,CROSS,DLON,DTR                  &
                        ,PHI1,PHI2,PI,PI_H
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      PI_H=ACOS(0.)
      PI=2.*PI_H
      DTR=PI/180.
!
!-----------------------------------------------------------------------
!
      PHI1=RLAT_1*DTR
      PHI2=RLAT_2*DTR
      DLON=(RLON_2-RLON_1)*DTR
!
      CROSS=ACOS(COS(DLON)*COS(PHI2))
      ARG=TAN(PHI2)/SIN(DLON)
      ALPHA=ATAN(ARG)
      IF(DLON<0.)ALPHA=-ALPHA
      BETA=PI_H-ALPHA
!
      DISTANCE=ACOS(COS(PHI1)*COS(PHI2)*COS(DLON)                       &
                   +SIN(PHI1)*SIN(CROSS)*COS(BETA))
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE DISTANCE_ON_SPHERE
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
      SUBROUTINE CENTER_NEST(SBD_DOMAIN                                 &
                            ,WBD_DOMAIN                                 &
                            ,SW_CORNER_LATD                             &
                            ,SW_CORNER_LOND                             &
                            ,TPH0D_DOMAIN                               &
                            ,TLM0D_DOMAIN )
!-----------------------------------------------------------------------
!***  GIVEN THE SOUTHERN AND WESTERN BOUNDARIES OF A ROTATED LAT/LON
!***  GRID AS WELL AS THE COORDINATES OF THE SOUTHWEST CORNER POINT,
!***  FIND THE CORRDINATES OF THE GRID'S CENTRAL POINT WITH RESPECT
!***  TO THE GRID UPON WHICH THE ROTATED GRID LIES.
!-----------------------------------------------------------------------
!
!---------------
!***  Arguments
!---------------
!
      REAL(KIND=KFPT),INTENT(IN) :: SBD_DOMAIN                          &  !<-- Latitude (deg) of domain's southern boundary
                                   ,WBD_DOMAIN                          &  !<-- Longitude (deg, +east) of domain's western boundary
                                   ,SW_CORNER_LATD                      &  !<-- Latitude (deg) of domain's southwest corner point
                                   ,SW_CORNER_LOND                         !<-- Longitude (deg, +east) of domain's southwest corner point
!
      REAL(KIND=KFPT),INTENT(OUT) :: TPH0D_DOMAIN                       &  !<-- Latitude (deg) of domain's center
                                    ,TLM0D_DOMAIN                          !<-- Longitude (deg) of domain's center
!
!-----------------------------------------------------------------------
!
!---------------------
!***  Local Variables
!---------------------
!
      REAL :: ALPHA,BETA,CENTRAL_LAT,CENTRAL_LON                        &
             ,DEG_RAD,DELTA,GAMMA                                       &
             ,PI_2,SB_R,SIDE1,SIDE2,SIDE3,SIDE4,SIDE5                   &
             ,SW_LAT,SW_LON,WB_R
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      PI_2=ACOS(0.)
      DEG_RAD=PI_2/90.
!
!-----------------------------------------------------------------------
!***  SOUTHERN AND WESTERN BOUNDARIES OF THE ROTATED DOMAIN IN RADIANS.
!-----------------------------------------------------------------------
!
      SB_R=-SBD_DOMAIN*DEG_RAD
      WB_R=-WBD_DOMAIN*DEG_RAD
!
!-----------------------------------------------------------------------
!***  SOUTHWEST CORNER OF THE DOMAIN IN RADIANS.
!-----------------------------------------------------------------------
!
      SW_LAT=SW_CORNER_LATD*DEG_RAD
      SW_LON=SW_CORNER_LOND*DEG_RAD
      write(0,*)' CENTER_NEST SBD_DOMAIN=',SBD_DOMAIN,' WBD_DOMAIN=',WBD_DOMAIN &
               ,' SW_CORNER_LATD=',SW_CORNER_LATD,' SW_CORNER_LOND=',SW_CORNER_LOND
!
!-----------------------------------------------------------------------
!***  SIDE1 IS THE ARC FROM THE SOUTHWEST CORNER TO THE CENTER 
!***  OF THE DOMAIN.
!-----------------------------------------------------------------------
!
      SIDE1=ACOS(COS(SB_R)*COS(WB_R))
!
!-----------------------------------------------------------------------
!***  ALPHA IS THE ANGLE BETWEEN SIDE1 AND THE DOMAIN'S EQUATOR WEST OF
!***  THE CENTRAL POINT.
!-----------------------------------------------------------------------
!
      ALPHA=ATAN(TAN(SB_R)/SIN(WB_R))
!
!-----------------------------------------------------------------------
!***  BETA IS THE ANGLE BETWEEN SIDE1 AND THE DOMAIN'S PRIME MERIDIAN
!***  SOUTH OF THE CENTRAL POINT.
!-----------------------------------------------------------------------
!
      BETA=PI_2-ALPHA
!
!-----------------------------------------------------------------------
!***  SIDE2 IS THE ARC FROM THE CENTRAL POINT SOUTHWARD ALONG THE 
!***  DOMAIN'S PRIME MERIDIAN TO THE GREAT CIRCLE THAT INTERSECTS 
!***  BOTH THE SW AND SE CORNERS OF THE DOMAIN.
!-----------------------------------------------------------------------
!
      SIDE2=ATAN(COS(BETA)*TAN(SIDE1))
!
!-----------------------------------------------------------------------
!***  SIDE3 IS THE ARC BETWEEN THE DOMAIN'S PRIME MERIDIAN AND THE SW
!***  CORNER ALONG THE GREAT CIRCLE THAT CONNECTS THE DOMAIN'S SW AND
!***  SE CORNERS.
!-----------------------------------------------------------------------
!
      SIDE3=ASIN(SIN(BETA)*SIN(SIDE1))
!
!-----------------------------------------------------------------------
!***  SIDE4 IS THE ARC ALONG THE OUTER GRID'S EQUATOR THAT LIES BETWEEN 
!***  ITS WESTERN INTERSECTION WITH THE ABOVE MENTIONED GREAT CIRCLE 
!***  AND THE OUTER GRID'S MERIDIAN THAT PASSES THROUGH THE DOMAIN'S 
!***  SW CORNER.
!-----------------------------------------------------------------------
!
      SIDE4=ACOS(SIN(SIDE3)/COS(SW_LAT))
!
!-----------------------------------------------------------------------
!***  GAMMA IS THE ANGLE BETWEEN THE OUTER GRID'S EQUATOR AND THE ARC 
!***  THAT CONNECTS THE DOMAIN'S SW CORNER WITH THE POINT WHERE THE 
!***  DOMAIN'S CENTRAL MERIDIAN CROSSES THE OUTER GRID'S EQUATOR.
!-----------------------------------------------------------------------
!
      GAMMA=ATAN(TAN(SW_LAT)/COS(SIDE4))
!
!-----------------------------------------------------------------------
!***  DELTA IS THE ANGLE BETWEEN THE ARC THAT CONNECTS THE DOMAIN'S SW
!***  CORNER WITH THE POINT WHERE THE DOMAIN'S CENTRAL MERIDIAN CROSSES
!***  THE OUTER GRID'S EQUATOR AND THE DOMAIN'S CENTRAL MERIDIAN ITSELF.
!-----------------------------------------------------------------------
!
      DELTA=PI_2-GAMMA
!
!-----------------------------------------------------------------------
!***  SIDE5 IS THE ARC ALONG THE DOMAIN'S CENTRAL MERIDIAN THAT LIES
!***  BETWEEN THE OUTER GRID'S EQUATOR AND THE GREAT CIRCLE THAT PASSES
!***  THROUGH THE SW AND SE CORNERS OF THE DOMAIN.
!-----------------------------------------------------------------------
!
      SIDE5=ASIN(TAN(SIDE3)/TAN(DELTA))
!
!-----------------------------------------------------------------------
!***  THE CENTRAL LATITUDE AND LONGITUDE OF THE DOMAIN IN TERMS OF
!***  THE COORDINATES OF THE OUTER GRID.
!-----------------------------------------------------------------------
!
      CENTRAL_LAT=SIDE2+SIDE5
      CENTRAL_LON=SW_LON+PI_2-SIDE4
!
      TPH0D_DOMAIN=CENTRAL_LAT/DEG_RAD
      TLM0D_DOMAIN=CENTRAL_LON/DEG_RAD
      write(0,*)' CENTER_NEST CENTRAL_LAT=',CENTRAL_LAT,' CENTRAL_LON=',CENTRAL_LON,' DEG_RAD=',DEG_RAD
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CENTER_NEST
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE SET_NEST_GRIDS(DOMAIN_ID_MINE                          &
                               ,TPH0D,TLM0D                             &
                               ,SBD_MINE,WBD_MINE                       &
                               ,DPHD_MINE,DLMD_MINE)
! 
!-----------------------------------------------------------------------
!***  BASIC GRID CHARACTERISTICS FOR NESTS ARE BASED UPON THOSE OF
!***  THE UPPERMOST PARENT GRID.  USE THOSE PARENT VALUES TO COMPUTE
!***  APPROPRIATE ANALOGS FOR THE NESTS.
!***  THIS SUBROUTINE IS RELEVANT ONLY TO GRID-ASSOCIATED NESTS.
!-----------------------------------------------------------------------
!
!---------------
!***  Arguments
!---------------
!
      INTEGER,INTENT(IN) :: DOMAIN_ID_MINE                                 !<-- Domain ID for this nested domain
!
      REAL,INTENT(OUT) :: DPHD_MINE                                     &  !<-- Delta phi of this nested domain (degrees)
                         ,DLMD_MINE                                     &  !<-- Delta lambda of this nested domain (degrees)
                         ,TLM0D                                         &  !<-- Central rotated longitude of all domains (degrees)
                         ,TPH0D                                         &  !<-- Central rotated latitude of all domains (degrees)
                         ,SBD_MINE                                      &  !<-- Southern boundary this nested domain (degrees)
                         ,WBD_MINE                                         !<-- Western boundary this nested domain (degrees)
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: MAX_DOMAINS=99
!
      INTEGER :: IM_1,JM_1                                              &
                ,ID_ANCESTOR,ID_DOMAIN                                  &
                ,IDE_1,JDE_1                                            &
                ,I_BOUND,J_BOUND                                        &
                ,I_PARENT_SW,J_PARENT_SW                                &
                ,I_START_SW,J_START_SW                                  &
                ,N,NUM_ANCESTORS
!
      INTEGER :: RC,RC_SET
!
      INTEGER,DIMENSION(MAX_DOMAINS) :: ID_ANCESTORS=0
      INTEGER,DIMENSION(MAX_DOMAINS) :: PARENT_CHILD_SPACE_RATIO
!
      INTEGER,DIMENSION(2,MAX_DOMAINS) :: SW_CORNER
!
      REAL :: DPHD_1,DLMD_1,TLM_BASE_1,TPH_BASE_1,SBD_1,WBD_1
      REAL :: DPHD_X,DLMD_X,TLM_BASE,TPH_BASE,SBD_X,WBD_X
!
      CHARACTER(2)  :: INT_TO_CHAR
      CHARACTER(6)  :: FMT='(I2.2)'
      CHARACTER(50) :: GLOBAL
      CHARACTER(99) :: CONFIG_FILE_NAME
!
      TYPE(ESMF_Config),DIMENSION(MAX_DOMAINS) :: CF
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_SET=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  FIRST LOAD ALL OF THE DOMAINS' CONFIGURE FILES.
!-----------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*******  See NMM_ATM_INIT where
!*******  CF(N) is made to be
!*******  CF(ID of domain).
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      DO N=1,MAX_DOMAINS
        CF(N)=ESMF_ConfigCreate(rc=RC)
!
        WRITE(INT_TO_CHAR,FMT)N
        CONFIG_FILE_NAME='configure_file_'//INT_TO_CHAR                    !<-- Prepare the config file names
!
        CALL ESMF_ConfigLoadFile(config  =CF(N)                         &
                                ,filename=CONFIG_FILE_NAME              &
                                ,rc      =RC)
        IF(RC/=0)EXIT                                                      !<-- Exit loop after running out of config files
      ENDDO
!
!
!-----------------------------------------------------------------------
!***  WE MUST LOOP THROUGH THE CONFIGURE FILES OF ALL OF THE CURRENT
!***  DOMAIN'S ANCESTORS TO COLLECT INFORMATION NEEDED TO PROPERLY
!***  DESCRIBE THE CURRENT GRID.  THIS IS NECESSARY BECAUSE ALL
!***  GRIDS' ROWS AND COLUMNS LIE PARALLEL TO THOSE OF THE UPPERMOST
!***  GRID.
!-----------------------------------------------------------------------
!
      ID_DOMAIN=DOMAIN_ID_MINE
!
      N=0
!
!-----------------------------------------------------------------------
      main_loop: DO
!-----------------------------------------------------------------------
!
        N=N+1
!
!-----------------------------
!***  Domain IDs of Ancestors
!-----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Domain ID of Ancestor"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOMAIN)               &  !<-- The config object
                                    ,value =ID_ANCESTOR                 &  !<-- The variable filled 
                                    ,label ='my_parent_id:'             &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------
!***  SW Corner Locations on Ancestors
!--------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get SW Corner I and J on Ancestor Grid"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOMAIN)               &  !<-- The config object
                                    ,value =I_START_SW                  &  !<-- The variable filled 
                                    ,label ='i_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOMAIN)               &  !<-- The config object
                                    ,value =J_START_SW                  &  !<-- The variable filled 
                                    ,label ='j_parent_start:'           &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------
!***  Parent-to-Child Ratios
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Child-to-Parent Ratio of Ancestor"  
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ConfigGetAttribute(config=CF(ID_DOMAIN)               &  !<-- The config object
                                    ,value =PARENT_CHILD_SPACE_RATIO(N) &  !<-- The variable filled 
                                    ,label ='parent_child_space_ratio:' &  !<-- Give this label's value to the previous variable
                                    ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
        ID_ANCESTORS(N)=ID_ANCESTOR                                        !<-- Store domain IDs of all ancestors
        SW_CORNER(1,N)=I_START_SW                                          !<-- Store parent I of SW corner of its child
        SW_CORNER(2,N)=J_START_SW                                          !<-- Store parent J of SW corner of its child
!
        IF(ID_ANCESTOR==1)EXIT                                             !<-- We have reached the uppermost domain
!
        ID_DOMAIN=ID_ANCESTOR
!
!-----------------------------------------------------------------------
!
      ENDDO main_loop
!
!-----------------------------------------------------------------------
!
      NUM_ANCESTORS=N                                                        !<-- How many ancestors are there?
!
!-----------------------------------------------------------------------
!***  ROWS AND COLUMNS OF ALL NESTS' GRIDS LIE PARALLEL TO THOSE OF 
!***  UPPERMOST PARENT GRID.  THUS THE CENTRAL ROTATED LATITUDE AND
!***  LONGITUDE OF ALL NESTS MUST BE THOSE OF THE UPPERMOST DOMAIN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Central Lat/Lon of Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =TPH0D                         &  !<-- The variable filled 
                                  ,label ='tph0d:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =TLM0D                         &  !<-- The variable filled 
                                  ,label ='tlm0d:'                      &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  GET DIMENSIONS OF UPPERMOST DOMAIN AS THE BASELINE.
!***  WE MUST ALSO KNOW SOUTHERN AND WESTERN BOUNDARY LOCATIONS
!***  AS WELL AS WHETHER IT IS GLOBAL OR NOT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Baseline Dimensions of Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =IM_1                          &  !<-- The variable filled 
                                  ,label ='im:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =JM_1                          &  !<-- The variable filled 
                                  ,label ='jm:'                         &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Southern/Western Boundary of Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =SBD_1                         &  !<-- The variable filled 
                                  ,label ='sbd:'                        &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =WBD_1                         &  !<-- The variable filled 
                                  ,label ='wbd:'                        &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Global Flag for Uppermost Domain"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF(ID_ANCESTORS(N))           &  !<-- The config object
                                  ,value =GLOBAL                        &  !<-- The variable filled 
                                  ,label ='global:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_SET)   
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  FULL GRID DIMENSIONS; DELTA PHI AND DELTA LAMBDA 
!***  FOR UPPERMOST DOMAIN.
!-----------------------------------------------------------------------
!
      IF(TRIM(GLOBAL)=='true')THEN                                         !<-- Uppermost domain is global
        IDE_1=IM_1+2
        JDE_1=JM_1+2
        DPHD_1=-SBD_1*2./REAL(JDE_1-3)
        DLMD_1=-WBD_1*2./REAL(IDE_1-3)
        TPH_BASE_1=SBD_1-2.*DPHD_1
        TLM_BASE_1=WBD_1-2.*DLMD_1
      ELSE                                                                 !<-- Uppermost domain is regional
        IDE_1=IM_1  
        JDE_1=JM_1  
        DPHD_1=-SBD_1*2./REAL(JDE_1-1)
        DLMD_1=-WBD_1*2./REAL(IDE_1-1)
        TPH_BASE_1=SBD_1-DPHD_1
        TLM_BASE_1=WBD_1-DLMD_1
      ENDIF
!
!-----------------------------------------------------------------------
!***  LOOP THROUGH THIS NEST'S ANCESTORS IN ORDER TO OBTAIN ITS:
!***  (1) DELTA PHI AND DELTA LAMBDA
!***  (2) SOUTHERN/WESTERN BOUNDARY LOCATIONS 
!
!***  WE MUST WORK DOWNWARD THROUGH THE ANCESTORS BECAUSE
!***  THE UPPERMOST DOMAIN IS THE FOUNDATION.
!
!-----------------------------------------------------------------------
!
      DPHD_X=DPHD_1
      DLMD_X=DLMD_1
      TPH_BASE=TPH_BASE_1
      TLM_BASE=TLM_BASE_1
!
!-----------------------------------------------------------------------
!
      work_loop: DO N=NUM_ANCESTORS,1,-1
!
        I_START_SW=SW_CORNER(1,N)
        J_START_SW=SW_CORNER(2,N)
!
        SBD_X=TPH_BASE+J_START_SW*DPHD_X                                   !<-- Southern boundary of ancestor N
        WBD_X=TLM_BASE+I_START_SW*DLMD_X                                   !<-- Western boundary of ancestor N
!
        DPHD_X=DPHD_X/REAL(PARENT_CHILD_SPACE_RATIO(N))                    !<-- Delta phi for child of ancestor N
        DLMD_X=DLMD_X/REAL(PARENT_CHILD_SPACE_RATIO(N))                    !<-- Delta lambda for child of ancestor N
!
        TPH_BASE=SBD_X-DPHD_X
        TLM_BASE=WBD_X-DLMD_X
!
      ENDDO work_loop
!
      DPHD_MINE=DPHD_X
      DLMD_MINE=DLMD_X
!!!   SBD_MINE=SBD_X
!!!   WBD_MINE=WBD_X
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE SET_NEST_GRIDS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PARENT_DATA_TO_DOMAIN(EXP_STATE_DYN                    &
                                      ,EXP_STATE_PHY                    &
                                      ,EXP_STATE_DOMAIN )
!
!-----------------------------------------------------------------------
!***  Transfer from Dynamics/Physics export states to DOMAIN export 
!***  state the data needed for parent generation of child boundary 
!***  data.
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      TYPE(ESMF_State),INTENT(IN)  :: EXP_STATE_DYN                     &   !<-- Dynamics export state
                                     ,EXP_STATE_PHY                         !<-- Physics export state
!
      TYPE(ESMF_State),INTENT(INOUT) :: EXP_STATE_DOMAIN                    !<-- DOMAIN export state into which fcst Arrays are transferred
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: ITS,ITE,JTS,JTE                                        &
                ,IDS,IDE,JDS,JDE                                        &
                ,INDX_CW,INDX_Q                                         &
                ,LNSH,LNSV                                              &      
                ,NHALO
!
      INTEGER :: RC,RC_TRANS
!
      REAL :: PDTOP,PT
!
      REAL,DIMENSION(:),ALLOCATABLE :: ARRAY_1D
!
      TYPE(ESMF_StateItemType) :: STATEITEMTYPE
!
      TYPE(ESMF_Field) :: HOLD_FIELD
!
      integer :: n
!
      CHARACTER(LEN=8), DIMENSION(5) :: EXP_FIELD
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC_TRANS=ESMF_SUCCESS
!
#if 1
      EXP_FIELD = (/ 'PD      '                                         &
                    ,'T       '                                         &
                    ,'U       '                                         &
                    ,'V       '                                         &
                    ,'TRACERS '                                         &
                               /)
!
!-----------------------------------------------------------------------
!
      item_loop: DO N=1,SIZE(EXP_FIELD)
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Test the presence of "//TRIM(EXP_FIELD(N))//" in Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state        =EXP_STATE_DYN                  &  !<-- The Dynamics export state
                          ,name         =TRIM(EXP_FIELD(N))             &  !<-- Check presence of this Field
                          ,stateItemType=STATEITEMTYPE                  &  !<-- ESMF Type of the Field
                          ,rc           =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        IF (STATEITEMTYPE /= ESMF_STATEITEM_NOTFOUND) THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract "//TRIM(EXP_FIELD(N))//" from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_StateGet(state   =EXP_STATE_DYN                     &  !<-- The Dynamics export state
                            ,itemName=TRIM(EXP_FIELD(N))                &  !<-- Extract this Field
                            ,field   =HOLD_FIELD                        &  !<-- Put the extracted Field here
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ELSE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract "//TRIM(EXP_FIELD(N))//" from Physics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_StateGet(state   =EXP_STATE_PHY                     &  !<-- The Physics export state
                            ,itemName=TRIM(EXP_FIELD(N))                &  !<-- Extract this Field
                            ,field   =HOLD_FIELD                        &  !<-- Put the extracted Field here
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        END IF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert PD into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateAdd(state=EXP_STATE_DOMAIN                       &  !<-- Insert PD into DOMAIN export state
                          ,field=HOLD_FIELD                             &  !<-- The Field to be inserted
                          ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      END DO  item_loop
#else
!-----------------
!***  Transfer PD 
!-----------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract PD from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DYN                         &  !<-- The Dynamics export state
                        ,itemName='PD'                                  &  !<-- Extract PD
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert PD into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=EXP_STATE_DOMAIN                         &  !<-- Insert PD into DOMAIN export state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------
!***  Transfer Temperature
!--------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract T from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DYN                         &  !<-- The Dynamics export state
                        ,itemName='T'                                   &  !<-- Extract T
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert T into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=EXP_STATE_DOMAIN                         &  !<-- Insert T into DOMAIN export state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------
!***  Transfer U Wind
!---------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract U from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DYN                         &  !<-- The Dynamics export state
                        ,itemName='U'                                   &  !<-- Extract U
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert U into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=EXP_STATE_DOMAIN                         &  !<-- Insert U into DOMAIN export state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------
!***  Transfer V Wind
!---------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract V from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DYN                         &  !<-- The Dynamics export state
                        ,itemName='V'                                   &  !<-- Extract V
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert V into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=EXP_STATE_DOMAIN                         &  !<-- Insert V into DOMAIN export state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------
!***  Transfer TRACER Array
!---------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract TRACERS from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DYN                         &  !<-- The Dynamics export state
                        ,itemName='TRACERS'                             &  !<-- Extract Tracers
                        ,field   =HOLD_FIELD                            &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert TRACERS into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=EXP_STATE_DOMAIN                         &  !<-- Insert Tracers into DOMAIN export state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
#endif
!---------------------
!***  Transfer INDX_Q
!---------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract INDX_Q from Physics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_PHY                        &  !<-- The Physics export state
                            ,name ='INDX_Q'                             &  !<-- Name of Attribute to extract
                            ,value=INDX_Q                               &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert INDX_Q into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='INDX_Q'                             &  !<-- The name of the Attribute to insert
                            ,value=INDX_Q                               &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------
!***  Transfer INDX_CW
!----------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract INDX_CW from Physics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_PHY                        &  !<-- The Physics export state
                            ,name ='INDX_CW'                            &  !<-- Name of Attribute to extract
                            ,value=INDX_CW                              &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert INDX_CW into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='INDX_CW'                            &  !<-- The name of the Attribute to insert
                            ,value=INDX_CW                              &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------
!***  Transfer FIS
!------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract FIS from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state   =EXP_STATE_DYN                     &  !<-- The Dynamics export state
                        ,itemName='FIS'                             &  !<-- Extract FIS
                        ,field   =HOLD_FIELD                        &  !<-- Put the extracted Field here
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert FIS into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateAdd(state=EXP_STATE_DOMAIN                         &  !<-- Insert FIS into DOMAIN export state
                        ,field=HOLD_FIELD                               &  !<-- The Field to be inserted
                        ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------
!***  LM is needed but not transferred
!--------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract LM from DYN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='LM'                                 &  !<-- Name of Attribute to extract
                            ,value=LM                                   &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------------
!***  Transfer PT,PDTOP,PSGML1,SG1,SG2,SGML2
!--------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert PT into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='PT'                                 &  !<-- Name of Attribute to extract
                            ,value=PT                                   &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='PT'                                 &  !<-- The name of the Attribute to insert
                            ,value=PT                                   &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='PDTOP'                              &  !<-- Name of Attribute to extract
                            ,value=PDTOP                                &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='PDTOP'                              &  !<-- Name of Attribute to extract
                            ,value=PDTOP                                &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      ALLOCATE(ARRAY_1D(1:LM))
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DYN                    &  !<-- The Dynamics export state
                            ,name     ='PSGML1'                         &  !<-- Extract PGMSL1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The DOMAIN export state
                            ,name     ='PSGML1'                         &  !<-- Extract PGMSL1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DYN                    &  !<-- The Dynamics export state
                            ,name     ='SGML2'                          &  !<-- Extract PGMSL1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The DOMAIN export state
                            ,name     ='SGML2'                          &  !<-- Extract PGMSL1
                            ,count    =LM                               &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      DEALLOCATE(ARRAY_1D)
      ALLOCATE(ARRAY_1D(1:LM+1))
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DYN                    &  !<-- The Dynamics export state
                            ,name     ='SG1'                            &  !<-- Extract PGMSL1
                            ,count    =LM+1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The DOMAIN export state
                            ,name     ='SG1'                            &  !<-- Extract PGMSL1
                            ,count    =LM+1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeGet(state    =EXP_STATE_DYN                    &  !<-- The Dynamics export state
                            ,name     ='SG2'                            &  !<-- Extract PGMSL1
                            ,count    =LM+1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      CALL ESMF_AttributeSet(state    =EXP_STATE_DOMAIN                 &  !<-- The DOMAIN export state
                            ,name     ='SG2'                            &  !<-- Extract PGMSL1
                            ,count    =LM+1                             &  !<-- # of words in data list
                            ,valueList=ARRAY_1D                         &  !<-- Put extracted values ehre
                            ,rc       =RC)
!
      DEALLOCATE(ARRAY_1D)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-------------------------------------------
!***  Transfer Subdomain Integration Limits
!***  and Full Domain Limits
!-------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Integration Limits from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='ITS'                                &  !<-- Name of Attribute to extract
                            ,value=ITS                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='ITE'                                &  !<-- Name of Attribute to extract
                            ,value=ITE                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='JTS'                                &  !<-- Name of Attribute to extract
                            ,value=JTS                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='JTE'                                &  !<-- Name of Attribute to extract
                            ,value=JTE                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='LM'                                 &  !<-- Name of Attribute to extract
                            ,value=LM                                   &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='NHALO'                              &  !<-- Name of Attribute to extract
                            ,value=NHALO                                &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='IDS'                                &  !<-- Name of Attribute to extract
                            ,value=IDS                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='IDE'                                &  !<-- Name of Attribute to extract
                            ,value=IDE                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='JDS'                                &  !<-- Name of Attribute to extract
                            ,value=JDS                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='JDE'                                &  !<-- Name of Attribute to extract
                            ,value=JDE                                  &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert Integration Limits into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='ITS'                                &  !<-- The name of the Attribute to insert
                            ,value=ITS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='ITE'                                &  !<-- The name of the Attribute to insert
                            ,value=ITE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JTS'                                &  !<-- The name of the Attribute to insert
                            ,value=JTS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JTE'                                &  !<-- The name of the Attribute to insert
                            ,value=JTE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='LM'                                 &  !<-- The name of the Attribute to insert
                            ,value=LM                                   &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='NHALO'                              &  !<-- The name of the Attribute to insert
                            ,value=NHALO                                &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='IDS'                                &  !<-- The name of the Attribute to insert
                            ,value=IDS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='IDE'                                &  !<-- The name of the Attribute to insert
                            ,value=IDE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JDS'                                &  !<-- The name of the Attribute to insert
                            ,value=JDS                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='JDE'                                &  !<-- The name of the Attribute to insert
                            ,value=JDE                                  &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------------------------
!***  Transfer Width of Blending Region Along Boundaries
!--------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract LNSH,LNSV from Dynamics Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='LNSH'                               &  !<-- Name of Attribute to extract
                            ,value=LNSH                                 &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
      CALL ESMF_AttributeGet(state=EXP_STATE_DYN                        &  !<-- The Dynamics export state
                            ,name ='LNSV'                               &  !<-- Name of Attribute to extract
                            ,value=LNSV                                 &  !<-- Put the extracted Attribute here
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Insert LNSH,LNSV into DOMAIN Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='LNSH'                               &  !<-- The name of the Attribute to insert
                            ,value=LNSH                                 &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
      CALL ESMF_AttributeSet(state=EXP_STATE_DOMAIN                     &  !<-- The DOMAIN export state
                            ,name ='LNSV'                               &  !<-- The name of the Attribute to insert
                            ,value=LNSV                                 &  !<-- The Attribute to be inserted
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_TRANS)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PARENT_DATA_TO_DOMAIN
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
!***  THIS ROUTINE MOVES NEW BOUNDARY DATA FOR NESTED DOMAINS FROM
!***  ONE IMPORT/EXPORT STATE TO ANOTHER.
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
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE SIDES_1D_REAL
        REAL,DIMENSION(:),ALLOCATABLE :: SOUTH
        REAL,DIMENSION(:),ALLOCATABLE :: NORTH
        REAL,DIMENSION(:),ALLOCATABLE :: WEST
        REAL,DIMENSION(:),ALLOCATABLE :: EAST
      END TYPE SIDES_1D_REAL
!
      INTEGER :: ISTAT,KOUNT,NTIMESTEP,RC,RC_BND_MV
!
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF
!
      TYPE(SIDES_1D_REAL),SAVE :: BOUNDARY_H                            &
                                 ,BOUNDARY_V
!
      integer :: nnnn
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  IF THE CLOCK WAS SENT IN THEN THIS TRANSFER OF DATA IS 
!***  DEPENDENT UPON THE TIMESTEP.  
!-----------------------------------------------------------------------
!
      IF(PRESENT(CLOCK))THEN
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="State_to_State Gets DOMAIN Clock"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_ClockGet(clock       =CLOCK                           &  !<-- The ESMF clock
                          ,advanceCount=NTIMESTEP_ESMF                  &  !<-- The number of times the clock has been advanced
                          ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        NTIMESTEP=NTIMESTEP_ESMF
!
        IF(MOD(NTIMESTEP,RATIO)/=0)RETURN
!
      ENDIF
!-----------------------------------------------------------------------
!***  CHECK EACH SIDE OF THE CHILD BOUNDARY.  IF DATA IS PRESENT FROM
!***  THAT SIDE IN THE INPUT STATE THEN MOVE IT TO THE OUTPUT STATE.
!-----------------------------------------------------------------------
!
!-------------
!***  South H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Input State for South H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &   !<-- Look at the input state
                            ,name ='SOUTH_H'                            &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      south_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => South boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%SOUTH))THEN
          ALLOCATE(BOUNDARY_H%SOUTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%SOUTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South H Data from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =STATE_IN                       &   !<-- Extract data from input state
                              ,name     ='SOUTH_H'                      &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert South H Data into Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =STATE_OUT                      &   !<-- Insert data into output state
                              ,name     ='SOUTH_H'                      &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF south_h
!
!-------------
!***  South V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Input State for South V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &   !<-- Look at the input state
                            ,name ='SOUTH_V'                            &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      south_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => South boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%SOUTH))THEN
          ALLOCATE(BOUNDARY_V%SOUTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%SOUTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract South V Data from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =STATE_IN                       &   !<-- Extract data from input state
                              ,name     ='SOUTH_V'                      &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert South V Data into Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =STATE_OUT                      &   !<-- Insert data into output state
                              ,name     ='SOUTH_V'                      &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%SOUTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF south_v
!
!-------------
!***  North H
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Input State for North H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &   !<-- Look at input state
                            ,name ='NORTH_H'                            &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      north_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => North boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%NORTH))THEN
          ALLOCATE(BOUNDARY_H%NORTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%NORTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North H Data from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =STATE_IN                       &   !<-- Extract data from input state
                              ,name     ='NORTH_H'                      &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert North H Data into Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =STATE_OUT                      &   !<-- Insert data into output state
                              ,name     ='NORTH_H'                      &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF north_h
!
!-------------
!***  North V
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Input State for North V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &   !<-- Look at the input state
                            ,name ='NORTH_V'                            &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      north_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                         !<-- True => North boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%NORTH))THEN
          ALLOCATE(BOUNDARY_V%NORTH(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%NORTH stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract North V Data from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =STATE_IN                       &   !<-- Extract data from input state
                              ,name     ='NORTH_V'                      &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert North V Data into Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =STATE_OUT                      &   !<-- Insert data into output state
                              ,name     ='NORTH_V'                      &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%NORTH               &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF north_v
!
!------------
!***  West H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Input State for West H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &   !<-- Look at the input state
                            ,name ='WEST_H'                             &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      west_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => West boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%WEST))THEN
          ALLOCATE(BOUNDARY_H%WEST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%WEST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West H Data from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =STATE_IN                       &   !<-- Extract data from input state
                              ,name     ='WEST_H'                       &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert West H Data into Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =STATE_OUT                      &   !<-- Insert data into output state
                              ,name     ='WEST_H'                       &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF west_h
!
!------------
!***  West V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Input State for West V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &   !<-- Look at the input state
                            ,name ='WEST_V'                             &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      west_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => West boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%WEST))THEN
          ALLOCATE(BOUNDARY_V%WEST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%WEST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract West V Data from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =STATE_IN                       &   !<-- Extract data input state
                              ,name     ='WEST_V'                       &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert West V Data into Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =STATE_OUT                      &   !<-- Insert data into output state
                              ,name     ='WEST_V'                       &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%WEST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF west_v
!
!------------
!***  East H
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Input State for East H Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &   !<-- Look at the input state
                            ,name ='EAST_H'                             &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      east_h: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => East boundary H point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_H%EAST))THEN
          ALLOCATE(BOUNDARY_H%EAST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_H%EAST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East H Data from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =STATE_IN                       &   !<-- Extract data from input state
                              ,name     ='EAST_H'                       &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert East H Data into Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =STATE_OUT                      &   !<-- Insert data into output state
                              ,name     ='EAST_H'                       &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_H%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF east_h
!
!------------
!***  East V
!------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Check Input State for East V Data"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_AttributeGet(state=STATE_IN                             &   !<-- Look at the input state
                            ,name ='EAST_V'                             &   !<-- Is this name present?
                            ,count=KOUNT                                &   !<-- How many items present?
                            ,rc   =RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!!!   CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      east_v: IF(KOUNT>0.AND.RC==ESMF_SUCCESS)THEN                          !<-- True => East boundary V point data is present
!
        IF(.NOT.ALLOCATED(BOUNDARY_V%EAST))THEN
          ALLOCATE(BOUNDARY_V%EAST(1:KOUNT),stat=ISTAT)
          IF(ISTAT/=0)WRITE(0,*)' Failed to allocate BOUNDARY_V%EAST stat=',ISTAT
        ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract East V Data from Input State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(state    =STATE_IN                       &   !<-- Extract data from input state
                              ,name     ='EAST_V'                       &   !<-- The name of the data
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Insert East V Data into Output State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeSet(state    =STATE_OUT                      &   !<-- Insert data into output state
                              ,name     ='EAST_V'                       &   !<-- The name of the data 
                              ,count    =KOUNT                          &   !<-- The data has this many items
                              ,valueList=BOUNDARY_V%EAST                &   !<-- The new combined boundary data
                              ,rc=RC )
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_BND_MV)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF east_v
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE BOUNDARY_DATA_STATE_TO_STATE
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
! 
      SUBROUTINE WATERFALLS(FIS                                         &
                           ,SEA_MASK                                    &
                           ,LOWER_TOPO                                  &
                           ,IDS,IDE,JDS,JDE)
!
!-----------------------------------------------------------------------
!***  WHEN A PARENT INITIALIZES ITS CHILD, THE SEA MASK HAD TO BE DONE
!***  WITH NEAREST NEIGHBOR LOGIC WHILE FIS SHOULD BE DONE BILINEARLY.
!***  THIS CAN LEAD TO ADJACENT WATER POINTS HAVING DIFFERENT VALUES
!***  OF FIS.  WHEN THAT IS THE CASE, MAKE THE ELEVATION OF ALL
!***  ADJACENT WATER POINTS EQUAL TO THE LOWEST OF THEIR VALUES.
!***  SAVE THE I,J OF ALL LOWERED POINTS SO THE ATMOSPHERIC COLUMN
!***  CAN ULTIMATELY BE ADJUSTED.
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),INTENT(IN) :: IDS,IDE,JDS,JDE                     !<-- Lateral dimensions of nest grid
!
      REAL(KIND=KFPT),DIMENSION(IDS:IDE,JDS:JDE),INTENT(IN) :: SEA_MASK    !<-- Sea mask of nest grid points
!
      REAL(KIND=KFPT),DIMENSION(IDS:IDE,JDS:JDE,1),INTENT(INOUT) :: FIS    !<-- Sfc geopotential on nest grid points
!
      LOGICAL,DIMENSION(IDS:IDE,JDS:JDE),INTENT(OUT) :: LOWER_TOPO         !<-- Flag points where topography is lowered
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT) :: I,J
      INTEGER(KIND=KINT) :: ITER,KOUNT_CHANGE
!
      REAL(KIND=KFPT) :: FIS_0                                          &
                        ,FIS_E,FIS_N,FIS_W,FIS_S                        &
                        ,FIS_NE,FIS_NW,FIS_SW,FIS_SE                    &
                        ,FIS_NEW
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      iter_loop: DO ITER=1,500
!
        KOUNT_CHANGE=0
!
!-----------------------------------------------------------------------
!
        DO J=JDS,JDE
        DO I=IDS,IDE
!
          IF(SEA_MASK(I,J)<0.01)CYCLE                                      !<-- We are adjusting only water points 
!
!-----------------------------------------------------------------------
!
          FIS_0=FIS(I,J,1)
!
!----------
!***  East
!----------
!
          FIS_E=FIS_0
!
          IF(I+1<=IDE)THEN
            IF(SEA_MASK(I+1,J)>0.99)FIS_E=FIS(I+1,J,1)
          ENDIF
!
!---------------
!***  Northeast
!---------------
!
          FIS_NE=FIS_0
!
          IF(I+1<=IDE.AND.J+1<=JDE)THEN
            IF(SEA_MASK(I+1,J+1)>0.99)FIS_NE=FIS(I+1,J+1,1)
          ENDIF
!
!-----------
!***  North
!-----------
!
          FIS_N=FIS_0
!
          IF(J+1<=JDE)THEN
            IF(SEA_MASK(I,J+1)>0.99)FIS_N=FIS(I,J+1,1)
          ENDIF
!
!---------------
!***  Northwest
!---------------
!
          FIS_NW=FIS_0
!
          IF(I-1>=IDS.AND.J+1<=JDE)THEN
            IF(SEA_MASK(I-1,J+1)>0.99)FIS_NW=FIS(I-1,J+1,1)
          ENDIF
!
!----------
!***  West
!----------
!
          FIS_W=FIS_0
!
          IF(I-1>=IDS)THEN
            IF(SEA_MASK(I-1,J)>0.99)FIS_W=FIS(I-1,J,1)
          ENDIF
!
!---------------
!***  Southwest
!---------------
!
          FIS_SW=FIS_0
!
          IF(I-1>=IDS.AND.J-1>=JDS)THEN
            IF(SEA_MASK(I-1,J-1)>0.99)FIS_SW=FIS(I-1,J-1,1)
          ENDIF
!
!-----------
!***  South
!-----------
!
          FIS_S=FIS_0
!
          IF(J-1>=JDS)THEN
            IF(SEA_MASK(I,J-1)>0.99)FIS_S=FIS(I,J-1,1)
          ENDIF
!
!---------------
!***  Southeast
!---------------
!
          FIS_SE=FIS_0
!
          IF(I+1<=IDE.AND.J-1>=JDS)THEN
            IF(SEA_MASK(I+1,J-1)>0.99)FIS_SE=FIS(I+1,J-1,1)
          ENDIF
!
!-----------------------------------------------------------------------
!***  LOWER THE POINT IN QUESTION TO THE LOWEST VALUE OF ITSELF AND 
!***  ITS NEIGHBORS IF IT IS A WATER POINT.
!***  ALSO SAVE ALL I,J LOCATIONS WHERE FIS IS CHANGED SO THAT WE
!***  CAN ADJUST THE ATMOSPHERIC COLUMN APPROPRIATELY LATER.
!-----------------------------------------------------------------------
!
          FIS_NEW=MIN(FIS_0                                             &
                     ,FIS_E,FIS_N,FIS_W,FIS_E                           &
                     ,FIS_NE,FIS_NW,FIS_SW,FIS_SE)
!
          IF(FIS_NEW+0.1<FIS_0)THEN
            KOUNT_CHANGE=KOUNT_CHANGE+1
            FIS(I,J,1)=FIS_NEW
            LOWER_TOPO(I,J)=.TRUE.
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO
        ENDDO
!
!
        IF(KOUNT_CHANGE==0)EXIT iter_loop
        IF(ITER==100)THEN
          WRITE(0,*)' Reached 100 iterations and KOUNT_CHANGE='         &
                   ,KOUNT_CHANGE
        ENDIF
!
      ENDDO iter_loop
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WATERFALLS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
! 
      SUBROUTINE ADJUST_COLUMNS(PD_NEAREST                              &
                               ,PD_BILINEAR                             &
                               ,LOWER_TOPO                              &
                               ,DUMMY_3D                                &
                               ,PT                                      &
                               ,PDTOP                                   &
                               ,SG1                                     &
                               ,SG2                                     &
                               ,IM_CHILD                                &
                               ,JM_CHILD                                &
                                         )
!
!-----------------------------------------------------------------------
!***  WHEN THE SURFACE ELEVATION OF A NESTED DOMAIN IS CHANGED DUE TO
!***  LEVELING OF ADJACENT WATER POINTS, ADJUST THE ATMOSPHERIC COLUMN
!***  AT EACH OF THOSE POINTS.
!-----------------------------------------------------------------------
!
!--------
!***  In
!--------
!
      INTEGER,INTENT(IN) :: IM_CHILD                                    &
                           ,JM_CHILD
!
      REAL,INTENT(IN) :: PDTOP,PT
!
      REAL,DIMENSION(1:LM+1),INTENT(IN) :: SG1,SG2
!
      REAL,DIMENSION(1:IM_CHILD,1:JM_CHILD,1),INTENT(IN) :: PD_NEAREST  &
                                                           ,PD_BILINEAR
!
      LOGICAL,DIMENSION(1:IM_CHILD,1:JM_CHILD),INTENT(IN) :: LOWER_TOPO
!
!-----------
!***  Inout
!-----------
!
      REAL,DIMENSION(1:IM_CHILD,1:JM_CHILD,LM),INTENT(INOUT) :: DUMMY_3D
!
!---------------------
!***  Local Variables
!---------------------
!
      INTEGER :: I,J,L,NUM_LEVS_SPLINE
!
      REAL :: COEFF_1,DELP_EXTRAP,PBOT_IN,PBOT_TARGET,PDTOP_PT          &
             ,PTOP_IN,PTOP_TARGET,R_DELP
!
      REAL,DIMENSION(1:LM) :: PMID_TARGET                               &
                             ,VBL_COLUMN
!
      REAL,DIMENSION(1:LM+1) :: PMID_IN                                 &
                               ,SEC_DERIV                               &
                               ,VBL_INPUT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      DO L=1,LM+1
        SEC_DERIV(L)=0.
      ENDDO
!
      NUM_LEVS_SPLINE=LM+1
!
!-----------------------------------------------------------------------
!***  COMPUTE THE INPUT AND TARGET MID-LAYER PRESSURES FOR THE
!***  SPLINE INTERPOLATION.
!-----------------------------------------------------------------------
!
      DO J=1,JM_CHILD
      DO I=1,IM_CHILD
!
!-----------------------------------------------------------------------
!
        adjust: IF(LOWER_TOPO(I,J))THEN
!
          PTOP_IN=SG2(1)*PD_BILINEAR(I,J,1)+SG1(1)*PDTOP+PT
          PTOP_TARGET=SG2(1)*PD_NEAREST(I,J,1)+SG1(1)*PDTOP+PT
!
          DO L=1,LM
!
            PDTOP_PT=SG1(L+1)*PDTOP+PT
!
            PBOT_IN   =SG2(L+1)*PD_BILINEAR(I,J,1)+PDTOP_PT
            PMID_IN(L)=0.5*(PTOP_IN+PBOT_IN)
            PTOP_IN   =PBOT_IN
            VBL_INPUT(L)=DUMMY_3D(I,J,L)
!
            PBOT_TARGET   =SG2(L+1)*PD_NEAREST(I,J,1)+PDTOP_PT
            PMID_TARGET(L)=0.5*(PTOP_TARGET+PBOT_TARGET)
            PTOP_TARGET   =PBOT_TARGET
!
          ENDDO
!
!***  WE KNOW THE TARGET MID-LAYER PRESSURE IS GREATER THAN THAT IN
!***  THE ORIGINAL COLUMN SINCE THE SFC ELEVATION HAS BEEN LOWERED.
!***  ADD A NEW INPUT LEVEL BY EXTRAPOLATING LINEARLY IN PRESSURE
!***  TO OBTAIN A VALUE AT THE LOWEST OUTPUT MID-LAYER THEN FILL
!***  IN ALL OUTPUT LEVELS WITH THE SPLINE.
!
          PMID_IN(LM+1)=PMID_TARGET(LM)
          R_DELP=1./(PMID_IN(LM)-PMID_IN(LM-1))
          DELP_EXTRAP=PMID_TARGET(LM)-PMID_IN(LM)
          COEFF_1=(VBL_INPUT(LM)-VBL_INPUT(LM-1))*R_DELP
          VBL_INPUT(LM+1)=VBL_INPUT(LM)+COEFF_1*DELP_EXTRAP                !<-- Create extrapolated value at nest's lowest mid-layer
!                                                                               in input array.
!-----------------------------------------------------------------------
!
          IF(ABS(PMID_IN(LM+1)-PMID_IN(LM))<10.)EXIT
!
          CALL SPLINE(NUM_LEVS_SPLINE                                   &  !<-- # of input levels
                     ,PMID_IN                                           &  !<-- Input mid-layer pressures
                     ,VBL_INPUT                                         &  !<-- Input mid-layer mass variable value
                     ,SEC_DERIV                                         &  !<-- Specified 2nd derivatives (=0) at input levels
                     ,LM                                                &  !<-- # of mid-layers to which to interpolate
                     ,PMID_TARGET                                       &  !<-- Mid-layer pressures to which to interpolate 
                     ,VBL_COLUMN )                                         !<-- Mid-layer variable value returned
!
          DO L=1,LM
            DUMMY_3D(I,J,L)=VBL_COLUMN(L)
          ENDDO
!
        ENDIF adjust
!
!-----------------------------------------------------------------------
!
      ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ADJUST_COLUMNS
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
      END MODULE MODULE_NESTING
!
!-----------------------------------------------------------------------
