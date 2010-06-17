!-----------------------------------------------------------------------
!
      MODULE module_NMM_CORE_SETUP
!
!-----------------------------------------------------------------------
!
!***  This module holds the Dynamics Register, Init, Run, and Finalize
!***  routines.  They are called from the DOMAIN component
!***  (DOMAIN_INITIALIZE calls DYNAMICS_INITIALIZE, etc.)
!***  IN module_DOMAIN_GRID_COMP.F90.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE module_DYNAMICS_INTERNAL_STATE                                  !<-- Horizontal loop limits obtained here
!
      USE module_GET_CONFIG_DYN
      USE module_GET_CONFIG_PHY
      USE module_GET_CONFIG_WRITE
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE module_INCLUDE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: NMM_SETUP       !<-- An NMM-specific routine to set up parallelism and ESMF Grid
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE NMM_SETUP(MYPE_IN                                      &
                          ,MPI_INTRA                                    &
                          ,CF                                           &
                          ,DOMAIN_GRID_COMP                             &
                          ,DOMAIN_INT_STATE                             &
                          ,GRID_DOMAIN)
! 
!-----------------------------------------------------------------------
!***  This routine contains NMM-specific code for the DOMAIN component:
!***    (1) Setting up distributed memory parallelism in the NMM;
!***    (2) Creating the ESMF Grid for the DOMAIN components;
!***    (3) Sharing local subdomain index limits among tasks.
!-----------------------------------------------------------------------
!
      USE module_DOMAIN_INTERNAL_STATE
!
      USE module_DM_PARALLEL,ONLY : DECOMP                              &
                                   ,LOCAL_ISTART,LOCAL_IEND             &
                                   ,LOCAL_JSTART,LOCAL_JEND             &
                                   ,SETUP_SERVERS
!
      USE module_INCLUDE
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: MYPE_IN                          &  !<-- Each MPI task's rank
                                      ,MPI_INTRA                           !<-- The communicator with the domain's fcst and quilt tasks.
!
      TYPE(ESMF_Config),INTENT(INOUT) :: CF                                !<-- This domain's configure object
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: DOMAIN_GRID_COMP                !<-- The DOMAIN component
!
      TYPE(DOMAIN_INTERNAL_STATE),INTENT(INOUT) :: DOMAIN_INT_STATE        !<-- The DOMAIN Internal State
!
      TYPE(ESMF_Grid),INTENT(OUT) :: GRID_DOMAIN                           !<-- The ESMF Grid for the NMM integration grid
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER(kind=KINT) :: I,IERR,J,K,N,NUM_PES,RC,RC_CORE
!
      INTEGER(kind=KINT) :: IM,JM                                       &  !<-- Horizontal dimensions of the full integration grid
                           ,INPES,JNPES                                 &  !<-- MPI tasks in I and J directions
                           ,LM                                          &  !<-- Number of atmospheric model layers
                           ,MPI_INTRA_B                                 &  !<-- The MPI intra-communicator
                           ,MYPE                                        &  !<-- My MPI task ID
                           ,NUM_PES_FCST                                &  !<-- Number of MPI tasks applied to the forecast
                           ,NUM_PES_TOT                                 &  !<-- Total # of MPI tasks in the job
                           ,WRITE_GROUPS                                &  !<-- Number of groups of write tasks
                           ,WRITE_TASKS_PER_GROUP                          !<-- #of tasks in each write group
!
      INTEGER(kind=KINT),DIMENSION(2) :: I1                             &  !<-- # of I and J points in each fcst task's subdomain
                                        ,MIN,MAX                        &  !<-- Set start/end of each Grid dimension
                                        ,NCOUNTS                           !<-- Array with I/J limits of MPI task subdomains
!
      CHARACTER(50) :: MODE                                                !<-- Flag for global or regional run
!
      LOGICAL(kind=KLOG) :: GLOBAL                                         !<-- .TRUE. => global ; .FALSE. => regional
!
      TYPE(ESMF_VM)       :: VM                                            !<-- The ESMF virtual machine.
      TYPE(ESMF_DELayout) :: MY_DE_LAYOUT                                  !<-- The ESMF layout type array (for tasks).
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      MYPE=MYPE_IN
!
!-----------------------------------------------------------------------
!***  Set up parameters for MPI communications on this domain's grid.
!-----------------------------------------------------------------------
!
      CALL MPI_COMM_SIZE(MPI_INTRA,NUM_PES_TOT,IERR)
!
      NUM_PES=NUM_PES_TOT
!
!-----------------------------------------------------------------------
!***  Establish the task layout including the Write tasks.
!***  The MPI communicator was provided as input and
!***  the forecast tasks in the I and J directions are
!***  extracted from a configure file.
!***  Give those to SETUP_SERVERS which will split the
!***  communicator between Forecast and Quilt/Write tasks.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Get INPES/JNPES from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL MPI_COMM_DUP(MPI_INTRA,MPI_INTRA_B,RC)                          !<-- Use a duplicate of the communicator for safety
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =INPES                         &  !<-- # of fcst tasks in I direction
                                  ,label ='inpes:'                      &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =JNPES                         &  !<-- # of fcst tasks in J direction
                                  ,label ='jnpes:'                      &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set up Quilt/Write task specifications.
!***  First retrieve the task and group counts from the config file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Get Write Task/Group Info from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_GROUPS                         &  !<-- Number of write groups from config file
                                  ,label ='write_groups:'               &
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(CF                                   &  !<-- The configure file
                                  ,WRITE_TASKS_PER_GROUP                &  !<-- Number of write tasks per group from config file
                                  ,label ='write_tasks_per_group:'      &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Segregate the Forecast tasks from the Quilt/Write tasks.
!-----------------------------------------------------------------------
!
      CALL SETUP_SERVERS(MYPE,INPES,JNPES,NUM_PES                       &
                        ,WRITE_GROUPS,WRITE_TASKS_PER_GROUP             &
                        ,MPI_INTRA_B)
!
!***
!***  NOTE: At this point, NUM_PES is the number of Forecast tasks only.
!***
!-----------------------------------------------------------------------
!
      NUM_PES_FCST=INPES*JNPES                                             !<-- Number of forecast tasks
      domain_int_state%NUM_PES_FCST=NUM_PES_FCST                           !<-- Save this for DOMAIN's Run step
!
!-----------------------------------------------------------------------
!***  Allocate and fill the task list that holds the IDs of
!***  the Forecast tasks.
!-----------------------------------------------------------------------
!
      ALLOCATE(domain_int_state%PETLIST_FCST(NUM_PES_FCST))                !<-- Task IDs of the forecast tasks
!
      DO N=0,NUM_PES_FCST-1
        domain_int_state%PETLIST_FCST(N+1)=N                               !<-- Collect just the forecast task IDs
      ENDDO
!
!-----------------------------------------------------------------------
!***  Retrieve the VM (Virtual Machine) of the DOMAIN component.
!***  We need VM now to set up the DE layout.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Retrieve VM from DOMAIN Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=DOMAIN_GRID_COMP                   &  !<-- The DOMAIN component
                           ,vm      =VM                                 &  !<-- The ESMF Virtual Machine
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create DE layout based on the I tasks by J tasks specified in
!***  the config file.
!***  This refers only to Forecast tasks.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                            !<-- Select only the forecast tasks
        MY_DE_LAYOUT=ESMF_DELayoutCreate(            VM                 &  !<-- The ESMF virtual machine
                                        ,deCountList=(/INPES,JNPES/)    &  !<-- User-specified I-task by J-task layout
                                        ,rc         =RC)
      ENDIF
!
!-----------------------------------------------------------------------
!***  Create the ESMF Grid.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Extract the dimensions of the domain from the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Get IM,JM,LM from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =IM                            &  !<-- I dimension of full domain
                                  ,label ='im:'                         &  !<-- The label in the configure file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =JM                            &  !<-- J dimension of full domain
                                  ,label ='jm:'                         &  !<-- The label in the configure file
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =LM                            &  !<-- Vertical dimension of full domain
                                  ,label ='lm:'                         &  !<-- The label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------------------
!***  Retrieve the forecast domain mode from the config file.
!------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Get GLOBAL/REGIONAL Mode from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure file
                                  ,value =MODE                          &  !<-- Flag for global (true) or regional (false) run
                                  ,label ='global:'                     &  !<-- The label in the configure file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(TRIM(MODE)=='true')THEN
        GLOBAL=.TRUE.
      ELSE
        GLOBAL=.FALSE.
      ENDIF
!
!-----------------------------------------------------------------------
!***  If this is a global mode forecast, extend IM and JM.
!***  The first dimension of NCOUNTS is the I dimension for parallelization.
!***  The second dimension of NCOUNTS is the J dimension.
!-----------------------------------------------------------------------
!
      IF(GLOBAL)THEN      !<-- Global mode horizontal dimensions.
        NCOUNTS(1)=IM+2
        NCOUNTS(2)=JM+2
      ELSE                !<-- Regional mode horizontal dimensions.
        NCOUNTS(1)=IM
        NCOUNTS(2)=JM
      ENDIF
!
      MAX(1)=NCOUNTS(1)
      MAX(2)=NCOUNTS(2)
!
      MIN(1)=1
      MIN(2)=1
!
!-----------------------------------------------------------------------
!***  Now create the DOMAIN component's ESMF Grid
!***  for the NMM's integration grid.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="NMM_SETUP: Create the ESMF Grid"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GRID_DOMAIN=ESMF_GridCreateShapeTile(regDecomp     =(/INPES,JNPES/)    &  !<-- I x J task layout
                                       ,minIndex      =(/MIN(1),MIN(2)/)  &  !<-- Min indices in I and J
                                       ,maxIndex      =(/MAX(1),MAX(2)/)  &  !<-- Max indices in I and J
                                       ,gridEdgeLWidth=(/0,0/)            &  !<-- Padding, lower edges for noncentered stagger
                                       ,gridEdgeUWidth=(/0,0/)            &  !<-- Padding, upper edges for noncentered stagger
                                       ,name          ="GRID"             &  !<-- Name of the Grid
                                       ,indexflag     =ESMF_INDEX_GLOBAL  &
                                       ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Get the local array sizes for the DOMAIN Grid.
!***  Only forecast tasks are relevant here.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                              !<-- Select only fcst tasks
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="NMM_SETUP: Get EMSF Sizes of Local Subdomains"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridGet(grid              =GRID_DOMAIN                &
                         ,localDe           =0                          &
                         ,staggerloc        =ESMF_STAGGERLOC_CENTER     &
                         ,computationalCount=I1                         & !<-- # of local points in I and J on each task
                         ,rc                =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!***  Using 'computationalCount' from array I1 obtained in the
!***  previous call, generate all of the local task index limits
!***  for all Forecast tasks.  
!***  The user, not ESMF, does this work.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                            !<-- Select only the forecast tasks
        CALL DECOMP(MYPE,INPES,JNPES,NUM_PES_FCST,IM,JM,LM,GLOBAL,I1)
!
        ALLOCATE(domain_int_state%LOCAL_ISTART(0:NUM_PES_FCST-1))
        ALLOCATE(domain_int_state%LOCAL_IEND  (0:NUM_PES_FCST-1))
        ALLOCATE(domain_int_state%LOCAL_JSTART(0:NUM_PES_FCST-1))
        ALLOCATE(domain_int_state%LOCAL_JEND  (0:NUM_PES_FCST-1))
!
        DO N=0,NUM_PES_FCST-1
          domain_int_state%LOCAL_ISTART(N)=LOCAL_ISTART(N)                 !<-- Starting I for all forecasts' subdomains
          domain_int_state%LOCAL_IEND  (N)=LOCAL_IEND  (N)                 !<-- Ending I for all forecasts' subdomains
          domain_int_state%LOCAL_JSTART(N)=LOCAL_JSTART(N)                 !<-- Starting J for all forecasts' subdomains
          domain_int_state%LOCAL_JEND  (N)=LOCAL_JEND  (N)                 !<-- Ending J for all forecasts' subdomains
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_SETUP
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE module_NMM_CORE_SETUP
!
!-----------------------------------------------------------------------

