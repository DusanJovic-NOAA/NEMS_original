!-----------------------------------------------------------------------
!
      MODULE MODULE_CORE_SETUP
!
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
      USE MODULE_DYNAMICS_INTERNAL_STATE      !<-- Horizontal loop limits obtained here
!
      USE MODULE_GET_CONFIG_DYN
      USE MODULE_GET_CONFIG_PHY
      USE MODULE_GET_CONFIG_WRITE
!
      USE MODULE_CLOCKTIMES
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
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
      PUBLIC :: NMM_SETUP,GFS_SETUP       !<-- An NMM-specific routine to set up parallelism and ESMF Grid
!
!-----------------------------------------------------------------------
      INCLUDE 'kind.inc'
!-----------------------------------------------------------------------
!
!
!
!
!-----------------------------------------------------------------------
!***  FOR DETERMINING CLOCKTIMES OF VARIOUS PIECES OF THE DYNAMICS.
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT) :: btim,btim0
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      SUBROUTINE NMM_SETUP(ATM_GRID_COMP,ATM_INT_STATE,GRID_ATM)
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
      USE MODULE_DM_PARALLEL,ONLY : DECOMP                              &
                                   ,LOCAL_ISTART,LOCAL_IEND             &
                                   ,LOCAL_JSTART,LOCAL_JEND             &
                                   ,SETUP_SERVERS
!
      USE MODULE_INCLUDE
!
!-----------------------------------------------------------------------
!***  SET UP THE ESMF GRID FOR THE NMM AND ESTABLISH THE
!***  DISTRIBUTED MEMORY PARALLELISM.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ATM_INTERNAL_STATE),INTENT(INOUT) :: ATM_INT_STATE             !<-- The ATM Internal State
      TYPE(ESMF_Grid)         ,INTENT(OUT)   :: GRID_ATM                  !<-- The ESMF GRID for the NMM integration grid
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)   :: CF                                           !<-- The config object
      TYPE(ESMF_VM)       :: VM                                           !<-- The ESMF virtual machine.
      TYPE(ESMF_DELayout) :: MY_DE_LAYOUT                                 !<-- The ESMF layout type array (for tasks).
!
      INTEGER(KIND=KINT)  :: INPES,JNPES                                  !<-- MPI tasks in I and J directions
      INTEGER             :: MPI_INTRA,MPI_INTRA_B                        !<-- The MPI intra-communicator
      INTEGER(KIND=KINT)  :: NUM_PES_FCST                                 !<-- Number of MPI tasks applied to the forecast
      INTEGER             :: NUM_PES_TOT                                  !<-- Total # of MPI tasks in the job
      INTEGER             :: WRITE_GROUPS                               & !<-- Number of groups of write tasks
                            ,WRITE_TASKS_PER_GROUP                        !<-- #of tasks in each write group
!
      INTEGER,DIMENSION(2)             :: NCOUNTS                         !<-- Parameter array to set up the
                                                                          !    size of the 2-D ESMF grid.
      INTEGER,DIMENSION(2)             :: I1                              !<-- # of I and J points in each fcst task's subdomain
!
      INTEGER,DIMENSION(2)  :: MIN,MAX                                    !<-- Parameter arrays to set up the
                                                                          !    start number and the end number of
                                                                          !    the ESMF grid in each dimension.
!
      CHARACTER(50)      :: MODE                                          !<-- Flag for global or regional run
!
      LOGICAL            :: GLOBAL                                        !<-- .TRUE. => global ; .FALSE. => regional
!
      INTEGER :: I,J,K,N,NUM_PES,RC,RC_CORE,IM,JM,LM,MYPE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE VM (VIRTUAL MACHINE) OF THE ATM GRIDDED COMPONENT.
!***  CALL ESMF_GridCompGet TO RETRIEVE THE VM ANYWHERE YOU NEED IT.
!***  WE NEED VM NOW TO SET UP THE DE LAYOUT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Retrieve VM from ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,vm      =VM                                 &  !<-- The ESMF Virtual Machine
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  SET UP PARAMETERS FOR MPI COMMUNICATIONS.
!***  USE ESMF UTILITY TO GET PE IDENTIFICATION AND TOTAL NUMBER OF PEs
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get Task IDs and Count"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Local PE rank
                     ,petcount=NUM_PES_TOT                              &  !<-- Total # of tasks (fcst + quilt)
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NUM_PES=NUM_PES_TOT
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIGURE OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Retrieve Config Object from ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ESTABLISH THE TASK LAYOUT INCLUDING THE WRITE TASKS.
!***  GET THE GLOBAL MPI COMMUNICATOR AND THE NUMBER OF
!***  OF FORECAST TASKS IN THE I AND J DIRECTIONS AND
!***  GIVE THOSE TO SETUP_SERVERS WHICH WILL SPLIT THE
!***  COMMUNICATOR BETWEEN FORECAST AND QUILT/WRITE TASKS.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Retrieve Global Communicator"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm             =VM                                &
                     ,mpiCommunicator=MPI_INTRA                         &  !<-- The global communicator
                     ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get INPES/JNPES from Config File"
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
!***  SET UP QUILT/WRITE TASK SPECIFICATIONS.
!***  FIRST RETRIEVE THE TASK AND GROUP COUNTS FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get Write Task/Group Info from Config File"
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
!***  SEGREGATE THE FORECAST TASKS FROM THE QUILT/WRITE TASKS.
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
      NUM_PES_FCST=INPES*JNPES                                           !<-- Number of forecast tasks
      atm_int_state%NUM_PES_FCST=NUM_PES_FCST                            !<-- Save this for the ATM's Run step
!
!-----------------------------------------------------------------------
!***  ALLOCATE AND FILL THE TASK LIST THAT HOLDS THE IDs OF
!***  THE FORECAST TASKS.
!-----------------------------------------------------------------------
!
      ALLOCATE(atm_int_state%PETLIST_FCST(NUM_PES_FCST))                   !<-- Task IDs of the forecast tasks
!
      DO N=0,NUM_PES_FCST-1
        atm_int_state%PETLIST_FCST(N+1)=N                                  !<-- Collect just the forecast task IDs
      ENDDO
!
!-----------------------------------------------------------------------
!***  CREATE DE LAYOUT BASED ON THE I TASKS BY J TASKS SPECIFIED IN
!***  THE CONFIG FILE.
!***  THIS REFERS ONLY TO FORECAST TASKS.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                            !<-- Select only the forecast tasks
        MY_DE_LAYOUT=ESMF_DELayoutCreate(            VM                 &  !<-- The ESMF virtual machine
                                        ,deCountList=(/INPES,JNPES/)    &  !<-- User-specified I-task by J-task layout
                                        ,rc         =RC)
      ENDIF
!
!-----------------------------------------------------------------------
!***  CREATE THE ESMF GRID.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  EXTRACT THE DIMENSIONS OF THE DOMAIN FROM THE CONFIGURE FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get IM,JM,LM from Config File"
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
!***  RETRIEVE THE FORECAST DOMAIN MODE FROM THE CONFIG FILE.
!------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Get GLOBAL/REGIONAL Mode from Config File"
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
        GLOBAL=.true.
        print *,'Initialized as a global run.'
      ELSE
        GLOBAL=.false.
        print *,'Initialized as a regional run.'
      ENDIF
!
!-----------------------------------------------------------------------
!***  IF THIS IS A GLOBAL MODE FORECAST, EXTEND IM AND JM.
!***  THE FIRST DIMENSION OF NCOUNTS IS THE I DIMENSION FOR PARALLELIZATION.
!***  THE SECOND DIMENSION OF NCOUNTS IS THE J DIMENSION.
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
!***  NOW CREATE THE ATM GRIDDED COMPONENT's ESMF GRID
!***  FOR THE NMM's INTEGRATION GRID.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="CORE_SETUP: Create the ESMF Grid"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      GRID_ATM=ESMF_GridCreateShapeTile(regDecomp     =(/INPES,JNPES/)    &  !<-- I x J task layout
                                       ,minIndex      =(/MIN(1),MIN(2)/)  &  !<-- Min indices in I and J
                                       ,maxIndex      =(/MAX(1),MAX(2)/)  &  !<-- Max indices in I and J
                                       ,gridEdgeLWidth=(/0,0/)            &  !<-- Padding, lower edges for noncentered stagger
                                       ,gridEdgeUWidth=(/0,0/)            &  !<-- Padding, upper edges for noncentered stagger
                                       ,name          ="GRID"             &  !<-- Name of the Grid
                                       ,rc            =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CORE)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  GET THE LOCAL ARRAY SIZES FOR THE ATM GRID.
!***  ONLY FORECAST TASKS ARE RELEVANT HERE.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                              !<-- Select only fcst tasks
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="CORE_SETUP: Get EMSF Sizes of Local Subdomains"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_GridGet(grid              =GRID_ATM                   &
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
!***  USING 'computationalCount' FROM ARRAY I1 OBTAINED IN THE
!***  PREVIOUS CALL, GENERATE ALL OF THE LOCAL TASK INDEX LIMITS
!***  FOR ALL FORECAST TASKS.
!***  THE USER, NOT ESMF, DOES THIS WORK.
!-----------------------------------------------------------------------
!
      IF(MYPE<NUM_PES_FCST)THEN                                            !<-- Select only the forecast tasks
        CALL DECOMP(MYPE,INPES,JNPES,NUM_PES_FCST,IM,JM,LM,GLOBAL,I1)
!
        ALLOCATE(atm_int_state%LOCAL_ISTART(0:NUM_PES_FCST-1))
        ALLOCATE(atm_int_state%LOCAL_IEND  (0:NUM_PES_FCST-1))
        ALLOCATE(atm_int_state%LOCAL_JSTART(0:NUM_PES_FCST-1))
        ALLOCATE(atm_int_state%LOCAL_JEND  (0:NUM_PES_FCST-1))
!
        DO N=0,NUM_PES_FCST-1
          atm_int_state%LOCAL_ISTART(N)=LOCAL_ISTART(N)                    !<-- Starting I for all forecasts' subdomains
          atm_int_state%LOCAL_IEND  (N)=LOCAL_IEND  (N)                    !<-- Ending I for all forecasts' subdomains
          atm_int_state%LOCAL_JSTART(N)=LOCAL_JSTART(N)                    !<-- Starting J for all forecasts' subdomains
          atm_int_state%LOCAL_JEND  (N)=LOCAL_JEND  (N)                    !<-- Ending J for all forecasts' subdomains
        ENDDO
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE NMM_SETUP

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      SUBROUTINE GFS_SETUP(gc_atm,grid_atmos)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE CONTAINS NMM-SPECIFIC CODE FOR THE ATM COMPONENT:
!***    (1) SETTING UP DISTRIBUTED MEMORY PARALLELISM IN THE NMM;
!***    (2) CREATING THE ESMF Grid FOR THE ATM COMPONENT;
!***    (3) SHARING LOCAL SUBDOMAIN INDEX LIMITS AMONG TASKS.
!-----------------------------------------------------------------------
!
!
      USE MODULE_INCLUDE
!


!
      type(ESMF_gridcomp),intent(inout) :: gc_atm
      type(ESMF_grid),intent(out)  :: grid_atmos    ! the ESMF grid for the integration attached to



!
!-----------------------------------------------------------------------
!***  SET UP THE ESMF GRID FOR THE NMM AND ESTABLISH THE
!***  DISTRIBUTED MEMORY PARALLELISM.
!-----------------------------------------------------------------------
!
      type(ESMF_config)            :: cf           ! the config object
!
      type(ESMF_DistGrid)          :: DistGrid_atmos
      type(ESMF_VM)                :: vm
!
      integer, dimension(2)         :: ncounts     ! parameter array to set up the
                                                   ! size of the 2-d ESMF grid.
      integer, dimension(2)         :: min,max     ! parameter arrays to set up the
                                                   ! start number and the end number of
                                                   ! the ESMF grid in each dimension.
      real(ESMF_kind_r8),dimension(ESMF_maxgriddim) :: mincoords,maxcoords
      integer,dimension(ESMF_maxgriddim) :: counts
      INTEGER , DIMENSION(2)             :: i1
      INTEGER , DIMENSION(:, :), POINTER :: i2
      integer                      :: num_pes_fcst,num_pes_tot,num_pes,im,jm,lm
      integer                      :: mpi_intra,mpi_intra_b     ! the mpi intra-communicator
      integer                      :: rc,irtn,mype
      integer                      :: RC_RUN
      integer                      :: inpes,jnpes  ! mpi tasks in i and j
      logical                      :: global
      character(50)                :: mode
      rc     =ESMF_success
      RC_RUN=ESMF_success


!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      CALL ESMF_GridCompGet(gridcomp=gc_atm                      &  !<-- The ATM gridded component
                           ,config  =cf                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)

!-----------------------------------------------------------------------
!***  RETRIEVE THE VM (VIRTUAL MACHINE) OF THE ATM GRIDDED COMPONENT.
!***  CALL ESMF_GridCompGet TO RETRIEVE THE VM ANYWHERE YOU NEED IT.
!***  WE NEED VM NOW TO SET UP THE DE LAYOUT.
!-----------------------------------------------------------------------
!
      CALL ESMF_logwrite("retrieve the config object and vm ",          &
                         ESMF_log_info,rc=RC)

!
      CALL ESMF_gridcompget(         gc_atm                             &
                           ,config  =cf                                 &
                           ,vm      =vm                                 &
                           ,rc      =RC)
!
!
!-----------------------------------------------------------------------
!***  set up parameters for mpi communications.
!***  use ESMF utility to get pe identification and total number of pes
!***  (referred to here as nodes.)
!-----------------------------------------------------------------------
!
      CALL ESMF_logwrite("get mype and nodes from vm",ESMF_log_info,rc=RC)
!
      CALL ESMF_vmget(vm                        &  !<-- the virtual machine
                     ,localpet=mype             &  !<-- local pe rank
                     ,petcount=num_pes          &  !<-- total # of tasks
                     ,rc      =RC)
!
      num_pes_tot=num_pes

! Allocate the local index array i2 to store the local size information of the
! ditributed grid.  Information is based per dimension and per De.
!-----------------------------------------------------------------------------
      ALLOCATE(i2(2, num_pes))

!
!***  note: at this point, num_pes is the total number of mpi tasks,
!***        i.e., forecast tasks + quilt tasks.
!
!
!-----------------------------------------------------------------------
!***  establish the task layout including the quilt servers
!***  here in the main gridded component.  get the global
!***  mpi communicator and give it to setup_servers who will
!***  split it between forecast and quilt tasks.
!-----------------------------------------------------------------------
!
      CALL ESMF_vmget(vm                                  &
                     ,mpicommunicator=mpi_intra           &  !<-- the global communicator
                     ,rc             =RC)
!
      CALL mpi_comm_dup(mpi_intra,mpi_intra_b,rc)
!
      CALL ESMF_configgetattribute(cf                      &
                                  ,value =inpes            &  !<-- # of fcst tasks in i direction
                                  ,label ='inpes:'         &
                                  ,rc    =RC)
!
      CALL ESMF_configgetattribute(cf                      &
                                  ,value =jnpes            &  !<-- # of fcst tasks in j direction
                                  ,label ='jnpes:'         &
                                  ,rc    =RC)

      num_pes_fcst=num_pes
!-----------------------------------------------------------------------
!***  create the ESMF grid.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  the first dimension of ncounts is the i dimension for parallelization.
!***  the second dimension of ncounts is the j dimension.
!-----------------------------------------------------------------------
!
      CALL ESMF_configgetattribute(       cf                            &
                                  ,value =im                            &
                                  ,label ='im:'                         &
                                  ,rc    =RC)
!
      CALL ESMF_configgetattribute(       cf                            &
                                  ,value =jm                            &
                                  ,label ='jm:'                         &
                                  ,rc    =RC)
!
!---------------------------------------------------------
!***  if this is a global mode forecast, extend im and jm.
!***  retrieve the mode from the config file.
!---------------------------------------------------------
!
      CALL ESMF_configgetattribute(       cf                            &
                                  ,value =mode                          &
                                  ,label ='global:'                     &
                                  ,rc    =RC)
!
      if(trim(mode)=='.true.')then
        GLOBAL=.true.
        print *,'Run as a global model.'
      else
        GLOBAL=.false.
        print *,'Run as a regional model.'
      endif
!
      if(GLOBAL)then      !<-- global mode horizontal dimensions.
        ncounts(1)=im+2
        ncounts(2)=jm+2
      else                !<-- regional mode horizontal dimensions.
        ncounts(1)=im
        ncounts(2)=jm
      endif
!
!
      max(1)=ncounts(1)
      max(2)=ncounts(2)
!
      min(1)=1
      min(2)=1
!
!
!
!-----------------------------------------------------------------------
!***  now create the main gridded component's ESMF grid.
!-----------------------------------------------------------------------
!
! Create the ESMF DistGrid_atmos.
!--------------------------------
!      CALL ESMF_LogWrite("Create DistGrid_atmos", ESMF_LOG_INFO, rc = rc)

      DistGrid_atmos = ESMF_DistGridCreate(minIndex  = min,              &
                                           maxIndex  = max,              &
                                           regDecomp = (/inpes, jnpes/), &
                                           rc        = rc)


! Create the ESMF grid_atmos based on the created ESMF DistGrid_atmos information.
!---------------------------------------------------------------------------------
      CALL ESMF_logwrite("create grid_atmos",ESMF_log_info,rc=RC)
!
      grid_atmos = ESMF_GridCreate(name     = "grid_atmos",   &
                                   distgrid = DistGrid_atmos, &
                                   rc       = rc)

!
!-----------------------------------------------------------------------
!***  attach the ESMF grid to the main gridded component.
!-----------------------------------------------------------------------
!
!      CALL ESMF_gridcompset(         gc_atm                          &
!                           ,grid    =grid_atmos                      &
!                           ,rc      =RC)
!
!-----------------------------------------------------------------------
!***  get the local array sizes for the main grid.
!***  again, only forecast tasks are relevant here.
!-----------------------------------------------------------------------
!
      if(mype<num_pes_fcst)then
          i2 = 0
          CALL ESMF_DistGridGet(DistGrid_atmos, indexCountPDimPDe = i2, rc = rc)
      endif
!
!-----------------------------------------------------------------------
!***  USING 'computationalCount' FROM ARRAY I1 OBTAINED IN THE
!***  PREVIOUS CALL, GENERATE ALL OF THE LOCAL TASK INDEX LIMITS
!***  FOR ALL FORECAST TASKS.
!***  THE USER, NOT ESMF, DOES THIS WORK.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE GFS_SETUP

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      END MODULE MODULE_CORE_SETUP
!
!-----------------------------------------------------------------------

