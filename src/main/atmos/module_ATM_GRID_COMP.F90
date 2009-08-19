!-----------------------------------------------------------------------
!
      MODULE MODULE_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
!
!***  THIS IS THE ATM (Atmosphere) GRIDDED COMPONENT MODULE.
!***  IT WILL SET UP GRIDDED AND COUPLER SUBCOMPONENTS
!***  AND RUN THEIR INITIALIZE, RUN, AND FINALIZE ROUTINES.
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2007-       Black - Modified from Wei-yu's version
!   2007-11-20  Black/Henderson - Created an ESMF Clock for the
!                                 ATM Component independent of
!                                 the Main Clock.
!   2007-12-11  Black - Generalized for easier use by any dynamics core.
!   2008-08     Colon - Added conditional checks multiple dynamics cores.
!   2008-10-14  Vasic - Added restart Alarm.
!   2009-05-29  Wang  - Added gfs write grid component
!   2009-07-23  Lu    - Added gocart grid component
!
! USAGE: ATM Gridded component parts CALLed from MAIN_ESMF.F90
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE MODULE_INCLUDE
!
      USE MODULE_ATM_INTERNAL_STATE,ONLY: ATM_INTERNAL_STATE            &
                                         ,WRAP_ATM_INTERNAL_STATE
      
!
      USE MODULE_DYNAMICS_GRID_COMP  
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,IHALO,JHALO                         &
                                   ,MPI_COMM_COMP                       &
                                   ,MPI_COMM_INTER_ARRAY                &
                                   ,MYPE_SHARE

!
      USE MODULE_PHYSICS_GRID_COMP,ONLY: PHY_REGISTER
!
      USE MODULE_DYN_PHY_CPL_COMP,ONLY: DYN_PHY_CPL_REGISTER 
!
      USE MODULE_GET_CONFIG_DYN
      USE MODULE_GET_CONFIG_PHY
      USE MODULE_GET_CONFIG_WRITE
      USE MODULE_CONTROL,ONLY: TIMEF
      USE MODULE_DIAGNOSE,ONLY : FIELD_STATS
      USE NEMSIO_MODULE
!
!      use Chem_RegistryMod                                      
!      use m_die, only: die
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
     
      
      USE gfs_dynamics_grid_comp_mod  ,only: gfs_dyn_setservices

      USE gfs_physics_grid_comp_mod   ,only: gfs_phy_setservices

!      USE GOCART_GridCompMod ,only: gocart_setservices => setservices 

      USE atmos_dyn_phy_cpl_comp_mod  ,only: atm_cpl_setservices

!
!-----------------------------------------------------------------------
!***  LIST OTHER MODULES WITH NON_GENERIC ROUTINES USED BY ATM.
!-----------------------------------------------------------------------
!
      USE MODULE_WRITE_ROUTINES ,ONLY: WRITE_INIT,WRITE_ASYNC             !<-- These are routines used only when asynchronous
      USE MODULE_WRITE_GRID_COMP,ONLY: WRITE_SETUP                      & !    quilting is specified by the user in the
                                      ,WRITE_DESTROY                      !    configure file for history output.
!
!jw for GFS
      USE MODULE_WRITE_GRID_COMP_GFS,ONLY: WRITE_SETUP_GFS                  & !    quilting is specified by the user in the
                                          ,WRITE_DESTROY_GFS
      USE MODULE_WRITE_ROUTINES_GFS ,ONLY: WRITE_INIT_GFS                 !<-- These are routines used only when asynchronous

!
      USE MODULE_CLOCKTIMES
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
#include "../../../inc/ESMF_LogMacros.inc"
!-----------------------------------------------------------------------
! 
      PRIVATE
!
      PUBLIC :: ATM_REGISTER

!
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT)       :: DT                                         !<-- The fundamental timestep (sec)
!
      INTEGER(KIND=KINT)    :: MYPE                                       !<-- Each MPI task ID

      INTEGER(KIND=KINT)    :: NPE_PRINT  
!
      CHARACTER(3)          :: CORE                                       !<-- The name of the selected dynamic core
!
      LOGICAL               :: RESTARTED_RUN                              !<-- Original/restarted run logical flag
!
      TYPE(ESMF_Clock),SAVE :: CLOCK_ATM                                  !<-- The ATM Component's ESMF Clock
      TYPE(ESMF_VM),SAVE    :: VM                                         !<-- The ESMF virtual machine.
!
      TYPE(ESMF_TimeInterval),SAVE :: TIMEINTERVAL_CLOCKTIME                 &
                              ,TIMEINTERVAL_HISTORY                          & !<--
                              ,TIMEINTERVAL_RESTART                           !<-- The ESMF Alarm for restart output
!jw
      TYPE(ESMF_TimeInterval),SAVE :: TIMEINTERVAL_gfs_output

!
      TYPE(ESMF_GridComp),SAVE :: GC_GFS_DYN
      TYPE(ESMF_GridComp),SAVE :: GC_GFS_PHY
      TYPE(ESMF_GridComp),SAVE :: GC_GOCART    
      TYPE(ESMF_CplComp), SAVE :: GC_ATM_CPL
      TYPE(ESMF_GridComp),DIMENSION(:),POINTER,SAVE :: WRT_COMPS
!
!      type(Chem_Registry), save :: reg
      character(len=*), parameter ::  myname = 'ut_Registry'
!
      TYPE(ESMF_State),SAVE :: IMP_GFS_DYN,EXP_GFS_DYN
      TYPE(ESMF_State),SAVE :: IMP_GFS_PHY,EXP_GFS_PHY
      TYPE(ESMF_State),SAVE :: IMP_GFS_WRT,EXP_GFS_WRT
!
      INTEGER :: INPES,JNPES                                               ! MPI tasks in i and j
      LOGICAL :: QUILTING                                                  !<-- Logical flag for quilting in
      LOGICAL :: STANDALONE_POST                                           !<-- Logical flag for running standalone post
      INTEGER :: WRITE_GROUP_READY_TO_GO                                   ! the write group to use
!
!      INTEGER(KIND=KINT),PUBLIC :: IM,JM,LM
!      INTEGER(KIND=KINT)        :: NUM_PES
!
!
!      LOGICAL :: ADVECT_TRACERS                                         & !<-- Flag for advecting tracers
!                ,OLD_PASSIVE                                            & !<-- Flag for old passive advection
!                ,OPERATIONAL_PHYSICS                                      !<-- Flag to designate use of operational physics suite
!
!-----------------------------------------------------------------------
!***  FOR DETERMINING CLOCKTIMES OF VARIOUS PIECES OF THE DYNAMICS.
!-----------------------------------------------------------------------
!
      REAL(KIND=KFPT) :: btim,btim0
!
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_REGISTER(ATM_GRID_COMP,RC_REG)
! 
!-----------------------------------------------------------------------
!***  Register the ATM gridded component's initialize, run, and finalize
!***  routines.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                !<-- ATM gridded component
!
      INTEGER,INTENT(OUT)               :: RC_REG                       !<-- Return code for register
!     
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)                :: CF                            !<-- The config object
      INTEGER :: RC
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC=ESMF_SUCCESS       ! Error signal variable
      RC_REG=ESMF_SUCCESS   ! Error signal variable

!-----------------------------------------------------------------------
!***  Register the ATM INITIALIZE subroutine.  Since it is just one 
!***  subroutine, use ESMF_SINGLEPHASE.  The second argument is
!***  a pre-defined subroutine type, such as ESMF_SETINIT, ESMF_SETRUN, 
!***  or ESMF_SETFINAL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Configure Object"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

!
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- ATM gridded component
                                     ,ESMF_SETINIT                      &  !<-- Subroutine type
                                     ,ATM_INITIALIZE                    &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the ATM RUN subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for ATM Run"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- ATM gridded component
                                     ,ESMF_SETRUN                       &  !<-- Subroutine type
                                     ,ATM_RUN                           &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)

 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Register the ATM FINALIZE subroutine.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      MESSAGE_CHECK="Set Entry Point for ATM Finalize"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
      CALL ESMF_GridCompSetEntryPoint(ATM_GRID_COMP                     &  !<-- ATM gridded component
                                     ,ESMF_SETFINAL                     &  !<-- Subroutine type
                                     ,ATM_FINALIZE                      &  !<-- User's subroutine name
                                     ,ESMF_SINGLEPHASE                  &
                                     ,RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ 
!
!-----------------------------------------------------------------------
!***  Check the error signal variable and print out the result.
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
!       WRITE(0,*)' ATM_SET_SERVICES SUCCEEDED'
      ELSE
        WRITE(0,*)' ATM_SET_SERVICES FAILED  RC_REG=',RC_REG
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_REGISTER
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_INITIALIZE(ATM_GRID_COMP                           &
                               ,IMP_STATE,EXP_STATE                     &
                               ,CLOCK_MAIN                              &
                               ,RC_INIT)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE SETS UP FUNDAMENTAL ASPECTS OF THE MODEL RUN.
!-----------------------------------------------------------------------
!
      USE MODULE_NMM_CORE_SETUP
      USE MODULE_GFS_CORE_SETUP
      use module_gfs_mpi_def, only: PETLIST_FCST,WRITE_GROUPS
!
!-----------------------------------------------------------------------
!***  ARGUMENT VARIABLES.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                !<-- The ATM gridded component
      TYPE(ESMF_State)   ,INTENT(INOUT) :: IMP_STATE,EXP_STATE          !<-- The ATM component's import/export states
      TYPE(ESMF_Clock)   ,INTENT(INOUT) :: CLOCK_MAIN                   !<-- The main program's ESMF Clock!
      INTEGER,OPTIONAL   ,INTENT(OUT)   :: RC_INIT                      !<-- Return code for Initialize step
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      CHARACTER(ESMF_MAXSTR)           :: INFILE
      TYPE(NEMSIO_GFILE) :: GFILE
      TYPE(WRAP_ATM_INTERNAL_STATE)    :: WRAP                          !<-- The F90 wrap of the ATM internal state
      TYPE(ATM_INTERNAL_STATE),POINTER :: ATM_INT_STATE                 !<-- The ATM internal state pointer
!
      TYPE(ESMF_Config)                :: CF                            !<-- The config object
!
      TYPE(ESMF_Grid)                  :: GRID_ATM                      !<-- The ESMF GRID for the integration attached to
                                                                        !    the ATM gridded component.
      TYPE(ESMF_Grid)                  :: GRID_DYN                      !<-- The ESMF GRID for the integration attached to
                                                                        !    the dynamics gridded component.
      TYPE(ESMF_Grid)                  :: GRID_PHY                      !<-- The ESMF GRID for the integration attached to
                                                                        !    the physics gridded component.
!
      INTEGER                          :: TIMESTEP_SEC_WHOLE            !<-- Integer part of timestep in seconds
      INTEGER                          :: TIMESTEP_SEC_NUMERATOR        !<-- Numerator of timestep fraction
      INTEGER                          :: TIMESTEP_SEC_DENOMINATOR      !<-- Denominator of timestep fraction
      TYPE(ESMF_TimeInterval)          :: TIMESTEP                      !<-- The ESMF timestep (s)
      TYPE(ESMF_TimeInterval)          :: RUNDURATION                   !<-- The ESMF simulation length (h)
      TYPE(ESMF_Time)                  :: CURRTIME                      !<-- The ESMF current time.
      TYPE(ESMF_Time)                  :: STARTTIME                     !<-- The ESMF start time.
      TYPE(ESMF_TimeInterval)          :: BCK_DURATION
      TYPE(ESMF_Grid)  :: grid_gfs_dyn     ! the ESMF grid for the integration attached to
                                                   ! the dynamics gridded component.
      TYPE(ESMF_Grid)  :: grid_gfs_phy     ! the ESMF grid for the integration attached to
!
      TYPE(ESMF_Grid)                  :: grid_atmos    ! the esmf grid for the integration attached to

      INTEGER                          :: NHOURS_CLOCKTIME            & !<-- Hours between clocktime prints
                                         ,NHOURS_HISTORY              & !<-- Hours between history output
                                         ,NHOURS_RESTART                !<-- Hours between restart output
!
      INTEGER                 :: FILTER_METHOD, HALF_BCK_DURATION

      INTEGER                          :: NFCST, NTSD
      INTEGER(ESMF_KIND_I8)            :: NTSD_START                    !<-- Timestep count (>0 for restarted runs
!
!
      LOGICAL                          :: DIST_MEM,PHYSICS_ON           !<-- Logical flag for distributed
      LOGICAL                          :: NEMSIO_INPUT
!
      LOGICAL                          :: OPENED
!
      INTEGER                          :: I,IERR,RC,IRTN,J,K,N,NN
      INTEGER                          :: RC_FINAL
      CHARACTER(50)                    :: MODE
!
      INTEGER                          :: IYEAR_FCST     &              !<-- Current year from restart file
                                         ,IMONTH_FCST    &              !<-- Current month from restart file
                                         ,IDAY_FCST      &              !<-- Current day from restart file
                                         ,IHOUR_FCST     &              !<-- Current hour from restart file
                                         ,IMINUTE_FCST   &              !<-- Current minute from restart file
                                         ,ISECOND_FCST                  !<-- Current second from restart file
!
      REAL                             :: SECOND_FCST                   !<-- Current second from restart file
!
      type(ESMF_DistGrid)          :: DistGrid_atmos

!
      integer, dimension(2)         :: ncounts     ! parameter array to set up the
                                                   ! size of the 2-d ESMF grid.
      integer, dimension(2)         :: min,max     ! parameter arrays to set up the
                                                   ! start number and the end number of
                                                   ! the ESMF grid in each dimension.
      integer,dimension(ESMF_maxgriddim) :: counts
      INTEGER , DIMENSION(2)             :: i1
      INTEGER , DIMENSION(:, :), POINTER :: i2
      integer                      :: im,jm,lm                  ! full grid dimensions
      integer                      :: num_pes,num_pes_fcst,num_pes_tot
      integer                      :: mpi_intra,mpi_intra_b     ! the mpi intra-communicator!
      logical                      :: global
      INTEGER                      :: NFHOUT     
      INTEGER :: FCSTDATE(7)
!
      RC     =ESMF_success


!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  INITIALIZE TIMING VARIABLES.
!-----------------------------------------------------------------------
!
      btim0=timef()
      total_tim=0.
      totalsum_tim=0.
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Config Object from ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      write(0,*)'in atm_init, after get cf'
!
!-----------------------------------------------------------------------
!***  EXTRACT THE DYNAMIC CORE NAME.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Dynamic Core Name"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      write(0,*)'in atm_init, core=',core
!-----------------------------------------------------------------------
!***  ALLOCATE THE ATM COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
      IF (CORE=='nmm') THEN

      ALLOCATE(ATM_INT_STATE,stat=RC)
!
      wrap%ATM_INT_STATE=>ATM_INT_STATE
!
!-----------------------------------------------------------------------
!***  ATTACH THE ATM INTERNAL STATE TO THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach ATM Internal State to Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSetInternalState(ATM_GRID_COMP                  &  !<-- The ATM gridded component
                                        ,WRAP                           &  !<-- Pointer to the ATM internal state
                                        ,RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!jw
      ENDIF
!-----------------------------------------------------------------------
!***  RETRIEVE THE VM (VIRTUAL MACHINE) OF THE ATM GRIDDED COMPONENT.
!***  CALL ESMF_GridCompGet TO RETRIEVE THE VM ANYWHERE YOU NEED IT.
!***  WE NEED VM NOW TO OBTAIN THE MPI TASK IDs.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve VM from ATM Gridded Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,vm      =VM                                 &  !<-- Get the Virtual Machine from the ATM component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Obtain MPI Task IDs from VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGet(vm      =VM                                       &  !<-- The virtual machine
                     ,localpet=MYPE                                     &  !<-- Each MPI task ID
                     ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!***  WILL THE WRITE COMPONENTS WITH ASYNCHRONOUS QUILTING BE USED?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Quilting Flag from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =QUILTING                      &  !<-- The variable filled (flag
                                  ,label ='quilting:'                   &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!jw
      IF (CORE=='nmm') THEN
!
      atm_int_state%QUILTING=QUILTING                                      !<-- Save this for the ATM's Run step
!
!-----------------------------------------------------------------------
!***  EXTRACT FUNDAMENTAL TIMESTEP INFORMATION FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Timestep Information from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =TIMESTEP_SEC_WHOLE            &  !<-- The variable filled (integer part of timestep (sec))
                                  ,label ='dt_int:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =TIMESTEP_SEC_NUMERATOR        &  !<-- The variable filled (numerator of timestep fraction)
                                  ,label ='dt_num:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =TIMESTEP_SEC_DENOMINATOR      &  !<-- The variable filled (denominator of timestep fraction)
                                  ,label ='dt_den:'                     &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  ESTABLISH THE TIMESTEP.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set up Time Step Interval"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP                   &  !<-- The model's fundamental timestep (sec) (ESMF)
                               ,s           =TIMESTEP_SEC_WHOLE         &
                               ,sn          =TIMESTEP_SEC_NUMERATOR     &
                               ,sd          =TIMESTEP_SEC_DENOMINATOR   &
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DT=TIMESTEP_SEC_WHOLE+                                            &  !<-- The model's fundamental timestep (sec) (REAL)
         REAL(TIMESTEP_SEC_NUMERATOR)/REAL(TIMESTEP_SEC_DENOMINATOR)
!
!-----------------------------------------------------------------------
!***  OBTAIN THE FORECAST START/END TIMES FROM THE MAIN CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Start Time and Fcst Duration from Main Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock      =CLOCK_MAIN                         &  !<-- The Main ESMF Clock
                        ,startTime  =STARTTIME                          &  !<-- The simulation start time
                        ,runDuration=RUNDURATION                        &  !<-- The simulation run duration
                        ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CURRTIME=STARTTIME                                                   !<-- Set the current forecast time
!
!-----------------------------------------------------------------------
!***  EXTRACT RESTART LOGICAL FROM CF;
!***  IF RESTARTED RUN:
!***     READ NTSD FROM restart_file
!***     SET CURRTIME TO TIME WHEN RESTART FILE WAS WRITTEN
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract NEMSIO FLAG Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =NEMSIO_INPUT                 &  !<-- The variable filled (logical restart or cold start)
                                  ,label ='nemsio_input:'                    &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!***  EXTRACT RESTART LOGICAL FROM CF;
!***  IF RESTARTED RUN:
!***     READ NTSD FROM restart_file
!***     SET CURRTIME TO TIME WHEN RESTART FILE WAS WRITTEN
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Restart Logical from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =RESTARTED_RUN                 &  !<-- The variable filled (logical restart or cold start)
                                  ,label ='restart:'                    &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      NTSD_START=0
!
      IF(RESTARTED_RUN) THEN                                               !<-- If this is a restarted run, set the current time
!
        if (NEMSIO_INPUT) then
!
          INFILE="restart_file_nemsio"
!
!----------------------------------------------------------------------
!*** read restart data
!-----------------------------------------------------------------------
!
          CALL NEMSIO_INIT()
!
          CALL NEMSIO_OPEN(gfile,INFILE,'read',iret=irtn)
!
!-----------------------------------------------------------------------
!              READ FROM RESTART FILE: INTEGER SCALARS
!-----------------------------------------------------------------------
!
          CALL NEMSIO_GETHEADVAR(gfile,'FCSTDATE',FCSTDATE,iret=irtn)
          IYEAR_FCST=FCSTDATE(1)
          IMONTH_FCST=FCSTDATE(2)
          IDAY_FCST=FCSTDATE(3)
          IHOUR_FCST=FCSTDATE(4)
          IMINUTE_FCST=FCSTDATE(5)
          SECOND_FCST=0.
          if(FCSTDATE(7)/=0) SECOND_FCST=FCSTDATE(6)/(FCSTDATE(7)*1.)
          CALL NEMSIO_GETHEADVAR(gfile,'NTIMESTEP',NTSD,iret=irtn)

          call nemsio_close(gfile,iret=ierr)
!
        else
!
          INFILE="restart_file"
!

          select_unit: DO N=51,59
            INQUIRE(N,OPENED=OPENED)
            IF(.NOT.OPENED)THEN
              NFCST=N
              EXIT select_unit
            ENDIF
          ENDDO select_unit
!
          OPEN(unit=NFCST,file=INFILE,status='old',form='unformatted')
!
          READ(NFCST) IYEAR_FCST                                             !<-- Read time form restart file
          READ(NFCST) IMONTH_FCST                                            !
          READ(NFCST) IDAY_FCST                                              !
          READ(NFCST) IHOUR_FCST                                             !
          READ(NFCST) IMINUTE_FCST                                           !
          READ(NFCST) SECOND_FCST                                            !<--
!
          READ(NFCST) NTSD                                                   !<-- Read timestep from restart file
!
          CLOSE(NFCST)

        endif
!
       ISECOND_FCST=NINT(SECOND_FCST)                                      !<-- ESMF clock needs integer seconds
        NTSD_START=NTSD
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="RESTART: Set the Forecast Current Time"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_TimeSet(time=CURRTIME                                   &  !<-- Current time of the forecast (ESMF)
                         ,yy  =IYEAR_FCST                                 &  !<-- Year from restart file
                         ,mm  =IMONTH_FCST                                &  !<-- Month from restart file
                         ,dd  =IDAY_FCST                                  &  !<-- Day from restart file
                         ,h   =IHOUR_FCST                                 &  !<-- Hour from restart file
                         ,m   =IMINUTE_FCST                               &  !<-- Minute from restart file
                         ,s   =ISECOND_FCST                               &  !<-- Second from restart file
                         ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
       ENDIF
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  WITH DATA FROM ABOVE, CREATE THE LOCAL ESMF CLOCK
!***  TO CONTROL THE TIMESTEPPING IN THE ATM COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Local Clock for the ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CLOCK_ATM=ESMF_ClockCreate(name       ='CLOCK_ATM'                &  !<-- The ATM Clock's name
                                ,timestep   =TIMESTEP                   &  !<-- The fundamental timestep in this component
                                ,starttime  =STARTTIME                  &  !<-- Start time of simulation
                                ,runduration=RUNDURATION                &  !<-- Duration of simulation
                                ,rc         =RC)
!
      CALL ESMF_ClockSet(clock       =CLOCK_ATM                         &  !<-- The ATM Component's Clock
                        ,currtime    =CURRTIME                          &  !<-- Current time of simulation
                        ,advanceCount=NTSD_START                        &  !<-- Timestep at this current time
                        ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =FILTER_METHOD                 &
                                  ,label ='filter_method:'              &
                                  ,rc    =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! 
      IF (.NOT. RESTARTED_RUN) THEN
      IF (FILTER_METHOD .eq. 2) THEN

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      MESSAGE_CHECK="DDFI employed so obtaining backward integration duration"
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value = HALF_BCK_DURATION             &
                                  ,label ='nsecs_bckddfi:'            &
                                  ,rc    =RC)
      CALL ESMF_TimeIntervalSet(timeinterval= BCK_DURATION         &  
                               ,s           =  2*HALF_BCK_DURATION    &  
                               ,rc          =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ELSE IF (FILTER_METHOD .eq. 3) THEN

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      MESSAGE_CHECK="TDFI employed so obtaining backward integration duration"
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value = HALF_BCK_DURATION             &
                                  ,label ='nsecs_bcktdfi:'            &
                                  ,rc    =RC)
      CALL ESMF_TimeIntervalSet(timeinterval= BCK_DURATION         &
                               ,s           =  2*HALF_BCK_DURATION    &
                               ,rc          =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      ELSE

      HALF_BCK_DURATION=0

      CALL ESMF_TimeIntervalSet(timeinterval= BCK_DURATION         &  
                               ,s           = HALF_BCK_DURATION     &  
                               ,rc          =RC)
      ENDIF



! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Set the Forecast Length"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

       STARTTIME=STARTTIME+BCK_DURATION
       CURRTIME=STARTTIME

       CALL ESMF_ClockSet(clock   =CLOCK_ATM                &  
                        ,currtime=CURRTIME                  &  
                        ,starttime=STARTTIME                &
                        ,rc      =RC)
      ENDIF

!-----------------------------------------------------------------------
!***  GET THE HISTORY OUTPUT INTERVAL (HOURS) FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Obtain History Interval from the Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NHOURS_HISTORY                &  !<-- Fill this variable
                                  ,label ='nhours_history:'             &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CREATE THE HISTORY FILE OUTPUT ALARM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the History Output Alarm"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMEINTERVAL_HISTORY           &  !<-- Time interval between
                               ,h           =NHOURS_HISTORY                 &  !<-- Hours between history
                               ,rc          =RC)
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  GET THE RESTART OUTPUT INTERVAL (HOURS) FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Obtain Restart Interval from the Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NHOURS_RESTART                &  !<-- Fill this variable
                                  ,label ='nhours_restart:'             &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CREATE THE RESTART FILE OUTPUT ALARM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Restart Output Alarm"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMEINTERVAL_RESTART           &  !<-- Time interval between restart writes (h) (ESMF)
                               ,h           =NHOURS_RESTART                 &  !<-- Hours between restart writes (REAL)
                               ,rc          =RC)
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  CREATE THE ALARM FOR PRINTING CLOCKTIMES USED BY MODEL PIECES.
!***  READ IN FORECAST TIME INTERVAL FOR CLOCKTIME OUTPUT AS WELL AS
!***  THE SELECTED TASK ID THAT WILL PROVIDE THE CLOCKTIMES.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Read Fcst Interval for Clocktime Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NHOURS_CLOCKTIME              &  !<-- Fill this variable
                                  ,label ='nhours_clocktime:'           &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
       
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Read MPI Task ID That Provides Clocktime Output"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =NPE_PRINT                     &  !<-- Fill this variable
                                  ,label ='npe_print:'                  &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create ESMF Clocktime Output Interval"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMEINTERVAL_CLOCKTIME         &  !<-- Time interval between
                               ,h           =NHOURS_CLOCKTIME               &  !<-- Hours between clocktime writes (REAL)
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
      ENDIF

      
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  MODEL-SPECIFIC ROUTINES MUST BE INVOKED IN ORDER TO ESTABLISH
!***  THE ESMF Grid.  THE DIFFERENT INTEGRATION GRIDS NECESSITATE
!***  DIFFERENT WAYS OF SETTING UP BOTH THE PARALLELISM FOR
!***  DISTRIBUTED MEMORY RUNS AND THE ESMF Grid ITSELF.
!***  WHEN THE PARALLELISM IS CONSTRUCTED, THE LOCAL DOMAIN LIMITS
!***  NEED TO BE INSERTED INTO THE ATM COMPONENT'S INTERNAL STATE
!***  IF QUILTING IS TO BE USED.  SEE 'IF(QUILTING)THEN' BELOW.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF(CORE=='nmm')THEN
        CALL NMM_SETUP(ATM_GRID_COMP,ATM_INT_STATE,GRID_ATM)
!-----------------------------------------------------------------------
!***  ATTACH THE CORE-SPECIFIC ESMF GRID TO THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the ESMF Grid to the ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompSet(gridcomp=ATM_GRID_COMP                      & !<-- The ATM gridded component
                           ,grid    =GRID_ATM                           & !<-- Attach the ESMF grid to the ATM component
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
!-----------------------------------------------------------------------

      ELSE IF (CORE=='gfs') THEN
        CALL GFS_SETUP(ATM_GRID_COMP,grid_atmos)
!jws
!jw        if(QUILTING ) then
        CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The configure object
                                  ,value =nfhout                          &  !<-- Fill this variable
                                  ,label ='nfhout:'                       &  !<-- Give the variable this label's value from the config file
                                  ,rc    =RC)
!
        CALL ESMF_TimeIntervalSet(timeinterval=TIMEINTERVAL_GFS_OUTPUT    &  !<-- Time interval between
                               ,h           =nfhout                       &  !<-- Hours between history
                               ,rc          =RC)
!jw        endif
!jwe

!  Read and print Chem_Registry (Sarah Lu)
!  ------------------------
!        print *, 'Read Chem_Registry'
!        reg = Chem_RegistryCreate ( ierr )
!        if ( ierr /= 0 ) call die ( myname, 'cannot create registry' )
!        CALL Chem_RegistryPrint ( reg )
!        call Chem_RegistryDestroy ( reg, ierr )
!        if ( ierr /= 0 ) call die ( myname, 'cannot destroy registry' )


      ENDIF
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CREATE THE DYNAMICS GRIDDED SUBCOMPONENT.
!***  REGISTER THE INITIALIZE, RUN, AND FINALIZE STEPS FOR IT.
!***  SINCE THERE IS ONLY A SINGLE INTEGRATION GRID, GIVE THE
!***  DYNAMICS THE ATM COMPONENT'S GRID.
!***  NOTE THAT THIS SUBCOMPONENT IS PART OF THE ATM COMPONENT'S
!***  INTERNAL STATE.  THIS WILL BE CONVENIENT IF WE NEED TO REACH
!***  THE Dynamics COMPONENT VIA THE ATM COMPONENT SUCH AS HAPPENS
!***  WHEN WRITE COMPONENTS ARE ESTABLISHED.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-------------------------------
!***  Create Dynamics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      atm_int_state%DYN_GRID_COMP=ESMF_GridCompCreate(                  &
                                  name      ="Dynamics component"       &  !<-- Name of the new Dynamics gridded component
                                 ,configFile='configure_file'           &  !<-- Attach this configure file to the component
                                 ,petList   =atm_int_state%PETLIST_FCST &  !<-- The forecast task IDs
                                 ,rc        =RC)
      ELSE IF (CORE=='gfs') THEN
      gc_gfs_dyn=esmf_gridcompcreate( name      ="dynamics component"   &
                                     ,configFile='dyn_namelist.rc'      &
!jw
                                     ,petList=PETLIST_FCST              &
                                     ,rc        =RC)
      write(0,*)'in atm_init, after dyn comp created'
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------
!***  Register the Init, Run, and Finalize steps
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register Dynamics Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
       
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompSetServices(atm_int_state%DYN_GRID_COMP         &  ! <-- The Dynamics gridded component
                                   ,DYN_REGISTER                        &  ! <-- The user's subroutineName for Register
                                   ,RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompsetservices(gc_gfs_dyn                   &
                                   ,gfs_dyn_setservices             &
                                   ,RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!----------------------------------------------------
!***  Attach the ESMF Grid to the Dynamics component
!----------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach ESMF Grid to the Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      
      IF (CORE=='nmm') THEN
      GRID_DYN=GRID_ATM
      CALL ESMF_GridCompSet(gridcomp=atm_int_state%DYN_GRID_COMP        &  !<-- The Dynamics component
                           ,grid    =GRID_DYN                           &  !<-- The Dynamics ESMF grid
                           ,rc      =RC)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create empty Import and Export states for the Dynamics component
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Dynamics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      atm_int_state%IMP_STATE_DYN=ESMF_StateCreate(stateName="Dynamics Import" &  !<-- The Dynamics import state name
                                                  ,statetype=ESMF_STATE_IMPORT &
                                                  ,rc       =RC)
!
      atm_int_state%EXP_STATE_DYN=ESMF_StateCreate(stateName="Dynamics Export" &  !<-- The Dynamics export state name
                                                  ,statetype=ESMF_STATE_EXPORT &
                                                  ,rc       =RC)
      ELSE IF (CORE=='gfs') THEN
      imp_gfs_dyn=esmf_statecreate(statename="dynamics import"        &
                                    ,statetype=esmf_state_import      &
                                    ,rc       =RC)
!
      exp_gfs_dyn=esmf_statecreate(statename="dynamics export"        &
                                    ,statetype=esmf_state_export      &
                                    ,rc       =RC)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------------------------
!***  CHECK CONDITION TO RUN ADIABATIC CASE (NO PHYSICS) 
!---------------------------------------------------------
!
      call ESMF_ConfigGetAttribute(       CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =rc)
      IF (CORE=='nmm') THEN
      if(trim(mode)=='true')then
        PHYSICS_ON=.false.
        print *,' initialize without physics coupling '
      else
        PHYSICS_ON=.true.
        print *,' initialize with physics coupling '
      endif
      ELSE IF (CORE=='gfs') THEN
      if(trim(mode)=='.true.')then
        PHYSICS_ON=.false.
        print *,' initialize without physics coupling '
      else
        PHYSICS_ON=.true.
        print *,' initialize with physics coupling '
      endif
      ENDIF
!
      IF (PHYSICS_ON) THEN
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CREATE THE PHYSICS GRIDDED SUBCOMPONENT.
!***  REGISTER THE INITIALIZE, RUN, AND FINALIZE STEPS FOR IT.
!***  SINCE THERE IS ONLY A SINGLE INTEGRATION GRID, GIVE THE
!***  PHYSICS THE ATM COMPONENT'S GRID.
!***  NOTE THAT THIS SUBCOMPONENT IS PART OF THE ATM COMPONENT'S
!***  INTERNAL STATE.  THIS WILL BE CONVENIENT IF WE NEED TO REACH
!***  THE Physics COMPONENT VIA THE ATM COMPONENT SUCH AS HAPPENS
!***  WHEN WRITE COMPONENTS ARE ESTABLISHED.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-------------------------------
!***  Create Physics component
!-------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      atm_int_state%PHY_GRID_COMP=ESMF_GridCompCreate(                  &
                                  name      ="Physics component"        &  !<-- Name of the new Physics gridded component
                                 ,configFile='configure_file'           &  !<-- Attach this configure file to the component
                                 ,petList   =atm_int_state%PETLIST_FCST &  !<-- The forecast task IDs
                                 ,rc        =RC)
      ELSE IF (CORE=='gfs') THEN
      gc_gfs_phy=esmf_gridcompcreate( name      ="physics component"    &
                                     ,configfile='phy_namelist.rc'      &
!jw
                                     ,petList=PETLIST_FCST              &
                                     ,rc        =rc)
      write(0,*)'in atm_init, after phys comp created,petlist_fcst=',petlist_fcst
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------
!***  Register the Init, Run, and Finalize steps
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register Physics Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompSetServices(atm_int_state%PHY_GRID_COMP         &  ! <-- The Physics gridded component
                                   ,PHY_REGISTER                        &  ! <-- The user's subroutineName
                                   ,RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompsetservices(gc_gfs_phy                   &
                                   ,gfs_phy_setservices             &
                                   ,rc)

! setservice for gocart
!      call esmf_gridcompsetservices(gc_gocart               &  ! <-- The GOCART gridded component
!                                   ,gocart_setservices      &  ! <-- The user's subroutineName
!                                   ,rc)

!
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------
!***  Attach the ESMF Grid to the Physics
!-----------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Attach the ESMF Grid to the Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      GRID_PHY=GRID_ATM
      CALL ESMF_GridCompSet(gridcomp=atm_int_state%PHY_GRID_COMP        &  !<-- The Physics component
                           ,grid    =GRID_PHY                           &  !<-- The Physics ESMF grid
                           ,rc      =RC)
!      ELSE IF (CORE=='gfs') THEN
!      grid_gfs_phy=grid_atmos
!      call esmf_gridcompset(gridcomp= gc_gfs_phy                         &
!                            ,grid    =grid_gfs_phy                       &
!                            ,rc	     =RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ENDIF !PHYSICS_ON
!
!------------------------------------------------------------------------
!***  Create empty Import and Export states for the Physics subcomponent
!------------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create Empty Import/Export States for Physics"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      atm_int_state%IMP_STATE_PHY=ESMF_StateCreate(stateName="Physics Import"  &  !<-- The Physics import
                                                  ,statetype=ESMF_STATE_IMPORT &
                                                  ,rc       =RC)
!
 
      atm_int_state%EXP_STATE_PHY=ESMF_StateCreate(stateName="Physics Export"  &  !<-- The Physics export
                                                  ,statetype=ESMF_STATE_EXPORT &
                                                  ,rc       =RC)
      
      ELSE IF (CORE=='gfs') THEN
      imp_gfs_phy=esmf_statecreate(statename="physics import"         &
                                    ,statetype=esmf_state_import      &
                                    ,rc       =rc)
!
      exp_gfs_phy=esmf_statecreate(statename="physics export"         &
                                    ,statetype=esmf_state_export      &
                                    ,rc       =rc)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  CREATE THE DYNAMICS-PHYSICS COUPLER SUBCOMPONENT.
!***  REGISTER THE INITIALIZE, RUN, AND FINALIZE STEPS FOR IT.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!----------------------------
!***  Create Dyn-Phy Coupler
!----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Dynamics-Physics Coupler Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      atm_int_state%COUPLER_DYN_PHY_COMP=ESMF_CplCompCreate                   &
                                         (name   ="Dyn-Phy coupler component" &
                                         ,petList=atm_int_state%PETLIST_FCST  &  !<-- The forecast task IDs
                                         ,rc     =RC)
      ELSE IF (CORE=='gfs') THEN
      gc_atm_cpl=esmf_cplcompcreate(name="coupler component"            &
!jw
                                     ,petList=PETLIST_FCST              &
                                     ,rc  =RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!------------------------------------------------
!***  Register the Init, Run, and Finalize steps
!------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register the Dyn-Phy Coupler's Init, Run, Finalize"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_CplCompSetServices(atm_int_state%COUPLER_DYN_PHY_COMP   &  ! <-- The Dyn-Phy coupler component
                                  ,DYN_PHY_CPL_REGISTER                 &  ! <-- The user's subroutineName
                                  ,RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_cplcompsetservices(gc_atm_cpl                &
                                  ,atm_cpl_setservices              &
                                  ,RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  IF QUILTING WAS SELECTED FOR THE GENERATION OF OUTPUT,
!***  SETUP THE WRITE COMPONENT(S).
!***  THIS MUST BE DONE PRIOR TO EXECUTING THE INITIALIZE STEPS
!***  OF THE DYNAMICS AND PHYSICS COMPONENTS BECAUSE THE WRITE
!***  COMPONENTS' IMPORT STATES ARE REQUIRED BY THOSE STEPS
!***  WHEN 'quilting' IS ACTIVE AND WRITE_SETUP IS THE ROUTINE
!***  IN WHICH THE WRITE COMPONENTS' IMPORT STATES ARE INSERTED
!***  INTO THE DYNAMICS AND PHYSICS EXPORT STATES WHICH ARE
!***  IN TURN PART OF THE ATM INTERNAL STATE.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      IF(QUILTING)THEN
        IF (CORE=='nmm') THEN
          CALL WRITE_SETUP(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
        ELSEIF (CORE=='gfs') THEN
          write(0,*)'before write_setup_gfs, allocate,write_groups=',write_groups
          allocate(WRT_COMPS(WRITE_GROUPS))
          write(0,*)'before write_setup_gfs call'
          CALL WRITE_SETUP_GFS(ATM_GRID_COMP,WRT_COMPS,          &
            EXP_GFS_DYN,EXP_GFS_PHY,IMP_GFS_WRT, EXP_GFS_WRT)
      write(0,*)'in atm_init, after writesetup'
        ENDIF
      ENDIF
      
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEPS FOR THE GRIDDED SUBCOMPONENTS.
!***  THESE ARE THE INITIALIZE SUBROUTINES SPECIFIED IN THE
!***  REGISTER ROUTINES CALLED IN ESMF_GridCompSetServices ABOVE.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!--------------
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompInitialize(gridcomp   =atm_int_state%DYN_GRID_COMP  &  !<-- The dynamics gridded component
                                  ,importState=atm_int_state%IMP_STATE_DYN  &  !<-- The dynamics import state
                                  ,exportState=atm_int_state%EXP_STATE_DYN  &  !<-- The dynamics export state
                                  ,clock      =CLOCK_ATM                    &  !<-- The ATM clock
                                  ,phase      =ESMF_SINGLEPHASE             &
                                  ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call ESMF_GridCompInitialize(gridcomp   =gc_gfs_dyn            &
                                  ,importstate=imp_gfs_dyn           &
                                  ,exportstate=exp_gfs_dyn           &
                                  ,clock      =CLOCK_MAIN            &
                                  ,phase      =esmf_singlephase      &
                                  ,rc         =RC)
      grid_gfs_dyn=grid_atmos
      call esmf_gridcompset(         gc_gfs_dyn      &
                           ,grid    =grid_gfs_dyn    &
                           ,rc      =rc)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (PHYSICS_ON) THEN
!-------------
!***  PHYSICS
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompInitialize(gridcomp   =atm_int_state%PHY_GRID_COMP  &  !<-- The physics gridded component
                                  ,importState=atm_int_state%IMP_STATE_PHY  &  !<-- The physics import state
                                  ,exportState=atm_int_state%EXP_STATE_PHY  &  !<-- The physics export state
                                  ,clock      =CLOCK_ATM                    &  !<-- The ATM clock
                                  ,phase      =ESMF_SINGLEPHASE             &
                                  ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompinitialize(gridcomp   =gc_gfs_phy            &
                                  ,importstate=imp_gfs_phy           &
                                  ,exportstate=exp_gfs_phy           &
                                  ,clock      =CLOCK_MAIN            &
                                  ,phase      =esmf_singlephase      &
                                  ,rc         =RC)
      grid_gfs_phy=grid_atmos
      call esmf_gridcompset(         gc_gfs_phy      &
                           ,grid    =grid_gfs_phy    &
                           ,rc      =RC)

      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ENDIF !PHYSICS_ON
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  INITIALIZE THE DYN-PHY COUPLER SUBCOMPONENT.
!***  THE CHOICE OF IMPORT AND EXPORT STATE DOES NOT MATTER
!***  FOR THE INITIALIZE STEP.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Initialize Dyn-Phy Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_CplCompInitialize(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &  !<-- The dyn_phy coupler component
                                 ,importState=atm_int_state%EXP_STATE_DYN        &  !<-- The dyn-phy coupler import state
                                 ,exportState=atm_int_state%IMP_STATE_PHY        &  !<-- The dyn-phy coupler export state
                                 ,clock      =CLOCK_ATM                          &  !<-- The ATM Clock
                                 ,rc         =RC)
      ELSE If (CORE=='gfs') THEN
      call esmf_cplcompinitialize(cplcomp    =gc_atm_cpl              &
                                 ,importstate=exp_gfs_dyn             &
                                 ,exportstate=imp_gfs_phy             &
                                 ,clock      =CLOCK_MAIN              &
                                 ,rc         =RC)
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  IF QUILTING WAS SELECTED FOR THE GENERATION OF OUTPUT,
!***  EXECUTE THE INITIALIZE STEP OF THE WRITE COMPONENT(S).
!***  THIS MUST BE DONE AFTER THE INITIALIZATION OF THE
!***  DYNAMICS AND PHYSICS COMPONENTS BECAUSE THOSE COMPONENTS'
!***  INTERNAL STATES CONTAIN THE HISTORY OUTPUT VARIABLES.
!***  POINTERS TO THOSE VARIABLES ARE SET DURING THEIR INITIALIZE
!***  STEPS AND ARE THEN LOADED INTO THE WRITE COMPONENTS'
!***  IMPORT STATES WHICH THEMSELVES RESIDE IN THE DYNAMICS/PHYSICS
!***  EXPORT STATES.
!-----------------------------------------------------------------------
!
      IF(QUILTING)THEN
        IF (CORE=='nmm') THEN
          CALL WRITE_INIT(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
        ELSEIF(CORE=='gfs') THEN
          CALL WRITE_INIT_GFS(ATM_GRID_COMP,WRT_COMPS,IMP_GFS_WRT,          &
            EXP_GFS_WRT,CLOCK_MAIN,WRITE_GROUP_READY_TO_GO)
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!***  WRITE THE FINAL ERROR SIGNAL.
!-----------------------------------------------------------------------
!
      IF(RC_INIT==ESMF_SUCCESS)THEN
!       WRITE(0,*)'ATM INITIALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'ATM INITIALIZE STEP FAILED RC_INIT=',RC_INIT
      ENDIF
!
!-----------------------------------------------------------------------
!
      total_tim=timef()-btim0
!
      if(mype==0)write(0,*)' atm_init_tim=',total_tim*1.e-3
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_INITIALIZE
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_RUN(ATM_GRID_COMP                                  &
                        ,IMP_STATE,EXP_STATE                            &
                        ,CLOCK_MAIN                                     &
                        ,RC_RUN)
!
!-----------------------------------------------------------------------
!***  RUN THE ATM (Atmosphere) GRIDDED COMPONENT.
!-----------------------------------------------------------------------
      USE MODULE_NMM_FWD_INTEGRATE
      USE MODULE_NMM_BCK_INTEGRATE
      USE MODULE_NMM_INTEGRATE
      USE MODULE_GFS_INTEGRATE
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State),   INTENT(IN)    :: IMP_STATE                       !<-- The ATM Run step's import
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM Run step's export
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_MAIN                      !<-- The main ESMF Clock
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_RUN                          !<-- Return code for the Run step
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)         :: CF  
      TYPE(WRAP_ATM_INTERNAL_STATE)    :: WRAP                             !<-- The F90 wrap of the ATM internal state
      TYPE(ATM_INTERNAL_STATE),POINTER :: ATM_INT_STATE                    !<-- The ATM internal state pointer
      TYPE(ESMF_Time)                  :: CURRTIME                         !<-- The current forecast time
      TYPE(ESMF_Time)                  :: STARTTIME                        !<-- The current forecast time
      TYPE(ESMF_Time)                  :: SDFITIME
      TYPE(ESMF_Time)                  :: DFITIME                          !<-- Digital filter time interval
      TYPE(ESMF_Time)                  :: HALFDFITIME                      !<-- Digital filter time interval
      TYPE(ESMF_TimeInterval)          :: HALFDFIINTVAL                    !<-- Digital filter time interval
      TYPE(ESMF_TimeInterval)          :: TIMESTEP                         !<-- Digital filter time interval
!
       TYPE(ESMF_Time)                  :: REFERENCE
       TYPE(ESMF_TimeInterval)          :: RUNDURATION
       TYPE(ESMF_TimeInterval)          :: DIFF
       INTEGER(KIND=KINT)    :: DFIHR                       !<-- Digital filter time interval

      INTEGER(KIND=KINT)         :: RC,NUM_PES_FCST                        !<-- Error signal variable.
      INTEGER(KIND=KINT)         :: I,IER,J
      INTEGER(KIND=KINT)         :: YY,MM,DD,H,M,S,NDFISTEP
      INTEGER(KIND=KINT),SAVE    :: NTIMESTEP                              !<-- The current forecast timestep (INT)
      INTEGER(KIND=ESMF_KIND_I8) :: NTIMESTEP_ESMF                         !<-- The current forecast timestep (ESMF_INT)
      CHARACTER(50)              :: MODE
      INTEGER(KIND=KINT)         :: NUM_TRACERS_MET,NUM_TRACERS_CHEM,MEAN_ON
      INTEGER(KIND=KINT)         :: FILTER_METHOD,HALF_FILTER_DURATION, HALF_BCK_DURATION
      LOGICAL                    :: PHYSICS_ON,FILTER_ON
      
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      btim0=timef()
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE ATM COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component                     
                            ,config  =CF                                 &  !<-- The config object (~namelist)
                            ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC_RUN)
     CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value  =NUM_TRACERS_MET               &  !<-- The variable filled (number of meteorological tracers)
                                  ,label ='num_tracers_met:'           &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC_RUN)
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =NUM_TRACERS_CHEM              &  !<-- The variable filled (number of chemical tracers)
                                  ,label ='num_tracers_met:'            &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC_RUN)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      IF (CORE=='nmm') THEN
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve ATM Component's Internal State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGetInternalState(ATM_GRID_COMP                  &  !<-- The ATM gridded component
                                        ,WRAP                           &  !<-- The F90 wrap of the ATM internal state
                                        ,RC)
!
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ATM_INT_STATE=>wrap%ATM_INT_STATE
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE TIMESTEP FROM THE ATM CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Timestep from the ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_ATM                         &
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- # of times the clock has advanced
                        ,rc          =RC)
!
      NTIMESTEP=NTIMESTEP_ESMF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ELSE IF (CORE=='gfs') THEN
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Timestep from the ATM MAIN"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock       =CLOCK_MAIN                         &
                        ,advanceCount=NTIMESTEP_ESMF                    &  !<-- # of times the clock has advanced
                        ,rc          =RC)
!
      NTIMESTEP=NTIMESTEP_ESMF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ENDIF 
      
!-----------------------------------------------------------------------
!***  RETRIEVE THE VM FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      IF (CORE=='nmm') THEN
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve VM from ATM Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- the ATM gridded component
                           ,vm      =VM                                 &
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_RUN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ELSE IF (CORE=='gfs') THEN
      MESSAGE_CHECK="Retrieve VM from ATM Component"
      call esmf_gridcompget(gridcomp=ATM_GRID_COMP                      &
                           ,vm      =VM                                 &
                           ,config  =CF                                 &
                           ,rc      =RC)
!
      call esmf_vmget(     vm       =VM                                 &
                           ,localpet=MYPE                               &
                     	   ,rc      =RC)
      ENDIF

!---------------------------------------------------------
!***  condition to run only adiabatic (dynamics only)
!---------------------------------------------------------
!
      CALL ESMF_ConfigGetAttribute(       CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =RC)
!
      IF (CORE=='nmm') THEN
      if(trim(mode)=='true')then
        PHYSICS_ON=.false.
        print *,' run without physics coupling '
      else
        PHYSICS_ON=.true.
        print *,' run with physics coupling.'
      endif
      ELSE IF (CORE=='gfs') THEN
      if(trim(mode)=='.true.')then
        PHYSICS_ON=.false.
        print *,' run without physics coupling.'
      else
        PHYSICS_ON=.true.
        print *,' run with physics coupling.'
      endif
      ENDIF

!
!
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
! ------------------------------------------------------------------
!
       IF (CORE=='nmm') THEN

       CALL ESMF_ClockGet(clock =CLOCK_ATM                          &
                          ,starttime =STARTTIME                         &
                          ,currtime =CURRTIME                           &
                          ,timestep =TIMESTEP                           &
                          ,runduration=RUNDURATION                      &
                          ,rc =RC)


       CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =FILTER_METHOD                 &
                                  ,label ='filter_method:'               &
                                  ,rc    =RC)

       IF (FILTER_METHOD .eq. 0) then

       print *,'no filtering.'


!-------------------------------------------------------------------

       ELSE IF (FILTER_METHOD .eq. 1) THEN


       CALL ESMF_ConfigGetAttribute(config=CF                     &  !<-- The config object
                                   ,value = DFIHR                 &
                                   ,label ='nsecs_dfl:'         &  !<-- Give this label's value
                                   ,rc    =RC)


       CALL ESMF_TimeIntervalSet(HALFDFIINTVAL                      &
                                 ,s=DFIHR,rc=RC)

       CALL NMM_FWD_INTEGRATE(ATM_GRID_COMP                        &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
                                  ,STARTTIME                            &
                                  ,HALFDFIINTVAL                        &
                                  ,FILTER_METHOD                        &
                                  ,TIMESTEP                             &
                                  ,MYPE                                 &
                                  ,NUM_TRACERS_MET                      &
                                  ,NUM_TRACERS_CHEM                     &
                                  ,NTIMESTEP                            &
                                  ,NPE_PRINT)
 

!----------------------------------------------------------------------

        ELSE IF (FILTER_METHOD .EQ. 2 ) THEN

         CALL ESMF_ConfigGetAttribute(config=CF                     &  !<-- The config object
                                  ,value = DFIHR                  &
                                  ,label ='nsecs_bckddfi:'         &  !<-- Give this label's value
                                  ,rc    =RC)

         CALL ESMF_TimeIntervalSet(HALFDFIINTVAL                      &
                                 ,s=DFIHR,rc=RC)


         CALL NMM_BCK_INTEGRATE(ATM_GRID_COMP                        &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
                                  ,STARTTIME                            &
				  ,HALFDFIINTVAL                        &
                                  ,FILTER_METHOD                        &
                                  ,TIMESTEP                             &
                                  ,MYPE                                 &
                                  ,NUM_TRACERS_MET                      &
                                  ,NUM_TRACERS_CHEM                     &
                                  ,NTIMESTEP                            &
                                  ,NPE_PRINT)



         CALL ESMF_ConfigGetAttribute(config=CF            &  !<-- The config object
                                  ,value = DFIHR           &
                                  ,label ='nsecs_fwdddfi:'      &  !<-- Give this label's value
                                  ,rc    =RC)

         CALL ESMF_TimeIntervalSet(HALFDFIINTVAL                      &
                                 ,s=DFIHR,rc=RC)


         CALL NMM_FWD_INTEGRATE(ATM_GRID_COMP                        &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
                                  ,STARTTIME                          &
                                  ,HALFDFIINTVAL                        &
				  ,FILTER_METHOD                        &
                                  ,TIMESTEP                             &
                                  ,MYPE                                 &
                                  ,NUM_TRACERS_MET                      &
                                  ,NUM_TRACERS_CHEM                     &
                                  ,NTIMESTEP                            &
                                  ,NPE_PRINT)

!----------------------------------------------------------------------

        ELSE IF (FILTER_METHOD .EQ. 3 ) THEN

         CALL ESMF_ConfigGetAttribute(config=CF                     &  !<-- The config object
                                  ,value = DFIHR                  &
                                  ,label ='nsecs_bcktdfi:'         &  !<-- Give this label's value
                                  ,rc    =RC)

         CALL ESMF_TimeIntervalSet(HALFDFIINTVAL                      &
                                 ,s=DFIHR,rc=RC)


         CALL NMM_BCK_INTEGRATE(ATM_GRID_COMP                        &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
                                  ,STARTTIME                            &
                                  ,HALFDFIINTVAL                        &
                                  ,FILTER_METHOD                        &
                                  ,TIMESTEP                             &
                                  ,MYPE                                 &
                                  ,NUM_TRACERS_MET                      &
                                  ,NUM_TRACERS_CHEM                     &
                                  ,NTIMESTEP                            &
                                  ,NPE_PRINT)


         CALL ESMF_ConfigGetAttribute(config=CF            &  !<-- The config object
                                  ,value = DFIHR           &
                                  ,label ='nsecs_fwdtdfi:'      &  !<-- Give this label's value
                                  ,rc    =RC)

         CALL ESMF_TimeIntervalSet(HALFDFIINTVAL                      &
                                 ,s=DFIHR,rc=RC)


         CALL NMM_FWD_INTEGRATE(ATM_GRID_COMP                        &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
                                  ,STARTTIME                            &
                                  ,HALFDFIINTVAL                        &
                                  ,FILTER_METHOD                        &
                                  ,TIMESTEP                             &
                                  ,MYPE                                 &
                                  ,NUM_TRACERS_MET                      &
                                  ,NUM_TRACERS_CHEM                     &
                                  ,NTIMESTEP                            &
                                  ,NPE_PRINT)

        ENDIF

            CALL NMM_INTEGRATE(ATM_GRID_COMP                        &
                                  ,ATM_INT_STATE                        &
                                  ,CLOCK_ATM                            &
                                  ,CURRTIME                             &
				  ,STARTTIME                            &
                                  ,TIMEINTERVAL_CLOCKTIME               &
                                  ,TIMEINTERVAL_HISTORY                 &
                                  ,TIMEINTERVAL_RESTART                 &
                                  ,MYPE                                 &
                                  ,NUM_TRACERS_MET                      &
                                  ,NUM_TRACERS_CHEM                     &
                                  ,NTIMESTEP                            &
                                  ,NPE_PRINT                            &
                                  ,PHYSICS_ON                           &
                                  ,RESTARTED_RUN)

       ELSE IF (CORE=='gfs') THEN

           CALL ESMF_ClockGet(clock =CLOCK_MAIN                          &
                          ,starttime =STARTTIME                         &
                          ,currtime =CURRTIME                           &
                          ,timestep =TIMESTEP                           &
                          ,runduration=RUNDURATION                      &
                          ,rc =RC)

           CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =DFIHR                         &
                                  ,label ='nhours_dfini:'               &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!jw
           CALL GFS_INTEGRATE(gc_gfs_dyn                            &
                                 ,gc_gfs_phy                            &
                                 ,gc_atm_cpl                            &
                                 ,wrt_comps                             &
                                 ,imp_gfs_dyn                           &
                                 ,exp_gfs_dyn                           &
                                 ,imp_gfs_phy                           &
                                 ,exp_gfs_phy                           &
                                 ,imp_gfs_wrt                           &
                                 ,exp_gfs_wrt                           &
                                 ,CLOCK_MAIN                            &
                                 ,TIMEINTERVAL_gfs_output               &
                                 ,quilting                              &
                                 ,WRITE_GROUP_READY_TO_GO               &
                                 ,CURRTIME                              &
                                 ,STARTTIME                             &
                                 ,NTIMESTEP                             &
                                 ,TIMESTEP                              &
                                 ,DFIHR                                 &
                                 ,MYPE                                  &
                                 ,PHYSICS_ON)
      ENDIF 



      CALL ESMF_ClockPrint(   clock=CLOCK_MAIN                      &
                              ,options="stoptime string"            &
                              ,rc=rc)


!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_RUN
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      SUBROUTINE ATM_FINALIZE(ATM_GRID_COMP                             &
                             ,IMP_STATE,EXP_STATE                       &
                             ,CLOCK_MAIN                                &
                             ,RC_FINALIZE)
!
!-----------------------------------------------------------------------
!***  THIS ROUTINE FINALIZES THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT) :: ATM_GRID_COMP                   !<-- The ATM gridded component
      TYPE(ESMF_State),   INTENT(INOUT) :: IMP_STATE                       !<-- The ATM finalize step's import state
      TYPE(ESMF_State),   INTENT(INOUT) :: EXP_STATE                       !<-- The ATM finalize step's export state
      TYPE(ESMF_Clock),   INTENT(INOUT) :: CLOCK_MAIN                      !<-- The main ESMF Clock
      TYPE(ESMF_Config)                 :: CF                               !<-- The config object
      INTEGER,OPTIONAL,   INTENT(OUT)   :: RC_FINALIZE                     !<-- Return code for the Finalize step
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(WRAP_ATM_INTERNAL_STATE)    :: WRAP                             !<-- The F90 wrap of the ATM internal state
      TYPE(ATM_INTERNAL_STATE),POINTER :: ATM_INT_STATE                    !<-- The ATM internal state pointer
!
      INTEGER :: I,J
      INTEGER :: RC,RC_FINAL                                                ! The final error signal variables.
      CHARACTER(50):: MODE
      LOGICAL    :: PHYSICS_ON     

!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC         =ESMF_SUCCESS
      RC_FINAL   =ESMF_SUCCESS
      RC_FINALIZE=ESMF_SUCCESS
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Retrieve Config Object from ATM Component"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Dynamic Core Name"
      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Give this label's value to the previous variable
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!---------------------------------------------------------
!***  condition to run only adiabatic (dynamics only)
!---------------------------------------------------------
!
      CALL ESMF_ConfigGetAttribute(       CF                            &
                                  ,value =MODE                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =RC)
!
      IF (CORE=='nmm') THEN
      if(trim(mode)=='true')then
        PHYSICS_ON=.false.
        print *,' finalize without physics coupling. '
      else
        PHYSICS_ON=.true.
        print *,' finalize with physics coupling. '
      endif
      ELSE IF (CORE=='gfs') THEN
      if(trim(mode)=='.true.')then
        PHYSICS_ON=.false.
        print *,' finalize without physics coupling. '
      else
        PHYSICS_ON=.true.
        print *,' finalize with physics coupling. '
      endif
      ENDIF


!-----------------------------------------------------------------------

!
!-----------------------------------------------------------------------
!***  RETRIEVE THE ATM GRIDDED COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompGetInternalState(ATM_GRID_COMP                  &  !<-- The ATM gridded component
                                        ,WRAP                           &  !<-- The F90 wrap of the ATM internal state
                                        ,RC)
!
      ATM_INT_STATE=>wrap%ATM_INT_STATE
      
      ENDIF 
!
!-----------------------------------------------------------------------
!***  FINALIZE EACH OF THE SUBCOMPONENTS.
!-----------------------------------------------------------------------
!
!--------------
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%DYN_GRID_COMP &
                                ,importState=atm_int_state%IMP_STATE_DYN &
                                ,exportState=atm_int_state%EXP_STATE_DYN &
                                ,clock      =CLOCK_ATM                   &
                                ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompfinalize(gridcomp   =gc_gfs_dyn                  &
                               ,importstate=imp_gfs_dyn                 &
                               ,exportstate=exp_gfs_dyn                 &
                               ,clock      =CLOCK_MAIN                  &
                               ,rc         =RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! - FIX later
 RC=ESMF_SUCCESS
 RC_FINAL=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (PHYSICS_ON) THEN
!-------------
!***  PHYSICS
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompFinalize(gridcomp   =atm_int_state%PHY_GRID_COMP &
                                ,importState=atm_int_state%IMP_STATE_PHY &
                                ,exportState=atm_int_state%EXP_STATE_PHY &
                                ,clock      =CLOCK_ATM                   &
                                ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_gridcompfinalize(gridcomp   =gc_gfs_phy                 &
                                ,importstate=imp_gfs_phy                &
                                ,exportstate=exp_gfs_phy                &
                                ,clock      =CLOCK_MAIN                 &
                                ,rc         =RC)
      ENDIF
      ENDIF !PHYSICS_ON
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!ratko    CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! - FIX later
 RC=ESMF_SUCCESS
 RC_FINAL=ESMF_SUCCESS
!ratko
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------
!***  DYNAMICS-PHYSICS COUPLER
!-----------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Finalize Dynamics-Physics Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_CplCompFinalize(cplcomp    =atm_int_state%COUPLER_DYN_PHY_COMP &
                               ,importState=atm_int_state%EXP_STATE_DYN        &
                               ,exportState=atm_int_state%IMP_STATE_PHY        &
                               ,clock      =CLOCK_ATM                          &
                               ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      call esmf_cplcompfinalize(cplcomp    = gc_atm_cpl                  &
                               ,importstate=exp_gfs_dyn                 &
                               ,exportstate=imp_gfs_phy                 &
                               ,clock      =CLOCK_MAIN                  &
                               ,rc         =RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  DESTROY ALL STATES.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy States"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_StateDestroy(state=atm_int_state%IMP_STATE_DYN          &
                            ,rc   =RC)
!
      CALL ESMF_StateDestroy(state=atm_int_state%EXP_STATE_DYN          &
                            ,rc   =RC)
!
      IF (PHYSICS_ON) THEN
      CALL ESMF_StateDestroy(state=atm_int_state%IMP_STATE_PHY          &
                            ,rc   =RC)
!
      CALL ESMF_StateDestroy(state=atm_int_state%EXP_STATE_PHY          &
                            ,rc   =RC)
      ENDIF
      ELSE IF (CORE=='gfs') THEN
      call esmf_statedestroy(state=imp_gfs_dyn, rc=RC)
      call esmf_statedestroy(state=exp_gfs_dyn, rc=RC)
      IF (PHYSICS_ON) THEN
      call esmf_statedestroy(state=imp_gfs_phy, rc=RC)
      call esmf_statedestroy(state=exp_gfs_phy, rc=RC)
      ENDIF
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  IF QUILTING WAS SELECTED FOR THE GENERATION OF OUTPUT,
!***  FINALIZE AND DESTROY OBJECTS RELATED TO THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      IF (CORE=='nmm') THEN
      IF(atm_int_state%QUILTING)THEN
        CALL WRITE_DESTROY(ATM_GRID_COMP,ATM_INT_STATE,CLOCK_ATM)
      ENDIF
!
!-----------------------------------------------------------------------
!***  DESTROY THE ATM CLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockDestroy(clock=CLOCK_ATM                            &
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       ENDIF
!
!jw for GFS destroy WRT grid component
      IF (CORE=='gfs') THEN
        IF(QUILTING) THEN
          CALL WRITE_DESTROY_GFS(ATM_GRID_COMP,WRT_COMPS,               &
            IMP_GFS_WRT,EXP_GFS_WRT,CLOCK_MAIN)
        ENDIF
      ENDIF
!
!-----------------------------------------------------------------------
!***  DESTROY ALL SUBCOMPONENTS.
!-----------------------------------------------------------------------
!
!--------------
!***  DYNAMICS
!--------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Dynamics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompDestroy(gridcomp=atm_int_state%DYN_GRID_COMP      &
                               ,rc      =RC)
      ELSE 
      call esmf_gridcompdestroy(gridcomp=gc_gfs_dyn, rc=RC)
      ENDIF
    
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (PHYSICS_ON) THEN
!-------------
!***  PHYSICS
!-------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Physics Component"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_GridCompDestroy(gridcomp=atm_int_state%PHY_GRID_COMP      &
                               ,rc      =RC)
      ELSE
      call esmf_gridcompdestroy(gridcomp=gc_gfs_phy, rc=RC)
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!------------------------------
!***  DYNAMICS-PHYSICS COUPLER
!------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy Dynamics-Physics Coupler"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN
      CALL ESMF_CplCompDestroy(cplcomp=atm_int_state%COUPLER_DYN_PHY_COMP &
                              ,rc     =RC)
      ELSE
      call esmf_cplcompdestroy(cplcomp=gc_atm_cpl, rc=RC)
      ENDIF
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_FINAL)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  THE FINAL ERROR SIGNAL INFORMATION.
!-----------------------------------------------------------------------
!
      IF(RC_FINAL==ESMF_SUCCESS)THEN
        WRITE(0,*)'ATM FINALIZE STEP SUCCEEDED'
      ELSE
        WRITE(0,*)'ATM FINALIZE STEP FAILED'
      ENDIF
!
      IF(PRESENT(RC_FINALIZE))THEN
        RC_FINALIZE=RC_FINAL
      ENDIF
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE ATM_FINALIZE
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
