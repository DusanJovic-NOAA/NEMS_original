!-----------------------------------------------------------------------
!
      PROGRAM MAIN_NEMS
!
!-----------------------------------------------------------------------
!***  Main Program for NEMS.
!***  Define ESMF data types and procedures.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
! PROGRAM HISTORY LOG:
!   2007-       Black - Modified from Wei-yu's version
!   2007-09     Black - Create the Clock here.
!   2009-08     Colon - Unified NEM-NMM & NEMS-GFS
!
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
!
!-----------------------------------------------------------------------
!***  USE the model's ATM gridded component module.  Although it
!***  contains the calls to Register and the top level Initialize,
!***  Run, and Finalize, only the Register routine is public.
!-----------------------------------------------------------------------
!
!      USE module_ATM_GRID_COMP,ONLY: ATM_REGISTER
       USE module_ATM_GRID_COMP
!
!-----------------------------------------------------------------------
!***  The following module contains error-checking.
!-----------------------------------------------------------------------
!
!      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
       USE MODULE_ERR_MSG
!
!-----------------------------------------------------------------------
!***  This is some information required by the ESMF Log error utility.
!-----------------------------------------------------------------------
!
#include "ESMF_LogMacros.inc"
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  Local variables.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      TYPE(ESMF_VM)       :: VM                                          !<-- The ESMF virtual machine,
                                                                         !    which contains and manages
                                                                         !    the computer CPU resource
                                                                         !    for the ESMF grid components.
      TYPE(ESMF_GridComp) :: ATM_GRID_COMP                               !<-- The ATM gridded component.
!
      TYPE(ESMF_State)    :: ATM_STATE                                   !<-- This state will be used
                                                                         !    as both import and export
                                                                         !    for the ATM component.
      
      TYPE(ESMF_GridComp), allocatable :: gc_atm(:)     ! grid component.
      TYPE(ESMF_State)                 :: imp_atm       ! import state
      TYPE(ESMF_State)                 :: exp_atm       ! import state




!
      TYPE(ESMF_Clock)    :: CLOCK_MAIN                                  !<-- The ESMF time management clock
!
      TYPE(ESMF_Clock)	  :: clock
      TYPE(ESMF_Config)   :: CF                                          !<-- The Configure object
!
      TYPE(ESMF_LOG)      :: LOG_ERROR                                   !<-- The ESMF Log Error object.
!
!-----------------------------------------------------------------------
!
      INTEGER                 :: YY,MM,DD                                !<-- Time variables for date
      INTEGER                 :: HH,MNS,SEC                              !<-- Time variables for time of day
      INTEGER                 :: NHOURS_FCST                             !<-- Length of forecast in hours
      INTEGER                 :: NSECONDS_FCST                           !<-- Length of forecast in seconds
!
      TYPE(ESMF_TimeInterval) :: RUNDURATION                             !<-- The ESMF time. The total forecast hours.
!
      TYPE(ESMF_TimeInterval) :: TIMESTEP                                !<-- The ESMF timestep length (we only need a dummy here)
!
      TYPE(ESMF_Time)         :: STARTTIME                               !<-- The ESMF start time.
!
      INTEGER                 :: RC=ESMF_SUCCESS                         !<-- The running error signal
!
      INTEGER                 :: RC_MAIN                                 !<-- The final value of the
      INTEGER,ALLOCATABLE     :: pe_member(:), petlist(:,:)
      INTEGER                 :: member_id, i, j, ij, tasks, pe_max, me, total_member, mype, num_pes
      INTEGER		      :: timestep_sec_whole, timestep_sec_numerator, timestep_sec_denominator
!
!-----------------------------------------------------------------------
!***  Declare the ATM gridded component and state names.
!-----------------------------------------------------------------------
!
      CHARACTER(ESMF_MAXSTR) :: ATM_GRID_COMP_NAME='ATM Gridded Component'  &
                               ,ATM_STATE_NAME    ='ATM State'              &
                               ,CF_NAME           ='Configure_File'
      CHARACTER(ESMF_MAXSTR),ALLOCATABLE :: gridcompname(:)
      CHARACTER(ESMF_MAXSTR) :: impstatename, expstatename
      CHARACTER*20           :: pelab
      CHARACTER(3)           :: CORE
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Initialize the final error signal.
!-----------------------------------------------------------------------
!
      RC_MAIN=ESMF_SUCCESS
!
!-----------------------------------------------------------------------
!***  Initialize the ESMF framework. 
!-----------------------------------------------------------------------
! 
      CALL ESMF_Initialize(VM             =VM                           & !<-- The ESMF Virtual Machine
                          ,defaultCalendar=ESMF_CAL_GREGORIAN           & !<-- Set up the default calendar.
                          ,defaultlogtype =ESMF_LOG_MULTI               & !<-- Define multiple log error output file;
                                                                          !    each task has its own log error output file.
                          ,rc             =RC)
!
!-----------------------------------------------------------------------
!***  Open the log error file and set log parameters.
!***  When the user does not want to use the ESMF default log,
!***  use it to open a user defined log.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Open the Log Error File"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_LogOpen(log=LOG_ERROR                                   &
                       ,filename='Error_Log.txt'                        &
                       ,logtype=ESMF_LOG_MULTI                          &
                       ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Set up the user defined log.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Set Up ESMF Log"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_LogSet(log        =LOG_ERROR                            &
                      ,verbose    =ESMF_TRUE                            &
                      ,flush      =ESMF_TRUE                            &
                      ,rootOnly   =ESMF_FALSE                           &
                      ,halt       =ESMF_LOG_HALTERROR                   & ! It means that the job will be stopped
                                                                          ! when error happens.
                      ,maxElements=1                                    & ! Maximum number of elements in the log
                                                                          ! before printing them to the log file.
                      ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Create and load the Configure object which will hold the contents
!***  of the configure file.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create/Load the Configure Object"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CF=ESMF_ConfigCreate(rc=RC)
!
      CALL ESMF_ConfigLoadFile(config  =CF                              & !<-- The Configure object
                              ,filename='configure_file'   & !<-- The name of the configure file 
                              ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the ATM ESMF gridded component.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      
      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The config object
                                  ,value =CORE                          &  !<-- The variable filled (dynamic core name)
                                  ,label ='core:'                       &  !<-- Give this label's value to
                                  ,rc    =RC)

      MESSAGE_CHECK="Create the ATM Gridded Component"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='gfs') THEN
      OPEN (unit=50, FILE = 'atm_namelist.rc',POSITION = 'append',STATUS = 'OLD')
      WRITE(50,*)'core: ',CORE
      CLOSE(50)
      CALL ESMF_ConfigDestroy(CF,RC)
      CF=ESMF_ConfigCreate(rc=RC)
      CALL ESMF_ConfigLoadFile(config  =CF                              & !<-- The Configure object
                              ,filename='atm_namelist.rc'   & !<-- The name of the configure file
                              ,rc      =RC)
      ENDIF

      IF(CORE=='nmm') THEN
      ATM_GRID_COMP=ESMF_GridCompCreate(name        =ATM_GRID_COMP_NAME & !<-- ATM gridded component name
                                       ,gridcomptype=ESMF_ATM           & !<-- A Type that the user can query but it has
                                                                          !    no meaning in the job run.  Can also be
                                                                          !    ESMF_OCEAN, ESMF_LAND, etc.
                                       ,configFile  ='configure_file'   & !<-- Link the user-created configure file to the 
                                                                          !    ATM gridded component.
                                       ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ELSE IF(CORE=='gfs') THEN
!--------------------------------------------------
! to get the total number of the ensemble members.
!--------------------------------------------------
      call ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =total_member                  &
                                  ,label = 'total_member:'              &
                                  ,rc    = RC)

!-----------------
! allocate arrays.
!-----------------
      ALLOCATE(gc_atm(total_member))
      ALLOCATE(gridcompname(total_member))
      ALLOCATE(pe_member(total_member))
      
      call ESMF_VmGet(VM                                                &
                     , pecount = tasks    			        &	 
                     , localpet = me                                    &
                     , rc = rc)

      do i = 1, total_member
        write(gridcompname(i), '("atm main grid component",i2.2)') i
        write(impstatename   , '("atm import state")')
        write(expstatename   , '("atm export state")')


        call ESMF_ConfigGetAttribute(config=CF                          &
			          ,value = pe_member(i)                 & 
                                  ,label = pelab                        & 
				  ,rc = RC)
        if (pe_member(i) == 0) pe_member(i) = tasks / total_member
      enddo
      pe_max = 1
      do i=1,total_member
        pe_max = max(pe_max,pe_member(i))
      enddo

!----------------------
!  set up the pet list.
!----------------------
      ALLOCATE(petlist(pe_max, total_member))
      ij = 0
      do j = 1, total_member
          do i = 1, pe_member(j)
              petlist(i, j) = ij
              if(me == ij) then
                  member_id = j
              end if
              ij = ij+1
          end do
      end do

!-----------------------------------------------------------------------
      do i = 1, total_member
      gc_atm(i) = ESMF_GridCompCreate (                                  &
                              name         = gridcompname(i)             &
                             ,gridcomptype = esmf_atm                    &
                             ,petlist      = petlist(1:pe_member(i),i)  &
                             ,config       = CF                          &
                             ,rc           = RC)
      end do
     CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
     ENDIF 
!
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Register ATM Init, Run, Finalize"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF (CORE=='nmm') THEN	
      CALL ESMF_GridCompSetServices(ATM_GRID_COMP                        &  !<-- The ATM gridded component
                                   ,ATM_REGISTER                         &  !<-- User's subroutineName 
                                   ,RC)
      ELSE IF (CORE=='gfs') THEN
      RC   = ESMF_SUCCESS
      do i = 1, total_member
      CALL ESMF_GridCompSetServices(gc_atm(i)                            &  !<-- The ATM gridded component
                                   ,ATM_REGISTER                         &  !<-- User's subroutineName
                                   ,RC)
      end do
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Create the ATM component's ESMF state.
!***  It will be used for both import and export because 
!***  it will not be required to perform an actual function.
!***  When ATM is coupled to another component then the 
!***  import/export states will actually be utilized.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the ATM Import/Export State"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

      IF (CORE=='nmm') THEN
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ATM_STATE=ESMF_StateCreate(statename=ATM_STATE_NAME               &
                                 ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      ELSE IF (CORE=='gfs') THEN

      imp_atm = ESMF_StateCreate (statename = impstatename  &
                                  ,statetype = esmf_state_import &
                                  ,rc        = RC)

      exp_atm = ESMF_StateCreate (statename = expstatename   &
                                  ,statetype = esmf_state_export &
                                  ,rc        = RC)
      ENDIF
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
!
!-----------------------------------------------------------------------
!***  Create the ESMF Clock which will control the timestepping.
!***  The Clock is needed for all calls to Init, Run, and Finalize.
!
!***  Obtain the simulation start and end times from the Configure File.
!***  Each of the primary Components will create and use their own
!***  clocks since the timestepping within those Components will in
!***  general be different from each other but they will read the
!***  start/end times from the main Clock.
!-----------------------------------------------------------------------
!
      IF (CORE=='nmm') THEN
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create dummy Timestep in MAIN"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=TIMESTEP                   &  !<-- Dummy fundamental timestep (sec) (ESMF)
                               ,s           =1                          &  !<-- Dummy whole seconds
                               ,sn          =0                          &  !<-- Dummy numerator of timestep fraction
                               ,sd          =1                          &  !<-- Dummy denominator of timestep fraction
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ELSE IF (CORE=='gfs')THEN
!-----------------------------------------------------------------------
!***  EXTRACT FUNDAMENTAL TIMESTEP INFORMATION FROM THE CONFIG FILE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Extract Timestep Information from Config File"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
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
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
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
     ENDIF
!---------------
!***  Start Time
!---------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Year from Config File"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =YY                            &
                                  ,label ='start_year:'                 &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Month from Config File"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =MM                            &
                                  ,LABEL ='start_month:'                &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Day from Config File"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =DD                            &
                                  ,label ='start_day:'                  &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Hour from Config File"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =HH                            &
                                  ,label ='start_hour:'                 &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Minute from Config File"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =MNS                           &
                                  ,label ='start_minute:'               &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Starting Second from Config File"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =SEC                           &
                                  ,label ='start_second:'               &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Set the Forecast Start Time"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeSet(time=STARTTIME                                  &  !<-- The start time of the forecast (ESMF)
                       ,yy  =YY                                         &  !<-- Year from config file
                       ,mm  =MM                                         &  !<-- Month from config file
                       ,dd  =DD                                         &  !<-- Day from config file
                       ,h   =HH                                         &  !<-- Hour from config file
                       ,m   =MNS                                        &  !<-- Minute from config file
                       ,s   =SEC                                        &  !<-- Second from config file
                       ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------
!***  Run Duration
!-----------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Extract Forecast Length from Config File"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ConfigGetAttribute(config=CF                            &
                                  ,value =NHOURS_FCST                   &
                                  ,label ='nhours_fcst:'                &
                                  ,rc    =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      NSECONDS_FCST=NHOURS_FCST*3600                                       !<-- The forecast length (s) (REAL)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="MAIN: Set the Forecast Length"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeIntervalSet(timeinterval=RUNDURATION                &  !<-- The forecast length (s) (ESMF)
                               ,s           =NSECONDS_FCST              &
                               ,rc          =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Now the Main Clock can be created.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Create the Main Clock"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

      CLOCK_MAIN=ESMF_ClockCreate(name       ='MAIN_CLOCK'              &  !<-- The top-level ESMF Clock
                                 ,timestep   =TIMESTEP                  &  !<-- A dummy timestep needed by the Clock
                                 ,starttime  =STARTTIME                 &  !<-- The forecast start time
                                 ,runduration=RUNDURATION               &  !<-- The forecast duration
                                 ,rc         =RC)

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Execute the INITIALIZE step for the ATM gridded component.
!***  The initialize routine that is called here as well as the 
!***  run and finalize routines invoked below are those specified
!***  in the Register routine called in ESMF_GridCompSetServices above.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the ATM Component Initialize Step"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(CORE=='nmm') THEN 
      CALL ESMF_GridCompInitialize(gridcomp   =ATM_GRID_COMP            &  !<-- the ATM gridded component
                                  ,importstate=ATM_STATE                &  !<-- the ATM component's import state
                                  ,exportstate=ATM_STATE                &  !<-- the ATM component's export state
                                  ,clock      =CLOCK_MAIN               &  !<-- the ESMF clock
                                  ,phase      =ESMF_SINGLEPHASE         &  
                                  ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
      do i = 1, total_member
      if(member_id == i) then
      CALL ESMF_GridCompInitialize(gridcomp   =gc_atm(i)                &  !<-- the ATM gridded component
                                  ,importstate=imp_atm                  &  !<-- the ATM component's import state
                                  ,exportstate=exp_atm                  &  !<-- the ATM component's export state
                                  ,clock      =CLOCK_MAIN               &  !<-- the ESMF clock
                                  ,phase      =ESMF_SINGLEPHASE         &
                                  ,rc         =RC)
      end if
      end do
      ENDIF
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Execute the RUN step for the ATM gridded component.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the ATM Component Run Step"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(CORE=='nmm') THEN
      CALL ESMF_GridCompRun(gridcomp   =ATM_GRID_COMP                   &  !<-- the ATM gridded component
                           ,importstate=ATM_STATE                       &  !<-- the ATM component's import state
                           ,exportstate=ATM_STATE                       &  !<-- the ATM component's export state
                           ,clock      =CLOCK_MAIN                      &  !<-- the ESMF clock
                           ,phase      =ESMF_SINGLEPHASE                &
                           ,rc         =RC)
       ELSE IF (CORE=='gfs') THEN
       do i = 1, total_member
       if(member_id == i) then
      CALL ESMF_GridCompRun(gridcomp   =gc_atm(i)                   &  !<-- the ATM gridded component
                           ,importstate=imp_atm                       &  !<-- the ATM component's import state
                           ,exportstate=exp_atm                       &  !<-- the ATM component's export state
                           ,clock      =CLOCK_MAIN                      &  !<-- the ESMF clock
                           ,phase      =ESMF_SINGLEPHASE                &
                           ,rc         =RC)
      end if
      end do
      ENDIF

!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  Execute the FINALIZE step for the ATM gridded component.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Execute the ATM Component Finalize Step"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      IF(CORE=='nmm') THEN
      CALL ESMF_GridCompFinalize(gridcomp   =ATM_GRID_COMP              &  !<-- the ATM gridded component
                                ,importstate=ATM_STATE                  &  !<-- the ATM component's import state
                                ,exportstate=ATM_STATE                  &  !<-- the ATM component's export state
                                ,clock      =CLOCK_MAIN                 &  !<-- the ESMF clock
                                ,phase      =ESMF_SINGLEPHASE           &
                                ,rc         =RC)
      ELSE IF (CORE=='gfs') THEN
       do i = 1, total_member
       if(member_id == i) then
      CALL ESMF_GridCompFinalize(gridcomp   =gc_atm(i)              &  !<-- the ATM gridded component
                                ,importstate=imp_atm                  &  !<-- the ATM component's import state
                                ,exportstate=exp_atm                  &  !<-- the ATM component's export state
                                ,clock      =CLOCK_MAIN                 &  !<-- the ESMF clock
                                ,phase      =ESMF_SINGLEPHASE           &
                                ,rc         =RC)
      end if
      end do
      ENDIF


!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!
!-----------------------------------------------------------------------
!***  Destroy the Main Clock.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Destroy the Main Clock"
!      CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockDestroy(clock=CLOCK_MAIN                           &
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_MAIN)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  Shut down the ESMF system.
!-----------------------------------------------------------------------
!
      CALL ESMF_Finalize()
!
      IF(RC_MAIN==ESMF_SUCCESS)THEN
        WRITE(0,*)'MODEL FORECAST RUN SUCCEEDED'
      ELSE
        WRITE(0,*)'MODEL FORECAST RUN FAILED  RC_MAIN=',RC_MAIN
      ENDIF
!
!-----------------------------------------------------------------------
!
      END PROGRAM MAIN_NEMS
!
!-----------------------------------------------------------------------
