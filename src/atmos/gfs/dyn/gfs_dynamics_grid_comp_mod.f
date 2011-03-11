#include "../../../ESMFVersionDefine.h"

! !module: gfs_dynamics_grid_comp_mod --- 
!                       esmf gridded component of gfs dynamics
!
! !description: gfs dynamics gridded component main module.
!
! !revision history:
!
!  january 2007     hann-ming henry juang initiated and wrote the code
!  March   2009     Weiyu Yang, modified for the ensemble NEMS run.
!  oct 4   2009     sarah lu, 3D Gaussian grid (DistGrid5) added
!  oct 5   2009     sarah lu, grid_gr unfolded from 2D to 3D
!  oct     2009     Jun Wang, add not quilting option, add run time variables into wrt inport state
!  oct 12 2009      Sarah Lu, set up the association between imp/exp state and grid_gr
!                   => export state is pointed to grid_gr in init step
!                   => grid_gr is pointed to import state in run step
!  oct 17 2009      Sarah Lu, add debug print to check imp/exp state
!  nov 09 2009      Jun Wang, add grid_gr_dfi for digital filter
!  Feb 05 2010      Jun Wang, add restart step
!  Feb 20 2011      Henry Juang, add non-iterating dimensional-splitting semi-Lagrangian (NDSL)
!                   advection with options of MASS_DP and NDSLFV
!  Feb    2011      Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!                   ESMF 5 library and the the ESMF 3.1.0rp2 library.
!
! !interface:
!
      module gfs_dynamics_grid_comp_mod
 
!!uses:
!------
      use esmf_mod

      use gfs_dynamics_err_msg_mod
      use gfs_dynamics_initialize_mod
      use gfs_dynamics_run_mod
      use gfs_dynamics_finalize_mod

      use gfs_dyn_mpi_def
      use gfs_dynamics_output, only : point_dynamics_output_gfs

      implicit none

      private   ! by default, data is private to this module

      public gfs_dyn_setservices	! only set service is public

!eop
!-------------------------------------------------------------------------


      contains


!----------------------------------------------------------------------
!bop
!
! !routine: gfs_dyn_setservices --- 
!           set services for gfs dynamics gridded component.
! 
! !interface:
!
      subroutine gfs_dyn_setservices (gc_gfs_dyn, rc)
 
! !arguments:
!------------

      type(esmf_gridcomp), intent(in)  :: gc_gfs_dyn 	! gridded component
      integer,             intent(out) :: rc    	! return code
     
! !description: set services (register) for the gfs dynamics grid component.
!         
!eop         
!----------------------------------------------------------------------
  
      integer                            :: rc1     = esmf_success

! initializing the error signal variable rc.
!-------------------------------------------
      rc = esmf_success

! register services for this component
! ------------------------------------

! register the initialize subroutine.  since it is just one subroutine
! for the initialize, use esmf_singlephase.  the second argument is
! a pre-defined subroutine type, such as esmf_setinit, esmf_setrun, 
! esmf_setfinal.
!---------------------------------------------------------------------
      call esmf_logwrite("set entry point for initialize",              &
                          esmf_log_info, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_dyn, 			&
                                       esmf_setinit,  			&
                                       gfs_dyn_initialize,              &
                                       esmf_singlephase, rc1)
      call gfs_dynamics_err_msg(rc1,'set entry point for initialize',rc)

! register the run subroutine.
!-----------------------------
      call esmf_logwrite("set entry point for run",              	&
                           esmf_log_info, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_dyn, 			&
                                       esmf_setrun,   			&
                                       gfs_dyn_run,	                &
                                       esmf_singlephase, rc1)
      call gfs_dynamics_err_msg(rc1,'set entry point for run',rc)


! register the finalize subroutine.
!----------------------------------
      call esmf_logwrite("set entry point for finalize",                &
                        esmf_log_info, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_dyn, 			&
                                       esmf_setfinal, 			&
                                       gfs_dyn_finalize, 	        &
                                       esmf_singlephase, rc1)
      call gfs_dynamics_err_msg(rc1,'set entry point for finalize',rc)

! check the error signal variable and print out the result.
!----------------------------------------------------------
      call gfs_dynamics_err_msg_final(rc1,				&
                        'setservice for gfs dynamics grid comp.',rc)

      end subroutine gfs_dyn_setservices





!----------------------------------------------------------------------
!bop
! !routine:  gfs_dyn_initialize --- initialize routine to initialize 
!                                   and set up the gfs running job.
!
! !description: this subroutine initializes the gfs running before
!               the main running loop.
!
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     h.-m. h. juang
!  oct 12 2009       Sarah Lu, export state is pointed to grid_gr once and for all
!  November 2009     Weiyu Yang, Ensemble GEFS.
!
! !interface:
!

! this argument list is a standard list for all the initialize,
! the run and finalize routines for an esmf system.
!--------------------------------------------------------------
      subroutine gfs_dyn_initialize(gc_gfs_dyn,		                &
                                   imp_gfs_dyn, exp_gfs_dyn, clock, rc)

! user code, for computations related to the esmf interface states.
!------------------------------------------------------------------
      use gfs_dynamics_states_mod, only : gfs_dynamics_import2internal, &
                                          gfs_dynamics_internal2export
      use gfs_dynamics_grid_create_mod
      USE GFS_AddParameterToStateMod
!
! !input/output variables and parameters:
!----------------------------------------

      type(esmf_gridcomp), intent(inout) :: gc_gfs_dyn 
      type(esmf_state),    intent(inout) :: imp_gfs_dyn
      type(esmf_state),    intent(inout) :: exp_gfs_dyn
      type(esmf_clock),    intent(inout) :: clock

!
! !output variables and parameters:
!----------------------------------

      integer, intent(out) :: rc  

! !eop
!------------------------------------------------------------------------- 
 
! !working arrays and local parameters.  
!--------------------------------------
      type(gfs_dyn_wrap)                :: wrap         
! this wrap is a derived type which contains
! only a pointer to the internal state.  it is needed
! for using different architectures or compliers.
      type(gfs_dynamics_internal_state), pointer  :: int_state    
      type(esmf_vm)                      :: vm_local     
      type(esmf_timeinterval)            :: timestep     
      type(esmf_timeinterval)            :: runduration  
      type(esmf_time)                    :: starttime    
      type(esmf_time)                    :: stoptime    
      type(esmf_time)                    :: currtime     
      type(esmf_timeinterval)            :: reftimeinterval 
!jw
      type(esmf_state)                   :: imp_state_write  !<-- The write gc import state

      integer(kind=esmf_kind_i4)         :: yy, mm, dd   ! time variables for date
      integer(kind=esmf_kind_i4)         :: hh, mns, sec ! time variables for time
      integer                            :: advancecount4, timestep_sec
      integer                            :: atm_timestep_s, dyn_timestep_s
      integer(esmf_kind_i8)              :: advancecount

      TYPE(ESMF_DistGrid)                :: DistGrid5    ! the ESMF DistGrid.

      integer                            :: rc1 
      integer                            :: rcfinal, grib_inp
      integer                            :: ifhmax
      integer                            :: runduration_hour 

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! allocate the internal state pointer.
!-------------------------------------
      call esmf_logwrite("allocate the dyn internal state",                 &
                         esmf_log_info, rc = rc1)

      allocate(int_state, stat = rc1)

      call gfs_dynamics_err_msg(rc1,' - allocate the internal state',rc)

      wrap%int_state => int_state

!jws
!-----------------------------------------------------------------------
!***  RETRIEVE THE IMPORT STATE OF THE WRITE GRIDDED COMPONENT
!***  FROM THE DYNAMICS EXPORT STATE.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("get write gc import state",                  &
                        esmf_log_info, rc = rc1)

      CALL ESMF_StateGet(state      =exp_gfs_dyn                        &  !<-- The Dynamics export state
                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Dynamics export state
                        ,nestedState=IMP_STATE_WRITE                    &  !<-- Extract write component import state from Dynamics export
                        ,rc         =RC)
      call gfs_dynamics_err_msg(rc1,'get write gc import state',rc)
!jwe

! attach internal state to the gfs dynamics grid component.
!-------------------------------------------------
      call esmf_logwrite("set up the internal state",                   &
                        esmf_log_info, rc = rc1)

      call esmf_gridcompsetinternalstate(gc_gfs_dyn, wrap, rc1)

      call gfs_dynamics_err_msg(rc1,'set up the internal state',rc)

! use esmf utilities to get information from the configuration file.
! the function is similar to reading the namelist in the original gfs.
!---------------------------------------------------------------------
      call esmf_logwrite("getting information from the configure file", &
                        esmf_log_info, rc = rc1)

      call gfs_dynamics_getcf(gc_gfs_dyn, int_state,  rc1)

      call gfs_dynamics_err_msg(rc1,'get configure file information',rc)

! get the start time from reading the sigma file.
!----------------------------------------------------------
      call esmf_logwrite("getting the start time",                      &
                         esmf_log_info, rc = rc1)

      call gfs_dynamics_start_time_get(					&
                        yy, mm, dd, hh, mns, sec, int_state%kfhour,     &
                        int_state%n1,int_state%n2,int_state%grib_inp,   &
                        int_state%nam_gfs_dyn%grid_ini,                 &
                        int_state%nam_gfs_dyn%grid_ini2, rc1)
 
      call gfs_dynamics_err_msg(rc1,'getting the start time',rc)
 
      advancecount4    = nint(real(int_state%kfhour) * 3600.0 /         &
                              int_state%nam_gfs_dyn%deltim)
      int_state%phour  = advancecount4 * 				&
                         int_state%nam_gfs_dyn%deltim / 3600.0
      int_state%kfhour = nint(int_state%phour)
!
      int_state%kdt    = advancecount4
!      print *,'in dyn_grid_comp,advancecount4=',advancecount4,          &
!        'phour=',int_state%phour,'kfhour=',int_state%kfhour,'kdt=',     &
!         int_state%kdt
!
! initialize the clock with the start time based on the information
! from calling starttimeget.
!------------------------------------------
      call esmf_logwrite("set up the esmf time",                        &
                         esmf_log_info, rc = rc1)

      call esmf_timeset(starttime, yy = yy, mm = mm,  dd = dd,          &
                              h  = hh, m  = mns, s  = sec, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'set up the esmf time',rc)

      call esmf_logwrite("set up the reference time interval",          &
                        esmf_log_info, rc = rc1)

      call esmf_timeintervalset(reftimeinterval, h = int_state%kfhour,  &
                           m = 0, rc = rc1)

! re-set up the start time based on the kfhour value in the sigma file.
!----------------------------------------------------------------------
!      starttime = starttime + reftimeinterval

!     call gfs_dynamics_err_msg(rc1,					&
!                        'set up the reference time interval',rc)

! set the esmf clock which will control the gfs run do loop.
!--------------------------------------------------------------

      currtime = starttime + reftimeinterval
      call esmf_clockset(clock, currtime = currtime,                    &
                         rc = rc1)
!
! get the grid component vm.
! this esmf_gridcompget vm can be used at any where you need it.
!---------------------------------------------------------------
!      call esmf_logwrite("get the local vm", esmf_log_info, rc = rc1)

      call esmf_vmgetcurrent(vm_local, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'get the vm',rc)


! set up parameters of mpi communications.
! use esmf utility to get pe identification and total number of pes.
!-------------------------------------------------------------------
!      call esmf_logwrite("get me and nodes from vm", 			&
!                          esmf_log_info, rc = rc1)

      call esmf_vmget(vm_local, localpet = int_state%me,    		&
                           mpicommunicator = mpi_comm_all,      	&
                           petcount = int_state%nodes,			&
                           rc       = rc1)

      call gfs_dynamics_err_msg(rc1,'get me and nodes from vm',rc)
!      write(0,*)'in dyn_gc,after vmget,npes=',int_state%nodes,'mpi_comm_all=',mpi_comm_all

! initialize the gfs, including set up the internal state
! variables and some local parameter short names, aloocate
! internal state arrays.
!---------------------------------------------------------
      call esmf_logwrite("run the gfs_dynamics_initialize", 		&
                         esmf_log_info, rc = rc1)

! ======================================================================
! ----------------- gfs dynamics related initialize --------------------
! ======================================================================
! grid_gr unfolded (sarah lu)
      call gfs_dynamics_initialize(int_state, rc1)
!      write(0,*)'in dyn_init, t=',maxval(int_state%grid_gr(:,int_state%g_t)), &
!       minval(int_state%grid_gr(:,int_state%g_t)),'quilting=',quilting
! ======================================================================
! ----------------------------------------------------------------------
! ======================================================================

      call gfs_dynamics_err_msg(rc1,'run the gfs_dynamics_initialize',rc)

      call esmf_clockget(clock, timestep    = timestep,            	&
                         runduration = runduration,              	&
                         starttime   = starttime,                	&
                         currtime    = currtime,                 	&
                         rc          = rc1)
!
!
      call esmf_timeintervalget(runduration,                            &
                                h = runduration_hour, rc = rc1)
!
!
!moor ifhmax = nint(int_state%nam_gfs_dyn%fhmax)
      ifhmax = nint(fhmax)
      if(runduration_hour <= 0    .or.                  		&
          ifhmax /= 0             .and.                 		&
          ifhmax <= int_state%kfhour + runduration_hour) then
          ifhmax            = nint(fhmax)
          runduration_hour  = nint(fhmax) - nint(fhini)
          call esmf_timeintervalset(runduration,                        &
                                    h = runduration_hour, rc = rc1)
      end if
      if (runduration_hour <= 0) then
        write(0,*)'WRONG: fhini=',fhini, ' >= fhmax=',fhmax,' job aborted'
        if(me.eq.0)  call mpi_quit(444)
      endif
      stoptime = currtime  + runduration
                           
      call esmf_clockset(clock, stoptime = stoptime,               	&
                         rc       = rc1)
!
      call esmf_timeintervalget(timestep, s = timestep_sec, rc = rc1)
                           
!!
      if (me.eq.0) then
        call out_para(real(timestep_sec))
      endif
!!
      if (me.eq.0) then
        print *,' the gsm will forecast ',runduration_hour,' hours',    &
                ' from hour ',int_state%kfhour,' to hour ',               &
                 runduration_hour+int_state%kfhour
      endif
!
!
      call synchro
!
! create 3D Gaussian grid  (sarah lu)
!-----------------------
!
      call gfs_dynamics_grid_create_Gauss3D(vm_local,int_state,DistGrid5,rc1)

      call gfs_dynamics_err_msg(rc1,'gfs_dynamics_grid_create_gauss3d',rc)

      int_state%fhour_idate(1,1)=fhour
      int_state%fhour_idate(1,2:5)=idate(1:4)
!
      IF(int_state%ENS) THEN
          int_state%end_step = .true.

          CALL AddParameterToState(exp_gfs_dyn, int_state, rc = rc1)

          call gfs_dynamics_err_msg(rc1,                                    &
               'Add Parameter To export State',rc)
      END IF

!
! Define Dynamics Export states    (Sarah Lu)
!
      call gfs_dynamics_internal2export(int_state, exp_gfs_dyn, rc1)

      call gfs_dynamics_err_msg(rc1,'gfs_dynamics_internal2export',rc)
!
!-------------------------------------------------------
! send all the head info to write tasks
!-------------------------------------------------------
!
      call point_dynamics_output_gfs(int_state,IMP_STATE_WRITE)

!
!*******************************************************************
! print out the final error signal variable and put it to rc.
!------------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,				&
                        'initialize from gfs dynamics grid comp.',rc)

      end subroutine gfs_dyn_initialize





!----------------------------------------------------------------------
!bop
!
! !routine: gfs_dyn_run --- 
!           main grid component routine to run the gfs dynamics.
!
! !description: this subroutine will run the most part computations 
!               of the gfs dynamics.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  july     2007     hann-ming henry juang
!  oct 12 2009       Sarah Lu, point grid_gr to import state once and for all
!  oct 17 2009       Sarah Lu, debug print added to track imp/exp states
!  November 2009     Weiyu Yang, Ensemble GEFS.
!  nov 09 2009       Jun Wang, get data from grid_gr_dfi to internal state for dfi
!  feb 05 2010       Jun Wang, set restart step
!
! !interface:
!

      subroutine gfs_dyn_run(gc_gfs_dyn,			          &
                            imp_gfs_dyn, exp_gfs_dyn, clock, rc)

      use gfs_dynamics_states_mod
      use gfs_dyn_date_def
!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp), intent(inout) :: gc_gfs_dyn   
      type(esmf_state),    intent(inout) :: imp_gfs_dyn 
 
! !output variables and parameters:
!----------------------------------
      type(esmf_clock),    intent(inout) :: clock
      type(esmf_timeinterval)            :: timestep, donetime    
      type(esmf_time)                    :: starttime    
      type(esmf_time)                    :: currtime     
      type(esmf_time)                    :: stoptime     
      type(esmf_time)                    :: dfitime     
      type(esmf_state),    intent(inout) :: exp_gfs_dyn
      integer,             intent(out)   :: rc   
!eop
!-------------------------------------------------------------------------

!
! !working arrays and local parameters.
!--------------------------------------
      type(gfs_dyn_wrap)                :: wrap         
! this wrap is a derived type which contains
! only a pointer to the internal state.  it is needed
! for using different architectures or compliers.
      type(gfs_dynamics_internal_state), pointer  :: int_state   
      integer                                     :: rc1          
      integer                                     :: rcfinal     
!
      type(esmf_state)                  :: imp_state_write  !<-- The write gc import state
      logical,save                           :: first_reset=.true.
      logical,save                           :: first_dfiend=.true.
      TYPE(ESMF_TimeInterval)                :: HALFDFIINTVAL
      integer                                :: DFIHR
#ifdef ESMF_3
      TYPE(ESMF_LOGICAL)                     :: Cpl_flag1  ! in ESMF 3.1.0rp2, to  get logical from state
                                                       ! must use the ESMF_LOGICAL type.
#endif

!! debug print for tracking import and export state (Sarah Lu)
      TYPE(ESMF_Field)                   :: ESMFField             !chlu_debug
      TYPE(ESMF_FieldBundle)             :: ESMFBundle            !chlu_debug
      REAL , DIMENSION(:,:,:), POINTER   :: fArr3D                !chlu_debug
      integer                            :: localPE,ii1,ii2,ii3   !chlu_debug
      integer                            :: n, k, rc2             !chlu_debug
      logical, parameter                 :: ckprnt = .false.      !chlu_debug
      integer, parameter                 :: item_count = 3        !chlu_debug
      character(5) :: item_name(item_count)                       !chlu_debug
      character(20) :: vname                                      !chlu_debug
      data item_name/'t','u','v'/                                 !chlu_debug

      localPE = 0                                                 !chlu_debug

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! retrieve the esmf internal state.
!---------------------------------- 
      call esmf_logwrite("get the internal state in the run routine", 	&
                        esmf_log_info, rc = rc1)

      call esmf_gridcompgetinternalstate(gc_gfs_dyn, wrap, rc1)

      call gfs_dynamics_err_msg(rc1,					&
                  'get the internal state in the run routine',rc)

! pointing the local internal state pointer to the esmf internal state pointer.
!------------------------------------------------------------------------------
      int_state => wrap%int_state
!
! get the esmf import state and over-write the gfs internal state.
! update the initial condition arrays in the internal state based on
! the information of the esmf import state.
!------------------------------------------------------------------
      call esmf_logwrite("esmf import state to internal state",         &
                        esmf_log_info, rc = rc1)
!
      int_state%reset_step = .false.
      if(int_state%restart_step ) first_reset=.false.
!        int_state%dfiend_step=.false.
      if( int_state%ndfi>0 .and. first_reset.and. int_state%kdt==int_state%ndfi) then
        if( first_dfiend ) then
! first go through dfi step
          int_state%dfiend_step=.true.
          first_dfiend=.false.
        else
! second go through reset step
          int_state%reset_step = .true.
          int_state%dfiend_step = .false.
          first_reset=.false.
        endif
      endif
!      print *,'in grid comp,ndfi=',int_state%ndfi,'kdt=',int_state%kdt,  &
!       'ndfi=',int_state%ndfi,'first_reset=',first_reset,'reset_step=',  &
!       int_state%reset_step,'dfiend_step=',int_state%dfiend_step
!
      IF(.NOT. int_state%restart_step .AND. .NOT. int_state%start_step ) THEN
          IF(.NOT. int_state%reset_step) THEN
              CALL gfs_dynamics_import2internal(imp_gfs_dyn, int_state, rc1)
          ELSE
              CALL gfs_dynamics_import2internal(imp_gfs_dyn,          &
                  int_state, rc = rc1, exp_gfs_dyn = exp_gfs_dyn)

          END IF 

          CALL gfs_dynamics_err_msg(rc1, 'esmf import state to internal state', rc)
          idate(1 : 4) = int_state%fhour_idate(1, 2 : 5)
      END IF
!
! get clock times
! ------------------
      call esmf_clockget(clock,            				&
                         timestep    = timestep,                	&
                         starttime   = starttime,                 	&
                         currtime    = currtime,                 	&
                         stoptime    = stoptime,                	&
                         rc          = rc1)

      call gfs_dynamics_err_msg(rc1,'esmf clockget',rc)

      donetime = currtime-starttime

      int_state%kdt = nint(donetime/timeStep) 

! Set up the ensemble coupling time flag.
!----------------------------------------
#ifdef ESMF_3

      CALL ESMF_AttributeGet(imp_gfs_dyn, 'Cpl_flag', Cpl_flag1, rc = rc1)
      IF(Cpl_flag1 == ESMF_TRUE) THEN
          int_state%Cpl_flag = .true.
      ELSE
          int_state%Cpl_flag = .false.
      END IF

#else

      CALL ESMF_AttributeGet(imp_gfs_dyn, 'Cpl_flag', int_state%Cpl_flag, rc = rc1)

#endif

      if( currtime .eq. stoptime ) then
          print *,' currtime equals to stoptime '
          int_state%end_step = .true.
      else
          int_state%end_step=.false.
      endif
!
! get nfcstdate
      call esmf_timeget(currtime,                                         &
                        yy=int_state%nfcstdate7(1),                       &
                        mm=int_state%nfcstdate7(2),                       &
                        dd=int_state%nfcstdate7(3),                       &
                        h =int_state%nfcstdate7(4),                       &
                        m =int_state%nfcstdate7(5),                       &
                        s =int_state%nfcstdate7(6),                       &
                        rc=rc1)
      call gfs_dynamics_err_msg(rc1,'esmf timeget',rc)
!      print *,'in gfs grid comp,currtime=',int_state%nfcstdate7(1:6)
!
! ======================================================================
! --------------- run the gfs dynamics related -------------------------
! ======================================================================
      call esmf_logwrite("run the gfs_dynamics_run", 			&
                         esmf_log_info, rc = rc1)

      call gfs_dynamics_run(int_state, imp_gfs_dyn, rc = rc1)
      call gfs_dynamics_err_msg(rc1,'run the gfs_dynamics_run',rc)
! ======================================================================
! ======================================================================

! transfer the gfs export fields in the internal state 
! to the esmf export state which is the public interface
! for other esmf grid components. link is done in initialize, so do not need.
!-------------------------------------------------------
     call esmf_logwrite("internal state to esmf export state", 	&
                       esmf_log_info, rc = rc1)

! Need to check it?  Should be removed?
!--------------------------------------
     call gfs_dynamics_internal2export(int_state, exp_gfs_dyn, rc = rc1)

     call gfs_dynamics_err_msg(rc1,'internal state to esmf export state',rc)

!! debug print starts here  (Sarah Lu) -----------------------------------
      lab_if_ckprnt_ex : if ( ckprnt .and. (int_state%me ==0) ) then      !chlu_debug
        do n = 1, item_count                                              !chlu_debug
            vname = trim(item_name(n))                                    !chlu_debug
            if(associated(fArr3D)) nullify(fArr3D)                        !chlu_debug
            CALL ESMF_StateGet(state = exp_gfs_dyn                      & !chlu_debug
                        ,itemName  = vname                              & !chlu_debug
                        ,field     = ESMFField                          & !chlu_debug
                        ,rc        = rc1)                                 !chlu_debug
            call gfs_dynamics_err_msg(rc1,'LU_DYN: get ESMFarray',rc)     !chlu_debug

#ifdef ESMF_3
            CALL ESMF_FieldGet(field=ESMFField, localDe=0, &              !chlu_debug
                               farray=fArr3D, rc = rc1)                   !chlu_debug
#else
            CALL ESMF_FieldGet(field=ESMFField, localDe=0, &              !chlu_debug
                               farrayPtr=fArr3D, rc = rc1)                !chlu_debug
#endif

            call gfs_dynamics_err_msg(rc1,'LU_DYN: get F90array',rc)      !chlu_debug
!            ii1 = size(fArr3D, dim=1)                                     !chlu_debug
!            ii2 = size(fArr3D, dim=2)                                     !chlu_debug
!            ii3 = size(fArr3D, dim=3)                                     !chlu_debug
!            if(n==1) print *, 'LU_DYN:',ii1, 'x', ii2, 'x', ii3           !chlu_debug
!            print *,' LU_DYN: exp_: ',vname,fArr3D(1,1,1),fArr3D(1,2,1),& !chlu_debug
!                         fArr3D(2,1,1),fArr3D(ii1,ii2,ii3)                !chlu_debug
        enddo                                                             !chlu_debug

        call ESMF_StateGet(state=exp_gfs_dyn, ItemName='tracers', &       !chlu_debug
                          fieldbundle=ESMFBundle, rc = rc1 )              !chlu_debug
        call gfs_dynamics_err_msg(rc1,'LU_DYN: get Bundle from exp',rc)   !chlu_debug
        do n = 1, int_state%ntrac                                         !chlu_debug
          vname = int_state%gfs_dyn_tracer%vname(n, 1)                    !chlu_debug
          print *,'LU_DYN:',trim(vname)                                   !chlu_debug
          CALL ESMF_FieldBundleGet(bundle=ESMFBundle, &                   !chlu_debug
                 name=vname, field=ESMFfield, rc = rc1)                   !chlu_debug

#ifdef ESMF_3
          CALL ESMF_FieldGet(field=ESMFField, localDe=0, &                !chlu_debug
                             farray=fArr3D, rc = rc1)                     !chlu_debug
#else
          CALL ESMF_FieldGet(field=ESMFField, localDe=0, &                !chlu_debug
                             farrayPtr=fArr3D, rc = rc1)                  !chlu_debug
#endif

!          if(n==1) then                                                   !chlu_debug
!            ii1 = size(fArr3D, dim=1)                                     !chlu_debug
!            ii2 = size(fArr3D, dim=2)                                     !chlu_debug
!            ii3 = size(fArr3D, dim=3)                                     !chlu_debug
!            print *,'LU_DYN:',ii1, 'x', ii2, 'x', ii3                     !chlu_debug
!          endif                                                           !chlu_debug
!          print *,'LU_DYN: exp_:',trim(vname), &                          !chlu_debug
!               fArr3D(1,1,1),fArr3D(1,2,1),    &                          !chlu_debug
!               fArr3D(2,1,1),fArr3D(ii1,ii2,ii3)                          !chlu_debug
        enddo                                                             !chlu_debug

      endif lab_if_ckprnt_ex                                              !chlu_debug
!! -------------------------------------- debug print ends here  (Sarah Lu)

! ======================================================================
!------------- put run level variables into write_imp_state---------
! ======================================================================
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
     call esmf_logwrite("get imp_state_write from esmf export state",   &
                       esmf_log_info, rc = rc1)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state          =exp_gfs_dyn    &  !<-- The Dyn component's export state
                        ,itemName       ="Write Import State"     &  !<-- Name of state to be extracted
                        ,nestedState    =IMP_STATE_WRITE  &  !<-- The extracted state
                        ,rc             =RC1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
     call gfs_dynamics_err_msg(rc1,'get imp_state_write from esmf export state',rc)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
     call esmf_logwrite("set pdryini in imp_state_write",               &
                       esmf_log_info, rc = rc1)
     CALL ESMF_AttributeSet(state    =IMP_STATE_WRITE                   &  !<-- The Write component import state
                            ,name     ='pdryini'                        &  !<-- Name of the var
                            ,value    =int_state%pdryini                &  !<-- The var being inserted into the import state
                            ,rc       =RC1)
     call gfs_dynamics_err_msg(rc1,'set pdryini in imp_state_write',rc)
!*******************************************************************
!
! print out the final error signal information and put it to rc.
!---------------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,				&
                        'run from gfs dynamics grid comp.',rc)
!      print *,'end of gfs_dyn_run'
      end subroutine gfs_dyn_run


!----------------------------------------------------------------------
!bop
!
! !routine: finalize --- finalizing routine to finish the 
!                        gfs running job.
!
! !description: this subroutine will finish the gfs computations,
! !             and will release the memory space.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  february 2007     juang for dynamics only
!
! !interface:

      subroutine gfs_dyn_finalize(gc_gfs_dyn,		                &
                                 imp_gfs_dyn, exp_gfs_dyn, clock, rc)

!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp), intent(inout)  :: gc_gfs_dyn
      type(esmf_state),    intent(inout)  :: imp_gfs_dyn
      type(esmf_state),    intent(inout)  :: exp_gfs_dyn
      type(esmf_clock),    intent(inout)  :: clock

! !output variables and parameters:
!----------------------------------
      integer,             intent(out)    :: rc

! !working arrays and local parameters.
!--------------------------------------
      type(gfs_dyn_wrap)                            :: wrap   
      type(gfs_dynamics_internal_state), pointer    :: int_state  
      integer                                       :: rc1        
      integer                                       :: rcfinal   

!eop
!-------------------------------------------------------------------------

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! retrieve the esmf internal state.
!----------------------------------
     call esmf_logwrite(						&
                      "get the internal state in the finalize routine", &
                       esmf_log_info, rc = rc1)

     call esmf_gridcompgetinternalstate(gc_gfs_dyn, wrap, rc1)

     call gfs_dynamics_err_msg(rc1,					&
              'get the internal state in the finalize routine',rc)

! point the local internal state pointer to the esmf internal state pointer.
!------------------------------------------------------------------------------
      int_state => wrap%int_state

! ======================================================================
! run the gfs finalize routine to release the memory space, etc. 
! ======================================================================
      call esmf_logwrite("run the gfs_dynamics_finalize", 		&
                         esmf_log_info, rc = rc1)

      call gfs_dynamics_finalize(int_state, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'run the gfs_dynamics_finalize',rc)
! ======================================================================
! ======================================================================

! print out the final error signal information and put it to rc.
!---------------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,				&
                        'finalize from gfs dynamics grid comp.',rc)

      end subroutine gfs_dyn_finalize

! end of the gfs esmf grid component module.
!-------------------------------------------
      end module gfs_dynamics_grid_comp_mod
