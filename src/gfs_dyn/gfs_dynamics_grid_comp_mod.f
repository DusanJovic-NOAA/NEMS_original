! !module: gfs_dynamics_grid_comp_mod --- 
!                       esmf gridded component of gfs dynamics
!
! !description: gfs dynamics gridded component main module.
!
! !revision history:
!
!  january 2007     hann-ming henry juang
!                           
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
!jws
      use gfs_dyn_mpi_def
      use gfs_dynamics_output, only : point_dynamics_output_gfs
!jwe
      implicit none

#include "../../inc/ESMF_LogMacros.inc"

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
                                       gfs_dyn_initialize,    		&
                                       esmf_singlephase, rc1)
      call gfs_dynamics_err_msg(rc1,'set entry point for initialize',rc)

! register the run subroutine.
!-----------------------------
      call esmf_logwrite("set entry point for run",              	&
                           esmf_log_info, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_dyn, 			&
                                       esmf_setrun,   			&
                                       gfs_dyn_run,           		&
                                       esmf_singlephase, rc1)
      call gfs_dynamics_err_msg(rc1,'set entry point for run',rc)


! register the finalize subroutine.
!----------------------------------
      call esmf_logwrite("set entry point for finalize",                &
                        esmf_log_info, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_dyn, 			&
                                       esmf_setfinal, 			&
                                       gfs_dyn_finalize,       		&
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
!
! !interface:
!

! this argument list is a standard list for all the initialize,
! the run and finalize routines for an esmf system.
!--------------------------------------------------------------
      subroutine gfs_dyn_initialize(gc_gfs_dyn, 			&
                                   imp_gfs_dyn, exp_gfs_dyn, clock, rc)

! user code, for computations related to the esmf interface states.
!------------------------------------------------------------------
      use gfs_dynamics_states_mod
      use gfs_dynamics_grid_create_mod
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

      INTEGER , DIMENSION(:, :), POINTER :: i2

      TYPE(ESMF_DistGrid)                :: DistGrid0    ! the ESMF DistGrid.
      TYPE(ESMF_DistGrid)                :: DistGrid1    ! the ESMF DistGrid.
      TYPE(ESMF_DistGrid)                :: DistGrid2    ! the ESMF DistGrid.
      TYPE(ESMF_DistGrid)                :: DistGrid3    ! the ESMF DistGrid.
      TYPE(ESMF_DistGrid)                :: DistGrid4    ! the ESMF DistGrid.

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
      call esmf_logwrite("allocate the internal state",                 &
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
!
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

! initialize time interval to the parameter from the configure file.
!-------------------------------------------------------------------
!     call esmf_logwrite("set up time step interval",                   &
!                       esmf_log_info, rc = rc1)

!     call esmf_clockget(clock,            				&
!                        timestep    = timestep,                	&
!                        rc          = rc1)

!     call esmf_timeintervalget(timestep,                               &
!                               s  = atm_timestep_s,                    &
!                               rc = rc1)

!     dyn_timestep_s = nint(int_state%nam_gfs_dyn%deltim)

!     int_state%nam_gfs_dyn%deltim = atm_timestep_s /			&
!           min( 1, atm_timestep_s / dyn_timestep_s  )

!     call gfs_dynamics_err_msg(rc1,'set up time step interval',rc)

! get the start time from reading the sigma file.
!----------------------------------------------------------
      call esmf_logwrite("getting the start time",                      &
                         esmf_log_info, rc = rc1)


      call gfs_dynamics_start_time_get(					&
                        yy, mm, dd, hh, mns, sec, int_state%kfhour,     &
                        int_state%n1,int_state%n2,int_state%grib_inp,   &
                        fhrot, int_state%nam_gfs_dyn%sig_ini,           &
                        int_state%nam_gfs_dyn%sig_ini2, rc1)
 
      call gfs_dynamics_err_msg(rc1,'getting the start time',rc)
 
      advancecount4    = nint(real(int_state%kfhour) * 3600.0 /         &
                              int_state%nam_gfs_dyn%deltim)
      int_state%phour  = advancecount4 * 				&
                         int_state%nam_gfs_dyn%deltim / 3600.0
      int_state%kfhour = nint(int_state%phour)

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
      starttime = starttime + reftimeinterval

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
      call esmf_logwrite("get the local vm", esmf_log_info, rc = rc1)

      call esmf_vmgetcurrent(vm_local, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'get the vm',rc)


! set up parameters of mpi communications.
! use esmf utility to get pe identification and total number of pes.
!-------------------------------------------------------------------
      call esmf_logwrite("get me and nodes from vm", 			&
                          esmf_log_info, rc = rc1)

      call esmf_vmget(vm_local, localpet = int_state%me,    		&
                           mpicommunicator = mpi_comm_all,      	&
                           petcount = int_state%nodes,			&
                           rc       = rc1)

      call gfs_dynamics_err_msg(rc1,'get me and nodes from vm',rc)
      write(0,*)'in dyn_gc,after vmget,npes=',int_state%nodes,'mpi_comm_all=',mpi_comm_all

! Allocate the local index array i2 to store the local size information of the
! ditributed grid1, grid3, etc..  Information is based per dimension and per De.
!-------------------------------------------------------------------------------
      ALLOCATE(i2(2, Int_State%NODES))

! initialize the gfs, including set up the internal state
! variables and some local parameter short names, aloocate
! internal state arrays.
!---------------------------------------------------------
      call esmf_logwrite("run the gfs_dynamics_initialize", 		&
                         esmf_log_info, rc = rc1)

! ======================================================================
! ----------------- gfs dynamics related initialize --------------------
! ======================================================================
      call gfs_dynamics_initialize(int_state, rc1)
      write(0,*)'in dyn_init, t=',maxval(int_state%grid_gr(:,int_state%g_t)), &
       minval(int_state%grid_gr(:,int_state%g_t)),'quilting=',quilting
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
      if (runduration_hour < 0) then
        print *,' fhini=',fhini, ' > fhmax=',fhmax,' job aborted'
        call mpi_quit(444)
      endif
      stoptime = currtime  + runduration
                           
      call esmf_clockset(clock, stoptime = stoptime,               	&
                         rc       = rc1)
!
      call esmf_timeintervalget(timestep, s = timestep_sec, rc = rc1)
                           
      print *,' timestep_sec=',timestep_sec,' rc1=',rc1
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
! create the esmf grids. 
!-----------------------
      call esmf_logwrite("creat the esmf grid and delayout.", 		&
                         esmf_log_info, rc = rc1)

      call gfs_dynamics_grid_create_spect(vm_local,int_state, &
                                          DistGrid0, DistGrid1, DistGrid2, rc1)

      call gfs_dynamics_err_msg(rc1,'gfs_dynamics_grid_create_spect',rc)

      call gfs_dynamics_grid_create_gauss(vm_local,int_state, &
                                          DistGrid0, DistGrid3, DistGrid4, rc1)

      call gfs_dynamics_err_msg(rc1,'gfs_dynamics_grid_create_gauss',rc)


! associate the grid3 with the esmf grid component gsgfs
! used at the begining of the run routine when read in
! the surface arrays of the esmf import state.
!-------------------------------------------------------
      call esmf_logwrite(						&
                   "attach the esmf grids to the esmf grid component.", &
                   esmf_log_info, rc = rc1)

      call esmf_gridcompset(gc_gfs_dyn, grid = grid3, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'esmf_gridcompset - set grid3',rc)

! get the local array size of the grid1, the single level spectral arrays.
!-------------------------------------------------------------------------
      i2 = 0
      CALL ESMF_DistGridGet(DistGrid1, indexCountPDimPDe = i2, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'get grid1 info for lnt2_s ',rc)

! put the grid1 local array size into the internal state and print it out.
!-------------------------------------------------------------------------
      Int_State%lnt2_s = i2(1, Int_State%me + 1)

      print*, 'local number of the grid1', i2(1, Int_State%me + 1)

! ======================================================================
! get the local array size of the grid3, the gaussian grid arrays.
!-----------------------------------------------------------------
      i2 = 0
      CALL ESMF_DistGridGet(DistGrid3, indexCountPDimPDe = i2, rc = rc1)

! put the grid3 local array size into the internal state and print it out.
!-------------------------------------------------------------------------
      int_state%llgg_s = i2(1, Int_State%me + 1)

      print*, 'local number of the grid3', i2(:, Int_State%me + 1)

      call gfs_dynamics_err_msg(rc1,'grid get info - llgg_s',rc)

! get the size of grid0.  it is just for testing and can be removed.
!-------------------------------------------------------------------
!     i2 = 0
!     CALL ESMF_DistGridGet(DistGrid0, indexCountPDimPDe = i2, rc = rc1)

!     print*, 'local number of the grid0', i2(:, Int_State%me + 1)

!     call gfs_dynamics_err_msg(rc1,'grid get info - grid0',rc)

! allocate points of import/export state in internal module
!-------------------------------------------------------
!     call esmf_logwrite("allocate internal state for import/export", 	&
!                       esmf_log_info, rc = rc1)

!     call gfs_dynamics_states_allocate(gc_gfs_dyn, int_state, 	        &
!                                        rc = rc1)

!     call gfs_dynamics_err_msg(rc1,                                    &
!          'allocate internal state to import/export',rc)

! transfer the gfs fields in the internal state 
! to the esmf export state which is the public interface
! for other esmf grid components.
!-------------------------------------------------------
!     call esmf_logwrite("transfor internal state to export state ", 	&
!                       esmf_log_info, rc = rc1)

      int_state%fhour_idate(1,1)=fhour
      int_state%fhour_idate(1,2:5)=idate(1:4)
!
!      call gfs_dynamics_internal2export(gc_gfs_dyn, int_state,         &
!                                       exp_gfs_dyn, rc1)

!     call gfs_dynamics_err_msg(rc1,                                    &
!          'transfor internal state to export state',rc)

      DEALLOCATE(i2)
!
!-------------------------------------------------------
!##jw send all the head info to write tasks
!-------------------------------------------------------
!
      if(quilting) then
        write(0,*)'before point_dynamics_output_gfs'
        call point_dynamics_output_gfs(int_state,IMP_STATE_WRITE)
      endif
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
!
! !interface:
!

      subroutine gfs_dyn_run(gc_gfs_dyn, 				&
                            imp_gfs_dyn, exp_gfs_dyn, clock, rc)

      use gfs_dynamics_states_mod
!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp), intent(inout) :: gc_gfs_dyn   
      type(esmf_state),    intent(in)    :: imp_gfs_dyn 
 
! !output variables and parameters:
!----------------------------------
      type(esmf_clock),    intent(inout) :: clock
      type(esmf_timeinterval)            :: timestep, donetime    
      type(esmf_time)                    :: starttime    
      type(esmf_time)                    :: currtime     
      type(esmf_time)                    :: stoptime     
      type(esmf_state),    intent(inout) :: exp_gfs_dyn
      integer,             intent(out)   :: rc   
!
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

! get the esmf import state and over-write the gfs internal state.
! update the initial condition arrays in the internal state based on
! the information of the esmf import state. 
!------------------------------------------------------------------
      call esmf_logwrite("esmf import state to internal state", 	&
                        esmf_log_info, rc = rc1)

      if( .not. int_state%start_step ) then
        call gfs_dynamics_import2internal(gc_gfs_dyn, imp_gfs_dyn, 	&
                                          int_state, rc1)
        idate(1:4)=int_state%fhour_idate(1,2:5)
      endif

      call gfs_dynamics_err_msg(rc1,'esmf import state to internal state',rc)

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

      if( currtime .eq. stoptime ) then
          print *,' currtime equals to stoptime '
          int_state%end_step=.true.
      else
          int_state%end_step=.false.
      endif

! ======================================================================
! --------------- run the gfs dynamics related -------------------------
! ======================================================================
      call esmf_logwrite("run the gfs_dynamics_run", 			&
                         esmf_log_info, rc = rc1)

      call gfs_dynamics_run(int_state, rc = rc1)

      call gfs_dynamics_err_msg(rc1,'run the gfs_dynamics_run',rc)
! ======================================================================
! ======================================================================

! transfer the gfs export fields in the internal state 
! to the esmf export state which is the public interface
! for other esmf grid components. link is done in initialize, so do not need.
!-------------------------------------------------------
     call esmf_logwrite("internal state to esmf export state", 	&
                       esmf_log_info, rc = rc1)

     call gfs_dynamics_internal2export(gc_gfs_dyn, int_state,  	&
                                         exp_gfs_dyn, rc1)

     call gfs_dynamics_err_msg(rc1,'internal state to esmf export state',rc)
 
!*******************************************************************
!
! print out the final error signal information and put it to rc.
!---------------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,				&
                        'run from gfs dynamics grid comp.',rc)

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

      subroutine gfs_dyn_finalize(gc_gfs_dyn, 				&
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
