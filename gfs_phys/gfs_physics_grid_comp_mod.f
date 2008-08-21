! !module: gfs_physics_grid_comp_mod --- 
!                       esmf gridded component of gfs physics
!
! !description: gfs physics gridded component main module.
!
! !revision history:
!
!  july     2007     shrinivas moorthi
!  november 2007     hann-ming henry juang 
!                           
!
! !interface:
!
      module gfs_physics_grid_comp_mod
 
!!uses:
!------
      use esmf_mod

      use gfs_physics_err_msg_mod
      use gfs_physics_initialize_mod
      use gfs_physics_run_mod
      use gfs_physics_finalize_mod

      implicit none

      private   ! by default, data is private to this module

      public gfs_phy_setservices	! only set service is public

!eop
!-------------------------------------------------------------------------


      contains


!----------------------------------------------------------------------
!bop
!
! !routine: gfs_phy_setservices --- 
!           set services for gfs physics gridded component.
! 
! !interface:
!
      subroutine gfs_phy_setservices (gc_gfs_phy, rc)
 
! !arguments:
!------------

      type(esmf_gridcomp), intent(in)  :: gc_gfs_phy 	! gridded component
      integer,             intent(out) :: rc    	! return code
     
! !description: set services (register) for the gfs physics grid component.
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
      call esmf_gridcompsetentrypoint (gc_gfs_phy, 			&
                                       esmf_setinit,  			&
                                       gfs_phy_initialize,    		&
                                       esmf_singlephase, rc1)
      call gfs_physics_err_msg(rc1,'set entry point for initialize',rc)

! register the run subroutine.
!-----------------------------
      call esmf_logwrite("set entry point for run",              	&
                           esmf_log_info, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_phy, 			&
                                       esmf_setrun,   			&
                                       gfs_phy_run,           		&
                                       esmf_singlephase, rc1)
      call gfs_physics_err_msg(rc1,'set entry point for run',rc)


! register the finalize subroutine.
!----------------------------------
      call esmf_logwrite("set entry point for finalize",                &
                        esmf_log_info, rc = rc1)
      call esmf_gridcompsetentrypoint (gc_gfs_phy, 			&
                                       esmf_setfinal, 			&
                                       gfs_phy_finalize,       		&
                                       esmf_singlephase, rc1)
      call gfs_physics_err_msg(rc1,'set entry point for finalize',rc)

! check the error signal variable and print out the result.
!----------------------------------------------------------
      call gfs_physics_err_msg_final(rc1,				&
                        'setservice for gfs physics grid comp.',rc)

      end subroutine gfs_phy_setservices





!----------------------------------------------------------------------
!bop
! !routine:  gfs_phy_initialize --- initialize routine to initialize 
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
      subroutine gfs_phy_initialize(gc_gfs_phy, 			&
                                   imp_gfs_phy, exp_gfs_phy, clock, rc)

! user code, for computations related to the esmf interface states.
!------------------------------------------------------------------
      use gfs_physics_states_mod
      use gfs_physics_grid_create_mod
!
! !input/output variables and parameters:
!----------------------------------------

      type(esmf_gridcomp), intent(inout) :: gc_gfs_phy 
      type(esmf_state),    intent(inout) :: imp_gfs_phy
      type(esmf_state),    intent(inout) :: exp_gfs_phy
      type(esmf_clock),    intent(inout) :: clock

!
! !output variables and parameters:
!----------------------------------

      integer, intent(out) :: rc  

! !eop
!------------------------------------------------------------------------- 
 
! !working arrays and local parameters.  
!--------------------------------------
      type(gfs_phy_wrap)                :: wrap         
! this wrap is a derived type which contains
! only a pointer to the internal state.  it is needed
! for using different architectures or compliers.
      type(gfs_physics_internal_state), pointer  :: int_state    
      type(esmf_vm)                      :: vm_local     
      type(esmf_timeinterval)            :: timestep     
      type(esmf_timeinterval)            :: runduration  
      type(esmf_time)                    :: starttime    
      type(esmf_time)                    :: stoptime    
      type(esmf_time)                    :: currtime     
      type(esmf_timeinterval)            :: reftimeinterval 
      type(esmf_delayout)                :: mydelayout   
      integer(kind=esmf_kind_i4)         :: yy, mm, dd   ! time variables for date
      integer(kind=esmf_kind_i4)         :: hh, mns, sec ! time variables for time
      integer                            :: advancecount4, timestep_sec
      integer                            :: atm_timestep_s, phy_timestep_s
      integer(esmf_kind_i8)              :: advancecount

      INTEGER , DIMENSION(:, :), POINTER :: i2

      TYPE(ESMF_DistGrid)                :: DistGrid0    ! the ESMF DistGrid.
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

      call gfs_physics_err_msg(rc1,' - allocate the internal state',rc)

      wrap%int_state => int_state

! attach internal state to the gfs physics grid component.
!-------------------------------------------------
      call esmf_logwrite("set up the internal state",                   &
                        esmf_log_info, rc = rc1)

      call esmf_gridcompsetinternalstate(gc_gfs_phy, wrap, rc1)

      call gfs_physics_err_msg(rc1,'set up the internal state',rc)

! use esmf utilities to get information from the configuration file.
! the function is similar to reading the namelist in the original gfs.
!---------------------------------------------------------------------
      call esmf_logwrite("getting information from the configure file", &
                        esmf_log_info, rc = rc1)

      call gfs_physics_getcf(gc_gfs_phy, int_state,  rc = rc1)

      call gfs_physics_err_msg(rc1,'get configure file information',rc)

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

!     phy_timestep_s = nint(int_state%nam_gfs_phy%deltim)

!     int_state%nam_gfs_phy%deltim = atm_timestep_s /			&
!           min( 1, atm_timestep_s / phy_timestep_s  )

!     call gfs_physics_err_msg(rc1,'set up time step interval',rc)

! get the start time from reading the surface file.
!----------------------------------------------------------
      call esmf_logwrite("getting the start time",                      &
                         esmf_log_info, rc = rc1)


      call gfs_physics_start_time_get(					&
                        yy, mm, dd, hh, mns, sec, int_state%kfhour,     &
                        int_state%n3, int_state%nam_gfs_phy%sfc_ini,rc1)
 
      call gfs_physics_err_msg(rc1,'getting the start time',rc)
 
      advancecount4    = nint(real(int_state%kfhour) * 3600.0 /         &
                              int_state%nam_gfs_phy%deltim)
      int_state%phour  = advancecount4 * 				&
                         int_state%nam_gfs_phy%deltim / 3600.0
      int_state%kfhour = nint(int_state%phour)

! initialize the clock with the start time based on the information
! from calling starttimeget.
!------------------------------------------
      call esmf_logwrite("set up the esmf time",                        &
                         esmf_log_info, rc = rc1)

! in dynamics, we had already timeset, this is redo, update later.

      call esmf_timeset(starttime, yy = yy, mm = mm,  dd = dd,          &
                              h  = hh, m  = mns, s  = sec, rc = rc1)

      call gfs_physics_err_msg(rc1,'set up the esmf time',rc)

      call esmf_logwrite("set up the reference time interval",          &
                        esmf_log_info, rc = rc1)

      call esmf_timeintervalset(reftimeinterval, h = int_state%kfhour,  &
                           m = 0, rc = rc1)

! re-set up the start time based on the kfhour value in the sigma file.
!----------------------------------------------------------------------
!starttime = starttime + reftimeinterval

      call gfs_physics_err_msg(rc1,					&
                         'set up the reference time interval',rc)

! set the esmf clock which will control the gfs run do loop.
!--------------------------------------------------------------

! in dynamics, clock is set, this is redo, update later

      currtime = starttime + reftimeinterval
      call esmf_clockset(clock, 					&
                         currtime = currtime,                      	&
                         rc = rc1)
!
! get the grid component vm.
! this esmf_gridcompget vm can be used at any where you need it.
!---------------------------------------------------------------
      call esmf_logwrite("get the local vm", esmf_log_info, rc = rc1)

      call esmf_vmgetcurrent(vm_local, rc = rc1)

      call gfs_physics_err_msg(rc1,'get the vm',rc)


! set up parameters of mpi communications.
! use esmf utility to get pe identification and total number of pes.
!-------------------------------------------------------------------
      call esmf_logwrite("get me and nodes from vm", 			&
                          esmf_log_info, rc = rc1)

      call esmf_vmget(vm_local, localpet = int_state%me,    		&
                           mpicommunicator = mpi_comm_all,      	&
                           petcount = int_state%nodes,			&
                           rc       = rc1)

      call gfs_physics_err_msg(rc1,'get me and nodes from vm',rc)

! Allocate the local index array i2 to store the local size information of the
! ditributed grid1, grid3, etc..  Information is based per dimension and per De.
!-------------------------------------------------------------------------------
      ALLOCATE(i2(2, Int_State%NODES))

! initialize the gfs, including set up the internal state
! variables and some local parameter short names, aloocate
! internal state arrays.
!---------------------------------------------------------
      call esmf_logwrite("run the gfs_physics_initialize", 		&
                         esmf_log_info, rc = rc1)

! ----------------- gfs physics related initialize --------------------
! ----------------------------------------------------------------------
      call gfs_physics_initialize(int_state, rc1)
! ----------------------------------------------------------------------

      call gfs_physics_err_msg(rc1,'run the gfs_physics_initialize',rc)

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
!moor ifhmax = nint(int_state%nam_gfs_phy%fhmax)
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
        print *,' gsm physics will forecast ',runduration_hour,' hours',  &
                ' from hour ',int_state%kfhour,' to hour ',               &
                 runduration_hour+int_state%kfhour
      endif
!
!
      call synchro
!
! create the esmf grids and distribute the grids into 
! the esmf delayout (mydelayout).
!-----------------------------------------------------------------
      call esmf_logwrite("creat the esmf grid and delayout.", 		&
                         esmf_log_info, rc = rc1)

      call gfs_physics_grid_create_gauss(vm_local,int_state, &
                                         DistGrid0, DistGrid3, DistGrid4, rc1)

      call gfs_physics_err_msg(rc1,'gfs_physics_grid_create_gauss',rc)


! associate the grid3 with the esmf grid component gsgfs
! used at the begining of the run routine when read in
! the surface arrays of the esmf import state.
!-------------------------------------------------------
      call esmf_logwrite(						&
                   "attach the esmf grids to the esmf grid component.", &
                   esmf_log_info, rc = rc1)

      call esmf_gridcompset(gc_gfs_phy, grid = grid3, rc = rc1)

      call gfs_physics_err_msg(rc1,'esmf_gridcompset - set grid3',rc)

! get the local array size of the grid3, the gaussian grid arrays.
!-----------------------------------------------------------------
      i2 = 0
      CALL ESMF_DistGridGet(DistGrid3, indexCountPDimPDe = i2, rc = rc1)

! put the grid3 local array size into the internal state and print it out.
!-------------------------------------------------------------------------
      int_state%llgg_s = i2(1, Int_State%me + 1)

      print*, 'local number of the grid3', i2(:, Int_State%me + 1)

      call gfs_physics_err_msg(rc1,'grid get info - llgg_s',rc)

! get the size of grid0.  it is just for testing and can be removed.
!-------------------------------------------------------------------
!     i2 = 0
!     CALL ESMF_DistGridGet(DistGrid0, indexCountPDimPDe = i2, rc = rc1)

!     print*, 'local number of the grid0', i2(:, Int_State%me + 1)

!     call gfs_physics_err_msg(rc1,'grid get info - grid0',rc)

! transfer the gfs fields in the internal state 
! to the esmf import state which is the public interface
! for other esmf grid components.
!-------------------------------------------------------
!     call esmf_logwrite("allocate internal state for import/export", 	&
!                       esmf_log_info, rc = rc1)

!     call gfs_physics_states_allocate(gc_gfs_phy, int_state, 	        &
!                                        rc = rc1)

!     call gfs_physics_err_msg(rc1,                                    &
!          'allocate internal state to import/export',rc)
!
! set pointer the gfs export fields in the internal state 
! to the esmf exprot state which is the public interface
! for other esmf grid components.
!-------------------------------------------------------
!     call esmf_logwrite("internal state link to esmf export state", 	&
!                       esmf_log_info, rc = rc1)

!     int_state%fhour_idate(1,1)=int_state%kfhour
!     int_state%fhour_idate(1,2:5)=idate(1:4)

!      call gfs_physics_internal2export(gc_gfs_phy, int_state,  	&
!                                          exp_gfs_phy, rc = rc1)

!     call gfs_physics_err_msg(rc1,'internal state to esmf export state',rc)

      DEALLOCATE(i2)

!
!*******************************************************************
! print out the final error signal variable and put it to rc.
!------------------------------------------------------------
      call gfs_physics_err_msg_final(rcfinal,				&
                        'initialize from gfs physics grid comp.',rc)

      end subroutine gfs_phy_initialize





!----------------------------------------------------------------------
!bop
!
! !routine: gfs_phy_run --- 
!           main grid component routine to run the gfs physics.
!
! !description: this subroutine will run the most part computations 
!               of the gfs physics.
!
! !revision history:
!
!  november 2004     weiyu yang initial code.
!  may      2005     weiyu yang for the updated gfs version.
!  february 2006     moorthi
!  december 2007     juang
!  April    2008     Shrinivas Moorthi
!                    Used ESMF alrarm to set the lsout and lssav
!
! !interface:
!

      subroutine gfs_phy_run(gc_gfs_phy, 				&
                            imp_gfs_phy, exp_gfs_phy, clock, rc)

      use gfs_physics_states_mod
      use module_alarms
!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp), intent(inout) :: gc_gfs_phy   
      type(esmf_state),    intent(in)    :: imp_gfs_phy 
 
! !output variables and parameters:
!----------------------------------
      type(esmf_clock),    intent(inout) :: clock
      type(esmf_timeinterval)            :: timestep, donetime    
      type(esmf_time)                    :: starttime    
      type(esmf_time)                    :: currtime     
      type(esmf_time)                    :: stoptime     
      type(esmf_state),    intent(inout) :: exp_gfs_phy
      integer,             intent(out)   :: rc   
      logical, save :: first
      data first/.true./
!
!eop
!-------------------------------------------------------------------------

!
! !working arrays and local parameters.
!--------------------------------------
      type(gfs_phy_wrap)                :: wrap         
! this wrap is a derived type which contains
! only a pointer to the internal state.  it is needed
! for using different architectures or compliers.
      type(gfs_physics_internal_state), pointer   :: int_state   
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

      call esmf_gridcompgetinternalstate(gc_gfs_phy, wrap, rc1)

      call gfs_physics_err_msg(rc1,					&
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

      call gfs_physics_import2internal(gc_gfs_phy, imp_gfs_phy, 	&
                                          int_state, rc = rc1)
      idate(1:4)=int_state%fhour_idate(1,2:5)
      fhour     =int_state%fhour_idate(1,1)

      call gfs_physics_err_msg(rc1,'esmf import state to internal state',rc)

!
! get clock times
! ------------------
      call esmf_clockget(clock,            				&
                         timestep    = timestep,                	&
                         starttime   = starttime,                 	&
                         currtime    = currtime,                 	&
                         stoptime    = stoptime,                	&
                         rc          = rc1)

      call gfs_physics_err_msg(rc1,'esmf clockget',rc)

      donetime = currtime-starttime
      int_state%kdt = nint(donetime/timeStep) + 1

      print *,' in physics kdt=',int_state%kdt

!     if( currtime .eq. stoptime ) then
!         print *,' currtime equals to stoptime '
!         int_state%end_step=.true.
!     endif

! run the gfs.
!--------------------------
      call esmf_logwrite("run the gfs_physics_run", 			&
                         esmf_log_info, rc = rc1)
!
      if (first) then
        int_state%lsoutd = .true.
        first = .false.
      endif
!
!
!     If half of digital filter is over, then set lssav to flase
!
      if (ESMF_AlarmIsRinging(alarm(1), rc)) then
        int_state%lssavd = .false.
        int_state%lsoutd = .false.
      endif

!     print *,' bef phys_run lsoutd=',int_state%lsoutd  &
!    ,' lssavd=',int_state%lssavd,' kdt=',int_state%kdt
!
!     At the end of  digital filter reset lssavd and lsoutd to true
!
      if (.not. int_state%lssavd) then
        if (ESMF_AlarmIsRinging(alarm(2), rc)) then
          int_state%lssavd = .true.
          int_state%lsoutd = .true.
          print *,' LSSAVD and LSOUTD set to TRUE'
        endif
!     print *,' after IF test lsoutd=',int_state%lsoutd  &
!    ,' lssavd=',int_state%lssavd,' kdt=',int_state%kdt
      endif

      call gfs_physics_run(int_state, rc = rc1)

      call gfs_physics_err_msg(rc1,'run the gfs_physics_run',rc)
!
! transfer the gfs export fields in the internal state 
! to the esmf exprot state which is the public interface
! for other esmf grid components.
!-------------------------------------------------------
     call esmf_logwrite("internal state to esmf export state", 	&
                       esmf_log_info, rc = rc1)

     call gfs_physics_internal2export(gc_gfs_phy, int_state,  	&
                                         exp_gfs_phy, rc = rc1)

     call gfs_physics_err_msg(rc1,'internal state to esmf export state',rc)
!
!*******************************************************************
!
! print out the final error signal information and put it to rc.
!---------------------------------------------------------------
      call gfs_physics_err_msg_final(rcfinal,				&
                        'run from gfs physics grid comp.',rc)

      end subroutine gfs_phy_run


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
!  july     2007     juang for physics only
!
! !interface:

      subroutine gfs_phy_finalize(gc_gfs_phy, 				&
                                 imp_gfs_phy, exp_gfs_phy, clock, rc)

!
! !input variables and parameters:
!---------------------------------
      type(esmf_gridcomp), intent(inout)  :: gc_gfs_phy
      type(esmf_state),    intent(inout)  :: imp_gfs_phy
      type(esmf_state),    intent(inout)  :: exp_gfs_phy
      type(esmf_clock),    intent(inout)  :: clock

! !output variables and parameters:
!----------------------------------
      integer,             intent(out)    :: rc

! !working arrays and local parameters.
!--------------------------------------
      type(gfs_phy_wrap)                            :: wrap   
      type(gfs_physics_internal_state), pointer     :: int_state  
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

     call esmf_gridcompgetinternalstate(gc_gfs_phy, wrap, rc1)

     call gfs_physics_err_msg(rc1,					&
              'get the internal state in the finalize routine',rc)

! point the local internal state pointer to the esmf internal state pointer.
!------------------------------------------------------------------------------
      int_state => wrap%int_state

! run the gfs finalize routine to release the memory space, etc. 
!----------------------------------------------------------------------------
      call esmf_logwrite("run the gfs_physics_finalize", 		&
                         esmf_log_info, rc = rc1)

      call gfs_physics_finalize(int_state, rc = rc1)

      call gfs_physics_err_msg(rc1,'run the gfs_physics_finalize',rc)

! print out the final error signal information and put it to rc.
!---------------------------------------------------------------
      call gfs_physics_err_msg_final(rcfinal,				&
                        'finalize from gfs physics grid comp.',rc)

      end subroutine gfs_phy_finalize

! end of the gfs esmf grid component module.
!-------------------------------------------
      end module gfs_physics_grid_comp_mod
