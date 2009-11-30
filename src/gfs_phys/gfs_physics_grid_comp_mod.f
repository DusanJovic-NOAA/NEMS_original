! !module: gfs_physics_grid_comp_mod --- 
!                       esmf gridded component of gfs physics
!
! !description: gfs physics gridded component main module.
!
! !revision history:
!
!  july     2007     shrinivas moorthi
!  november 2007     hann-ming henry juang 
!  may      2009     jun wang, add write grid component
!  october  2009     jun wang, output every time step, add no quilting option 
!  oct 09   2009     sarah lu, 3D Gaussian grid (DistGrid5) added
!  oct 12   2009     sarah lu, set the association between imp/exp states
!                    and internal state grid_fld; reset start_step
!  oct 17 2009      Sarah Lu, add debug print to check imp/exp state
!                           
!
! !interface:
!
      module gfs_physics_grid_comp_mod
 
!!uses:
!------
      use esmf_mod

      use gfs_physics_err_msg_mod,        ONLY: gfs_physics_err_msg,        &
                                                gfs_physics_err_msg_final
      use gfs_physics_initialize_mod,     ONLY: gfs_physics_initialize
      use gfs_physics_run_mod,            ONLY: gfs_physics_run
      use gfs_physics_finalize_mod,       ONLY: gfs_physics_finalize
      USE gfs_physics_getcf_mod,          ONLY: gfs_physics_getcf
      USE gfs_physics_internal_state_mod, ONLY: gfs_physics_internal_state, &
                                                gfs_phy_wrap
      USE mpi_def,                        ONLY: mpi_comm_all,quilting
      USE layout1,                        ONLY: me
      USE date_def,                       ONLY: idate, fhour
      USE namelist_physics_def,           ONLY: fhini, fhmax
!jw
      USE gfs_physics_output,             ONLY: point_physics_output_gfs
!
      use GFS_Phy_States_Mod,             ONLY: gfs_physics_import2internal_mgrid, &  
                                                gfs_physics_internal2export_mgrid  
!
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
!  may      2009     j. wang
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
!*    use gfs_physics_states_mod
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
!jw
      type(esmf_state)                   :: imp_wrt_state
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
!jws
!-----------------------------------------------------------------------
!***  retrieve the import state of the write gridded component
!***  from the physics export state.
!-----------------------------------------------------------------------
      call esmf_logwrite("Retrieve Write Import State from Physics Export State", &
                        esmf_log_info, rc = rc1)
 
      CALL ESMF_StateGet(state      =exp_gfs_phy                        &  !<-- The Physics export state
                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Physics export state
                        ,nestedState=imp_wrt_state                      &  !<-- Extract Write component import state from Physics export
                        ,rc         =RC1)
!
      CALL gfs_physics_err_msg(rc1,"Retrieve Write Import State from Physics Export State",RC)
!jwe
!
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
                           mpicommunicator = mpi_comm_all,        	&
                           petcount = int_state%nodes,			&
                           rc       = rc1)
      me = int_state%me

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
!      write(0,*)'in after init, size fhour_idate=',size(int_state%fhour_idate,1), &
!        size(int_state%fhour_idate,2),'idate size=',size(idate)
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

!
! create 3D Gaussian grid  (sarah lu)                                          
!-----------------------
!
      call gfs_physics_grid_create_Gauss3D(vm_local,int_state,DistGrid5,rc1) 

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

!      write(0,*)'in grid, comp, size fhour_idate=',size(int_state%fhour_idate,1), &
!        size(int_state%fhour_idate,2),'idate size=',size(idate)
      int_state%fhour_idate(1,1)=int_state%kfhour
      int_state%fhour_idate(1,2:5)=idate(1:4)

!      call gfs_physics_internal2export(gc_gfs_phy, int_state,  	&
!                                          exp_gfs_phy, rc = rc1)

!     call gfs_physics_err_msg(rc1,'internal state to esmf export state',rc)

      DEALLOCATE(i2)
!
!-------------------------------------------------------
!jw send all the head info to write tasks
!-------------------------------------------------------
!
        call point_physics_output_gfs(int_state,imp_wrt_state)
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
!  oct 12 2009       Sarah Lu, call gfs_physics_import2internal_mgrid and
!                    gfs_physics_internal2export_mgrid to associate imp/exp
!                    states with internal state grid_fld
!  oct 17 2009       Sarah Lu, debug print added to track imp/exp states

!
! !interface:
!

      subroutine gfs_phy_run(gc_gfs_phy, 				&
                            imp_gfs_phy, exp_gfs_phy, clock, rc)

!*     use gfs_physics_states_mod
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
!
!jw
      type(esmf_state)                   :: imp_wrt_state
!
! these logic flags are used to handle pointer/copy options (Sarah Lu)
      logical       :: imp2int, flag1, flag2 

!! debug print for tracking import and export state (Sarah Lu)
      TYPE(ESMF_Field)                   :: ESMFField             !chlu_debug
      TYPE(ESMF_FieldBundle)             :: ESMFBundle            !chlu_debug
      REAL , DIMENSION(:,:), POINTER     :: fArr2D                !chlu_debug
      REAL , DIMENSION(:,:,:), POINTER   :: fArr3D                !chlu_debug
      integer                            :: localPE,ii1,ii2,ii3   !chlu_debug
      integer                            :: n, k, rc2             !chlu_debug
      integer                            :: exp_item              !chlu_debug
      logical, parameter                 :: ckprnt = .false.      !chlu_debug
      integer, parameter                 :: item_count = 3        !chlu_debug
      integer, parameter                 :: nfld_2d = 16          !chlu_debug
      character(20)                      :: exp_item_name(50)     !chlu_debug
      character(20)                      :: item_name(item_count) !chlu_debug
      character(8)                       :: vname_2d(nfld_2d)*8   !chlu_debug
      character(20)                      :: vname                 !chlu_debug


      data item_name/'t','u','v'/                                 !chlu_debug

      data vname_2d /'slmsk', 'fice', 'hpbl', 'smc1',     &       !chlu_debug
                     'stype', 'vtype', 'vfrac', 'rain',   &       !chlu_debug
                     'rainc', 'dtsfci', 'tsea', 'stc1',   &       !chlu_debug
                     'u10m', 'v10m',  'ustar','zorl'/             !chlu_debug


      localPE = 0                                                 !chlu_debug
!
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
!
! the pointer/copy option (Sarah Lu)
!  get the esmf import state and over-write the gfs internal state         
!  for one-copy option, import2internal is called once and for all        
!  for two-copy option, import2internal is called every time step        
!
!*    call gfs_physics_import2internal(gc_gfs_phy, imp_gfs_phy, 	&
!*                                        int_state, rc = rc1)
!

      flag1  =  int_state%grid_aldata                                     
      flag2  = .not. int_state%grid_aldata .and. int_state%start_step   
      imp2int = flag1 .or. flag2                                    

      if ( imp2int ) then
        call gfs_physics_import2internal_mgrid( imp_gfs_phy,    &         
                                            int_state, rc = rc1)      
        call gfs_physics_err_msg(rc1,'import2internal_mgrid',rc)     
      endif                                                         

      idate(1:4)=int_state%fhour_idate(1,2:5)
      fhour     =int_state%fhour_idate(1,1)

      call gfs_physics_err_msg(rc1,'esmf import state to internal state',rc)

!! debug print starts here  (Sarah Lu) -----------------------------------
      lab_if_ckprnt_im : if ( ckprnt .and. (int_state%me==0) ) then       !chlu_debug
        do n = 1, item_count                                              !chlu_debug
            vname = trim(item_name(n))                                    !chlu_debug
            if(associated(fArr3D)) nullify(fArr3D)                        !chlu_debug
            CALL ESMF_StateGet(state = imp_gfs_phy                      & !chlu_debug
                        ,itemName  = vname                              & !chlu_debug
                        ,field     = ESMFField                          & !chlu_debug
                        ,rc        = rc1)                                 !chlu_debug
            call gfs_physics_err_msg(rc1,'LU_PHY: get ESMFarray',rc)      !chlu_debug
            CALL ESMF_FieldGet(field=ESMFField, localDe=0, &              !chlu_debug
                               farray=fArr3D, rc = rc1)                   !chlu_debug
            call gfs_physics_err_msg(rc1,'LU_PHY: get F90array',rc)       !chlu_debug
            ii1 = size(fArr3D, dim=1)                                     !chlu_debug
            ii2 = size(fArr3D, dim=2)                                     !chlu_debug
            ii3 = size(fArr3D, dim=3)                                     !chlu_debug
            if(n==1) print *, 'LU_PHY:',ii1, 'x', ii2, 'x', ii3           !chlu_debug
            print *,' LU_PHY: imp_: ',vname,fArr3D(1,1,1),fArr3D(1,2,1),& !chlu_debug
                         fArr3D(2,1,1),fArr3D(ii1,ii2,ii3)                !chlu_debug
        enddo                                                             !chlu_debug

        call ESMF_StateGet(state=imp_gfs_phy, ItemName='tracers', &       !chlu_debug
                         fieldbundle=ESMFBundle, rc = rc1)                !chlu_debug
        call gfs_physics_err_msg(rc1,'LU_PHY: get Bundle from imp',rc)    !chlu_debug
        do n = 1, int_state%ntrac                                         !chlu_debug
          vname = int_state%gfs_phy_tracer%vname(n)                       !chlu_debug
          print *,'LU_PHY: ',trim(vname)                                  !chlu_debug
          CALL ESMF_FieldBundleGet(bundle=ESMFBundle, &                   !chlu_debug
                       name= vname, field=ESMFField, rc = rc1)            !chlu_debug
          CALL ESMF_FieldGet(field=ESMFField, localDe=0, &                !chlu_debug
                            farray=fArr3D, rc = rc1)                      !chlu_debug
          if(n==1) then                                                   !chlu_debug
             ii1 = size(fArr3D, dim=1)                                    !chlu_debug
             ii2 = size(fArr3D, dim=2)                                    !chlu_debug
             ii3 = size(fArr3D, dim=3)                                    !chlu_debug
             print *,'LU_PHY:',ii1, 'x', ii2, 'x', ii3                    !chlu_debug
          endif                                                           !chlu_debug
          print *,'LU_PHY: imp_:',trim(vname),&                           !chlu_debug
                 fArr3D(1,1,1),fArr3D(1,2,1), &                           !chlu_debug
                 fArr3D(2,1,1),fArr3D(ii1,ii2,ii3)                        !chlu_debug
        enddo                                                             !chlu_debug
      endif lab_if_ckprnt_im                                              !chlu_debug
!! -------------------------------------- debug print ends here  (Sarah Lu)


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
!
!-----------------------------------------------------------------------
!***  retrieve the import state of the write gridded component
!***  from the physics export state.
!-----------------------------------------------------------------------
      call esmf_logwrite("Retrieve Write Import State from Physics Export State", &
                        esmf_log_info, rc = rc1)

      CALL ESMF_StateGet(state      =exp_gfs_phy                        &  !<-- The Physics export state
                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Physics export state
                        ,nestedState=imp_wrt_state                      &  !<-- Extract Write component import state from Physics export
                        ,rc         =RC1)
!
      CALL gfs_physics_err_msg(rc1,"Retrieve Write Import State from Physics Export State",RC)
!-----------------------------------------------------------------------
      CALL ESMF_AttributeSet(state    =imp_wrt_state                    &  !<-- The Write component import state
                            ,name     ='zhour'                          &  !<-- Name of the var
                            ,value    =int_state%zhour                  &  !<-- The var being inserted into the import state
                            ,rc       =RC)

!
! run the gfs.
!--------------------------
      call esmf_logwrite("run the gfs_physics_run", 			&
                         esmf_log_info, rc = rc1)

      call gfs_physics_run(int_state, rc = rc1)

      call gfs_physics_err_msg(rc1,'run the gfs_physics_run',rc)

! transfer the gfs export fields in the internal state 
! to the esmf exprot state which is the public interface
! for other esmf grid components.
!-------------------------------------------------------
     call esmf_logwrite("internal state to esmf export state", 	&
                       esmf_log_info, rc = rc1)

! the pointer/copy option (Sarah Lu)
!  point export state to internal state grid_fld            
!  internal2export is called once and all                  

!*   call gfs_physics_internal2export(gc_gfs_phy, int_state,  	&
!*                                       exp_gfs_phy, rc = rc1)

     if ( int_state%start_step ) then                             
       call gfs_physics_internal2export_mgrid( int_state,        & 
                                          exp_gfs_phy, rc = rc1)  
       call gfs_physics_err_msg(rc1,'internal2export_mgrid',rc)  

       int_state%start_step = .false. 
     endif                                                     

     call gfs_physics_err_msg(rc1,'internal state to esmf export state',rc)

!! debug print starts here  (Sarah Lu) -----------------------------------
      lab_if_ckprnt_ex : if ( ckprnt .and. (int_state%me==0) ) then       !chlu_debug

        if (int_state%kdt == 1) then                                      !chlu_debug
          call esmf_stateget(exp_gfs_phy                               &  !chlu_debug
                            ,itemcount = exp_item                      &  !chlu_debug
                            ,itemnamelist = exp_item_name              &  !chlu_debug
                            ,rc   =rc)                                    !chlu_debug

          print *,'LU_PHY: export item count:',exp_item                   !chlu_debug
          print *,'LU_PHY: export item name :',(exp_item_name(n),n=1,exp_item)!chlu_debug
        endif                                                             !chlu_debug

        do n = 1, nfld_2d                                                 !chlu_debug
            if(associated(fArr2D)) nullify(fArr2D)                        !chlu_debug
            vname = trim(vname_2d(n))                                     !chlu_debug
            CALL ESMF_StateGet(state = exp_gfs_phy                      & !chlu_debug
                        ,itemName  = vname                              & !chlu_debug
                        ,field     = ESMFField                          & !chlu_debug
                        ,rc        = rc1)                                 !chlu_debug
            call gfs_physics_err_msg(rc1,'LU_PHY: get ESMFField',rc)      !chlu_debug
            CALL ESMF_FieldGet(field=ESMFfield, localDe=0, &              !chlu_debug
                          farray=fArr2D, rc = rc1)                        !chlu_debug
            call gfs_physics_err_msg(rc1,'LU_PHY: get F90array',rc)       !chlu_debug
            if ( n == 1 ) then                                            !chlu_debug
             ii1 = size(fArr2D, dim=1)                                    !chlu_debug
             ii2 = size(fArr2D, dim=2)                                    !chlu_debug
             print *, 'LU_PHY:',ii1, 'x', ii2                             !chlu_debug
            endif                                                         !chlu_debug
            print *,' LU_PHY: exp_: ',vname,fArr2D(1,1),fArr2D(1,2),    & !chlu_debug
                         fArr2D(2,1),fArr2D(ii1,ii2)                      !chlu_debug
        enddo                                                             !chlu_debug

        do n = 1, item_count                                              !chlu_debug
            vname = trim(item_name(n))                                    !chlu_debug
            if(associated(fArr3D)) nullify(fArr3D)                        !chlu_debug
            CALL ESMF_StateGet(state = exp_gfs_phy                      & !chlu_debug
                        ,itemName  = vname                              & !chlu_debug
                        ,field     = ESMFField                          & !chlu_debug
                        ,rc        = rc1)                                 !chlu_debug
            call gfs_physics_err_msg(rc1,'LU_PHY: get ESMFarray',rc)      !chlu_debug
            CALL ESMF_FieldGet(field=ESMFField, localDe=0, &              !chlu_debug
                               farray=fArr3D, rc = rc1)                   !chlu_debug
            call gfs_physics_err_msg(rc1,'LU_PHY: get F90array',rc)       !chlu_debug
            ii1 = size(fArr3D, dim=1)                                     !chlu_debug
            ii2 = size(fArr3D, dim=2)                                     !chlu_debug
            ii3 = size(fArr3D, dim=3)                                     !chlu_debug
            if(n==1) print *, 'LU_PHY:',ii1, 'x', ii2, 'x', ii3           !chlu_debug
            print *,' LU_PHY: exp_: ',vname,fArr3D(1,1,1),fArr3D(1,2,1),& !chlu_debug
                         fArr3D(2,1,1),fArr3D(ii1,ii2,ii3)                !chlu_debug
        enddo                                                             !chlu_debug

        call ESMF_StateGet(state=exp_gfs_phy, ItemName='tracers', &       !chlu_debug
                          fieldbundle=ESMFBundle, rc = rc1 )              !chlu_debug
        call gfs_physics_err_msg(rc1,'LU_PHY: get Bundle from exp',rc)    !chlu_debug
        do n = 1, int_state%ntrac                                         !chlu_debug
          vname = int_state%gfs_phy_tracer%vname(n)                       !chlu_debug
          print *,'LU_PHY:',trim(vname)                                   !chlu_debug
          CALL ESMF_FieldBundleGet(bundle=ESMFBundle, &                   !chlu_debug
                 name=vname, field=ESMFfield, rc = rc1)                   !chlu_debug
          CALL ESMF_FieldGet(field=ESMFfield, localDe=0, &                !chlu_debug
                          farray=fArr3D, rc = rc1)                        !chlu_debug
          if(n==1) then                                                   !chlu_debug
            ii1 = size(fArr3D, dim=1)                                     !chlu_debug
            ii2 = size(fArr3D, dim=2)                                     !chlu_debug
            ii3 = size(fArr3D, dim=3)                                     !chlu_debug
            print *,'LU_PHY:',ii1, 'x', ii2, 'x', ii3                     !chlu_debug
          endif                                                           !chlu_debug
          print *,'LU_PHY: exp_:',trim(vname), &                          !chlu_debug
               fArr3D(1,1,1),fArr3D(1,2,1),    &                          !chlu_debug
               fArr3D(2,1,1),fArr3D(ii1,ii2,ii3)                          !chlu_debug
        enddo                                                             !chlu_debug

      endif lab_if_ckprnt_ex                                              !chlu_debug
!! -------------------------------------- debug print ends here  (Sarah Lu)

!
!-----------------------------------------------------------------------
!***  retrieve the import state of the write gridded component
!***  from the physics export state.
!-----------------------------------------------------------------------
!      call esmf_logwrite("Retrieve Write Import State from Physics Export State", &
!                        esmf_log_info, rc = rc1)
!
!      CALL ESMF_StateGet(state      =exp_gfs_phy                        &  !<-- The Physics export state
!                        ,itemName   ='Write Import State'               &  !<-- Name of the state to get from Physics export state
!                        ,nestedState=imp_wrt_state                      &  !<-- Extract Write component import state from Physics export
!                        ,rc         =RC1)
!!
!      CALL gfs_physics_err_msg(rc1,"Retrieve Write Import State from Physics Export State",RC)
!-----------------------------------------------------------------------
!      CALL ESMF_AttributeSet(state    =imp_wrt_state                    &  !<-- The Write component import state
!                            ,name     ='zhour'                          &  !<-- Name of the var
!                            ,value    =int_state%zhour                  &  !<-- The var being inserted into the import state
!                            ,rc       =RC)
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
