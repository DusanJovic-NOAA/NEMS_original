
!include "options.h"
!
      module atmos_grid_comp_mod
!
!-----------------------------------------------------------------------
!
!***  this is the main (top-level) gridded component module.
!***  it will set up the gridded and coupler subcomponents
!***  and run their initialize, run, and finalize routines.
!
!  February 2007	Hann-Ming Henry Juang
!			options to run eithor dynamics, physics or couple.
!  February 2008        Weiyu Yang, updated to use the ESMF 3.1.0 library.
!
!-----------------------------------------------------------------------
!
      use esmf_mod
!
      use gfs_dynamics_grid_comp_mod  ,only: gfs_dyn_setservices

      use gfs_physics_grid_comp_mod   ,only: gfs_phy_setservices

      use atmos_dyn_phy_cpl_comp_mod  ,only: atm_cpl_setservices
      use module_include
      use atmos_err_msg_mod
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
#include "esmf_log_macros.inc"
!-----------------------------------------------------------------------
! 
      private
!
      public :: atmos_setservices,atmos_initialize,atmos_run,atmos_finalize
!
!-----------------------------------------------------------------------
!
      type(esmf_gridcomp),save :: gc_gfs_dyn
      type(esmf_gridcomp),save :: gc_gfs_phy
      type(esmf_cplcomp), save :: gc_atm_cpl
!
      type(esmf_state),save :: imp_gfs_dyn,exp_gfs_dyn
      type(esmf_state),save :: imp_gfs_phy,exp_gfs_phy
!
      type(esmf_vm),save :: vm     ! the esmf virtual machine.
!
      integer :: inpes,jnpes  ! mpi tasks in i and j
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!#######################################################################
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine atmos_setservices(gc_atm,rc_reg)
! 
!-----------------------------------------------------------------------
!***  register the atmos gridded component's initialize, run, and finalize
!***  routines.
!-----------------------------------------------------------------------
!
      type(esmf_gridcomp),intent(inout) :: gc_atm ! gridded component
!
      integer,intent(out) :: rc_reg    ! return code for register
!     
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer :: rc
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc    =esmf_success  ! the error signal variable
      rc_reg=esmf_success  ! the error signal variable

!-----------------------------------------------------------------------
!***  register the main initialize subroutine.  since it is just one 
!***  subroutine, use esmf_singlephase.  the second argument is
!***  a pre-defined subroutine type, such as esmf_setinit, esmf_setrun, 
!***  or esmf_setfinal.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("set entry point for main initialize"          &
                        ,esmf_log_info,rc=rc)
!
      call esmf_gridcompsetentrypoint(gc_atm                &  !<-- gridcomp
                                     ,esmf_setinit          &  !<-- subroutinetype
                                     ,atmos_initialize      &  !<-- user's routine
                                     ,esmf_singlephase      &
                                     ,rc)
!
      call atmos_err_msg(rc,'set entry point for main initialize',rc_reg)
!
!-----------------------------------------------------------------------
!***  register the main run subroutine.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("set entry point for main run"                 &
                        ,esmf_log_info,rc=rc)
!
      call esmf_gridcompsetentrypoint(gc_atm             &  !<-- gridcomp
                                     ,esmf_setrun        &  !<-- subroutinetype
                                     ,atmos_run          &  !<-- user's routine
                                     ,esmf_singlephase   &
                                     ,rc)
!
      call atmos_err_msg(rc,'set entry point for main run',rc_reg)
!
!-----------------------------------------------------------------------
!***  register the main finalize subroutine.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("set entry point for main finalize"            &
                        ,esmf_log_info,rc=rc)
!
      call esmf_gridcompsetentrypoint(gc_atm             &  !<-- gridcomp
                                     ,esmf_setfinal      &  !<-- subroutinetype
                                     ,atmos_finalize     &  !<-- user's routine
                                     ,esmf_singlephase   &
                                     ,rc)
!
      call atmos_err_msg(rc,'set entry point for main finalize',rc_reg)
!
!-----------------------------------------------------------------------
!***  check the error signal variable and print out the result.
!-----------------------------------------------------------------------
!
      call atmos_err_msg_final(rc_reg,'main_set_services',rc)

!-----------------------------------------------------------------------
!
      end subroutine atmos_setservices
!
!-----------------------------------------------------------------------
!#######################################################################
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine atmos_initialize(gc_atm,imp_state,exp_state       &
                                ,clock,rc_init)
!
!-----------------------------------------------------------------------
!***  this routine sets up fundamental aspects of the model run.
!-----------------------------------------------------------------------
!
!ifdef NMM_B
!     use module_dm_parallel,only: decomp,setup_servers
!endif
!
!-----------------------------------------------------------------------
!***  argument variables.
!-----------------------------------------------------------------------
!
      type(esmf_gridcomp),intent(inout) :: gc_atm
      type(esmf_state),   intent(inout) :: imp_state,exp_state
      type(esmf_clock),   intent(inout) :: clock
!
      integer,            intent(out)   :: rc_init
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      type(esmf_config)            :: cf           ! the config object
!
      type(esmf_grid)              :: grid_atmos    ! the esmf grid for the integration attached to 
                                                   ! the main gridded component.
      type(esmf_grid)              :: grid_gfs_dyn     ! the esmf grid for the integration attached to 
                                                   ! the dynamics gridded component.
      type(esmf_grid)              :: grid_gfs_phy     ! the esmf grid for the integration attached to 
                                                   ! the physics gridded component.

      TYPE(ESMF_DistGrid)          :: DistGrid_atmos

!
      integer, dimension(2)        :: ncounts      ! parameter array to set up the
                                                   ! size of the 2-d esmf grid.
      integer, dimension(2)        :: min,max      ! parameter arrays to set up the
                                                   ! start number and the end number of
                                                   ! the esmf grid in each dimension.
      real(esmf_kind_r8),dimension(esmf_maxgriddim) :: mincoords,maxcoords
      integer,dimension(esmf_maxgriddim) :: counts
      INTEGER , DIMENSION(2)             :: i1
      INTEGER , DIMENSION(:, :), POINTER :: i2
      integer                      :: im,jm,lm                  ! full grid dimensions
      integer                      :: mype,num_pes,num_pes_fcst,num_pes_tot
      integer                      :: mpi_intra,mpi_intra_b     ! the mpi intra-communicator
!
      integer                      :: rc,irtn
      integer                      :: rcfinal
!
      character(50) :: mode
      logical :: global, physics_on
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  initialize the error signal variables.
!-----------------------------------------------------------------------
!
      rc     =esmf_success
      rcfinal=esmf_success
!
!-----------------------------------------------------------------------
!***  retrieve the config object cf from the main gridded component.
!-----------------------------------------------------------------------
!***  retrieve the vm (virtual machine) of the main gridded component.
!***  call esmf_gridcompget to retrieve the vm anywhere you need it.
!***  we need vm now to set up the de layout.
!-----------------------------------------------------------------------
!
      call esmf_logwrite("retrieve the config object and vm ",          &
                         esmf_log_info,rc=rc)
!
      call esmf_gridcompget(         gc_atm                             &  
                           ,config  =cf                                 &  
                           ,vm      =vm                                 &
                           ,rc      =rc)
!
      call atmos_err_msg(rc,'retrieve the config object and vm ',rcfinal)
! ---
!
!     call esmf_logwrite("get the vm",esmf_log_info,rc=rc)
!
!     call esmf_gridcompget(         gc_atm                          &
!                          ,vm      =vm                              &
!                          ,rc      =rc)
! 
!     call atmos_err_msg(rc,'get the vm',rcfinal)
!
!-----------------------------------------------------------------------
!***  set up parameters for mpi communications.
!***  use esmf utility to get pe identification and total number of pes
!***  (referred to here as nodes.)
!-----------------------------------------------------------------------
!
      call esmf_logwrite("get mype and nodes from vm",esmf_log_info,rc=rc)
!
      call esmf_vmget(vm                        &  !<-- the virtual machine
                     ,localpet=mype             &  !<-- local pe rank
                     ,petcount=num_pes          &  !<-- total # of tasks
                     ,rc      =rc)
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
      call atmos_err_msg(rc,'get mype and num_pes from vm',rcfinal)
!
!-----------------------------------------------------------------------
!***  establish the task layout including the quilt servers
!***  here in the main gridded component.  get the global
!***  mpi communicator and give it to setup_servers who will
!***  split it between forecast and quilt tasks.
!-----------------------------------------------------------------------
!
      call esmf_vmget(vm                                  &
                     ,mpicommunicator=mpi_intra           &  !<-- the global communicator
                     ,rc             =rc)
!
      call mpi_comm_dup(mpi_intra,mpi_intra_b,rc)
!
      call esmf_configgetattribute(cf                      &
                                  ,value =inpes            &  !<-- # of fcst tasks in i direction
                                  ,label ='inpes:'         &
                                  ,rc    =rc)
!
      call esmf_configgetattribute(cf                      &
                                  ,value =jnpes            &  !<-- # of fcst tasks in j direction
                                  ,label ='jnpes:'         &
                                  ,rc    =rc)
!
!ifdef NMM_B
!     call setup_servers(mype,inpes,jnpes,num_pes,mpi_intra_b)
!endif
!
      num_pes_fcst=num_pes
!
!***  note: at this point, num_pes is the number of forecast tasks only.
!
!-----------------------------------------------------------------------
!***  create the esmf grid.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  the first dimension of ncounts is the i dimension for parallelization.
!***  the second dimension of ncounts is the j dimension.
!-----------------------------------------------------------------------
!
      call esmf_configgetattribute(       cf                            &
                                  ,value =im                            &
                                  ,label ='im:'                         &
                                  ,rc    =rc)
!
      call esmf_configgetattribute(       cf                            &
                                  ,value =jm                            &
                                  ,label ='jm:'                         &
                                  ,rc    =rc)
!
!---------------------------------------------------------
!***  if this is a global mode forecast, extend im and jm.
!***  retrieve the mode from the config file.
!---------------------------------------------------------
!
      call esmf_configgetattribute(       cf                            &
                                  ,value =mode                          &
                                  ,label ='global:'                     &
                                  ,rc    =rc)
!
      if(trim(mode)=='.true.')then
        global=.true.
      else
        global=.false.
      endif
!
      if(global)then      !<-- global mode horizontal dimensions.
        ncounts(1)=im+2
        ncounts(2)=jm+2
      else                !<-- regional mode horizontal dimensions.
        ncounts(1)=im 
        ncounts(2)=jm
      endif
!
      max(1)=ncounts(1)
      max(2)=ncounts(2)
!
      min(1)=1
      min(2)=1
!
!---------------------------------------------------------
!***  condition to run only adiabatic (dynamics only)
!---------------------------------------------------------
!
      call esmf_configgetattribute(       cf                            &
                                  ,value =mode                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =rc)
!
      if(trim(mode)=='.true.')then
        physics_on=.false.
        print *,' initialize without physics coupling '
      else
        physics_on=.true.
        print *,' initialize with physics coupling '
      endif
!
!
!-----------------------------------------------------------------------
!***  now create the main gridded component's esmf grid.
!-----------------------------------------------------------------------
!
! Create the ESMF DistGrid_atmos.
!--------------------------------
      CALL ESMF_LogWrite("Create DistGrid_atmos", ESMF_LOG_INFO, rc = rc)

      DistGrid_atmos = ESMF_DistGridCreate(minIndex  = min,              &
                                           maxIndex  = max,              &
                                           regDecomp = (/inpes, jnpes/), &
                                           rc        = rc)

      call atmos_err_msg(rc, 'Create DistGrid_atmos', rcfinal)

! Create the ESMF grid_atmos based on the created ESMF DistGrid_atmos information.
!---------------------------------------------------------------------------------
      call esmf_logwrite("create grid_atmos",esmf_log_info,rc=rc)
!
      grid_atmos = ESMF_GridCreate(name     = "grid_atmos",   &
                                   distgrid = DistGrid_atmos, &
                                   rc       = rc)

      call atmos_err_msg(rc, 'create grid_atmos', rcfinal)
!
!-----------------------------------------------------------------------
!***  attach the esmf grid to the main gridded component.
!-----------------------------------------------------------------------
!
      call esmf_gridcompset(         gc_atm                          &
                           ,grid    =grid_atmos                      &
                           ,rc      =rc)
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
!***  using 'localcellcountperdim' from array i2 obtained in the
!***  previous call, generate all of the local task index limits.
!-----------------------------------------------------------------------
!
!ifdef NMM_B
!     i1 = i2(:, mype + 1)
!     if(mype<num_pes_fcst)then                                              !<-- this excludes quilt tasks.
!       call decomp(inpes,jnpes,im,jm,lm,global,i1,gc_atm)
!     endif
!endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  create the dynamics and physics import/export states
!***  with empty sizes and without considering any options
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!--------------
!***  dynamics --- create empty import/export states
!--------------
!
      imp_gfs_dyn=esmf_statecreate(statename="dynamics import"        &  
                                    ,statetype=esmf_state_import      &
                                    ,rc       =rc)
!
      exp_gfs_dyn=esmf_statecreate(statename="dynamics export"        &  
                                    ,statetype=esmf_state_export      &
                                    ,rc       =rc)
!
!-------------
!***  physics --- create empty import/export states
!-------------
!
      imp_gfs_phy=esmf_statecreate(statename="physics import"         &  
                                    ,statetype=esmf_state_import      &
                                    ,rc       =rc)
!
      exp_gfs_phy=esmf_statecreate(statename="physics export"         & 
                                    ,statetype=esmf_state_export      &
                                    ,rc       =rc)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  create the dynamics and physics gridded subcomponents.
!***  associate the config file with them.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!

! ====================================================================
!--------------
!***  dynamics --- grid comp create
!--------------
!
      call esmf_logwrite("create the dynamics gridded component"        &
                        ,esmf_log_info,rc=rc)
!
      gc_gfs_dyn=esmf_gridcompcreate( name      ="dynamics component"   &
                                     ,configfile='dyn_namelist.rc'      &
                                     ,rc        =rc)
!
      call atmos_err_msg(rc,'create the dynamics component',rcfinal)
!
!--------------
!***  dynamics --- register
!--------------
!
      call esmf_logwrite("register the dynamics gridded component"      &
                        ,esmf_log_info,rc=rc)
!
      call esmf_gridcompsetservices(gc_gfs_dyn                   & 
                                   ,gfs_dyn_setservices             & 
                                   ,rc)                                 
!
      call atmos_err_msg(rc,'register gfs dyn grid component',rcfinal)
!
!--------------
!***  dynamics --- run initialize step
!--------------
!
      call esmf_logwrite("initialize dynamics",esmf_log_info,rc=rc)
!
      call esmf_gridcompinitialize(gc_gfs_dyn                        &  
                                  ,importstate=imp_gfs_dyn           &  
                                  ,exportstate=exp_gfs_dyn           &  
                                  ,clock      =clock                 &  
                                  ,phase      =esmf_singlephase      &  
                                  ,rc         =rc)
!
      call atmos_err_msg(rc,'initialize dynamics',rcfinal)
!
!-----------------------------------------------------------------------
!***  both gridded subcomponents will simply be given esmf grids
!***  equivalent to that created for the main gridded component
!***  because the dynamics and physics are computed on the same grid. 
!-----------------------------------------------------------------------
! ---------------
!     dynamics  ---- attach grid as the same as atmos
! --------------
!
      call esmf_logwrite("attach the subcomponent grids"                &
                        ,esmf_log_info,rc=rc)
!
      grid_gfs_dyn=grid_atmos                                                   
      call esmf_gridcompset(         gc_gfs_dyn      &  
                           ,grid    =grid_gfs_dyn    &  
                           ,rc      =rc)
!
      call atmos_err_msg(rc,'attach dyn subcomponent grids',rcfinal)
!
! ======================================================================




      if( physics_on ) then
! ======================================================================
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  create the physics subcomponent, register, create states and initialize.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-------------
!***  physics --- grid comp create
!-------------
!
      call esmf_logwrite("create the physics gridded component"         &
                        ,esmf_log_info,rc=rc)
!
      gc_gfs_phy=esmf_gridcompcreate( name      ="physics component"    &
                                     ,configfile='phy_namelist.rc'      &
                                     ,rc        =rc)
!
      call atmos_err_msg(rc,'create the physics component',rcfinal)
!
!-------------
!***  physics --- register
!-------------
!
      call esmf_logwrite("register the physics gridded component"       &
                        ,esmf_log_info,rc=rc)
!
      call esmf_gridcompsetservices(gc_gfs_phy                   &
                                   ,gfs_phy_setservices             &
                                   ,rc)                                
!
      call atmos_err_msg(rc,'register gfs phy grid component',rcfinal)
!
!-------------
!***  physics --- run initialize step
!-------------
!
      call esmf_logwrite("initialize physics",esmf_log_info,rc=rc)
!
      call esmf_gridcompinitialize(            gc_gfs_phy            &  
                                  ,importstate=imp_gfs_phy           &  
                                  ,exportstate=exp_gfs_phy           &  
                                  ,clock      =clock                 &  
                                  ,phase      =esmf_singlephase      &  
                                  ,rc         =rc)
!
      call atmos_err_msg(rc,'initialize physics',rcfinal)
!
!-----------------------------------------------------------------------
!***  both gridded subcomponents will simply be given esmf grids
!***  equivalent to that created for the main gridded component
!***  because the dynamics and physics are computed on the same grid. 
!-----------------------------------------------------------------------
! ---------------
!     physics  ---- attach grid as the same as atmos
! --------------
!
      call esmf_logwrite("attach the subcomponent grids"             &
                        ,esmf_log_info,rc=rc)
!
      grid_gfs_phy=grid_atmos                                                  
      call esmf_gridcompset(         gc_gfs_phy      &  
                           ,grid    =grid_gfs_phy    &  
                           ,rc      =rc)
!
      call atmos_err_msg(rc,'process the subcomponent grids',rcfinal)
!
! =====================================================================
      endif 	! physics_on
!-----------------


! =====================================================================
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  create the coupler subcomponent, register, and initialize.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-------------
!***  coupler --- create grid comp.
!-------------
!
      call esmf_logwrite("create the two-way coupler component"         &
                        ,esmf_log_info,rc=rc)
!
      gc_atm_cpl=esmf_cplcompcreate(name="coupler component"            &
                                     ,rc  =rc)
!
      call atmos_err_msg(rc,'create the coupler component',rcfinal)
!
!-------------
!***  coupler --- register
!-------------
!
      call esmf_logwrite("register the coupler component"    &
                        ,esmf_log_info,rc=rc)
!
      call esmf_cplcompsetservices(gc_atm_cpl                &  
                                  ,atm_cpl_setservices              &  
                                  ,rc)                         
!
      call atmos_err_msg(rc,'register cpl grid component',rcfinal)
!
!-------------
!***  coupler --- run initialize step
!-------------
!
      call esmf_logwrite("initialize coupler",esmf_log_info,rc=rc)
!
      call esmf_cplcompinitialize(            gc_atm_cpl              &  
                                 ,importstate=exp_gfs_dyn             &  
                                 ,exportstate=imp_gfs_phy             &  
                                 ,clock      =clock                   &  
                                 ,rc         =rc)
!
      call atmos_err_msg(rc,'initialize coupler',rcfinal)
!
 

!-----------------------------------------------------------------------
!***  write the final error signal.
!-----------------------------------------------------------------------
!
      call atmos_err_msg_final(rcfinal,'main initialize step',rc)
!
      rc_init=rcfinal
!
!-----------------------------------------------------------------------
!
      end subroutine atmos_initialize
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      subroutine atmos_run(gc_atm,imp_state,exp_state              &
                              ,clock,rc_run)
!
!-----------------------------------------------------------------------
!***  run the main (top-level) gridded component.
!-----------------------------------------------------------------------
!
!ifdef NMM_B
!     use module_dm_parallel,only: iquilt_group,mpi_comm_inter_array
!     use module_quilt
!endif
      use module_digital_filter_gfs
!
!-----------------------------------------------------------------------
!
      type(esmf_gridcomp),intent(inout) :: gc_atm
      type(esmf_state),   intent(in)    :: imp_state
      type(esmf_state),   intent(inout) :: exp_state
      type(esmf_clock),   intent(inout) :: clock
      integer,            intent(out)   :: rc_run   
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),save :: mype,num_pes_fcst
      integer(kind=kint)      :: rc      ! error signal variable.
      integer(kind=kint)      :: rcfinal ! the final error signal variable.
      integer(kind=kint)      :: i,ier,ndfistep
      integer(kind=kint)      :: yy,mm,dd,hh,mns,sec
      integer(kind=kint),save :: dfihr=0

      type(esmf_config) 	:: cf           ! the config object
      type(esmf_time)		:: starttime
      type(esmf_time)		:: currtime
      type(esmf_time)		:: dfitime
      type(esmf_time)		:: halfdfitime
      type(esmf_timeinterval)	:: halfdfiintval
      type(esmf_timeinterval)	:: timestep

      character(50)             :: mode
      logical                   :: physics_on
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc     =esmf_success
      rcfinal=esmf_success
!
!-----------------------------------------------------------------------
!***  retrieve the vm to get the local pe identities.
!-----------------------------------------------------------------------
!
      call esmf_gridcompget(         gc_atm                             &  
                           ,vm      =vm                                 &
                           ,config  =cf                                 &  
                           ,rc      =rc)
!
      call esmf_vmget(vm                                                &
                     ,localpet=mype                                     & 
                     ,rc      =rc)
!
!-----------------------------------------------------------------------
!***  what is the total number of forecast tasks?
!-----------------------------------------------------------------------
!
      num_pes_fcst=inpes*jnpes
!
!ifdef NMM_B
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  send the quilt tasks on their way and let the forecast tasks
!***  proceed with the integration.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     quilt_fcst_split: if(mype>=num_pes_fcst)then 
!
!       call quilt(inpes,jnpes)
!
!     else
!
!endif
!                                                                                          
!-----------------------------------------------------------------------
!***  the tasks that reach this point will execute the forecast.
!-----------------------------------------------------------------------
!
!---------------------------------------------------------
!***  condition to run only adiabatic (dynamics only)
!---------------------------------------------------------
!
      call esmf_configgetattribute(       cf                            &
                                  ,value =mode                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =rc)
!
      if(trim(mode)=='.true.')then
        physics_on=.false.
        print *,' run without physics coupling '
      else
        physics_on=.true.
        print *,' run with physics coupling '
      endif
!
!
!-----------------------------------------------------------------------
!*** do digital filter initialize when it is necessary
!-----------------------------------------------------------------------
        call esmf_configgetattribute(       cf                          &
                                  ,value =dfihr                         &
                                  ,label ='nhours_dfini:'               &
                                  ,rc    =rc)
        if( dfihr.gt.0 ) then

! -------------------- initial stage -------------------------------
           print *,' start to set up dfini  '
           call esmf_timeintervalset(halfdfiintval			&
                                 ,h=dfihr,rc=rc)
           call esmf_clockget(clock					&
                          ,starttime=starttime				&
                          ,timestep =timestep				&
                          ,rc=rc)
           ndfistep = halfdfiintval / timestep
           halfdfitime = starttime + halfdfiintval
           dfitime = halfdfitime + halfdfiintval

           call digital_filter_dyn_init_gfs(imp_gfs_dyn,ndfistep)
! -------------------- inital summation  ------------------------------
           call digital_filter_dyn_sum_gfs(imp_gfs_dyn, mype)


           if( physics_on ) then
             call digital_filter_phy_init_gfs(imp_gfs_phy)
           endif	! physics_on

! -----------------------------------------------------------------------

        endif
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  the main integration time loop.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        timeloop: do while (.not.esmf_clockisstoptime(clock,rc))
!
!-----------------------------------------------------------------------
!
          call esmf_clockprint(      clock                              &
                              ,options="currtime string"                &
                              ,rc=rc)
!

!-----------------------------------------------------------------------
!***  execute the run step of the dynamics.
!***  this is the run subroutine specified in
!***  the dynamics register routine called in
!***  esmf_gridcompsetservices above.
!-----------------------------------------------------------------------
!
          call esmf_logwrite("execute dynamics",esmf_log_info,rc=rc)
!
          call esmf_gridcomprun(            gc_gfs_dyn          &  
                               ,importstate=imp_gfs_dyn         &  
                               ,exportstate=exp_gfs_dyn         &  
                               ,clock      =clock               &  
                               ,rc         =rc)        
!
          call atmos_err_msg(rc,'execute dynamics',rcfinal)
!
!
!-----------------------------------------------------------------------
!***  bring export data from the dynamics into the coupler
!***  and export it to the physics.
!-----------------------------------------------------------------------
!
          call esmf_logwrite("couple dyn_exp-to-phy_imp",      &
                             esmf_log_info,rc=rc)
!
          call esmf_cplcomprun(            gc_atm_cpl          &  
                              ,importstate=exp_gfs_dyn         &  
                              ,exportstate=imp_gfs_phy         &  
                              ,clock      =clock               &  
                              ,rc         =rc)        

          call atmos_err_msg(rc,'couple dyn-to-phy',rcfinal)



!-------------------------------
        if( physics_on ) then
!-----------------------------------------------------------------------
!***  execute the run step of the physics.
!-----------------------------------------------------------------------
!
          call esmf_logwrite("execute physics",esmf_log_info,rc=rc)
!
          call esmf_gridcomprun(            gc_gfs_phy            &  
                               ,importstate=imp_gfs_phy           &  
                               ,exportstate=exp_gfs_phy           &  
                               ,clock      =clock                 &  
                               ,rc         =rc)        
!
          call atmos_err_msg(rc,'execute physics',rcfinal)
!
!----------------
        else
!----------------
!
          call esmf_logwrite("pass phy_imp to phy_exp ",       &
                             esmf_log_info,rc=rc)
!
          call esmf_cplcomprun(            gc_atm_cpl          &  
                              ,importstate=imp_gfs_phy         &  
                              ,exportstate=exp_gfs_phy         &  
                              ,clock      =clock               &  
                              ,rc         =rc)        
!
          call atmos_err_msg(rc,'pass phy_imp-to-phy_exp',rcfinal)
!
!----------------
        endif	! physics_on
!----------------
!
!
!-----------------------------------------------------------------------
!***  bring export data from the physics into the coupler
!***  and export it to the dynamics.
!-----------------------------------------------------------------------
!
          call esmf_logwrite("couple phy_exp-to-dyn_imp",         &
                             esmf_log_info,rc=rc)

          call esmf_cplcomprun(            gc_atm_cpl             &  
                              ,importstate=exp_gfs_phy            &  
                              ,exportstate=imp_gfs_dyn            &  
                              ,clock      =clock                  & 
                              ,rc         =rc)        

          call atmos_err_msg(rc,'couple phy_exp-to-dyn_imp',rcfinal)
!
!-----------------------------------------------------------------------
!***  increment the timestep
!-----------------------------------------------------------------------
!
          call esmf_clockadvance(clock=clock                            &
                                ,rc   =rc)
!
!         call esmf_clockprint(      clock                              &
!                             ,options="currtime string"                &
!                             ,rc=rc)
!
          if(rc>0)then
            write(0,*)' problem with integration clock advance rc=',rc
          endif
!
!-----------------------------------------------------------------------
!*** do digital filter when it is necessary
!-----------------------------------------------------------------------
          call esmf_clockget(         clock				&
                            ,currtime=currtime				&
                            ,rc      =rc)

          if( dfihr.gt.0 ) then

! -------------------- summation stage ------------------------------
            call digital_filter_dyn_sum_gfs(imp_gfs_dyn, mype)

            if( physics_on ) then
              if( currtime .eq. halfdfitime ) then
                call digital_filter_phy_save_gfs(imp_gfs_phy)
              endif
            endif

! --------------------- final stage -----------------------------------
            if( currtime .eq. dfitime ) then
              print *,' dfi at finaldfitime '

              call digital_filter_dyn_average_gfs(imp_gfs_dyn, mype)
           
              if( physics_on )  then
                call digital_filter_phy_restore_gfs(imp_gfs_phy)
              endif

              call esmf_clockset(         clock				&
                                ,currtime=halfdfitime			&
                                ,rc      =rc)
              dfitime = starttime
              dfihr = 0
              print *,' dfi reset time to '
              call esmf_clockprint(      clock                              &
                              ,options="currtime string"                &
                              ,rc=rc)
            endif
! -----------------------------------------------------------------------

          endif
!
        enddo timeloop
!
!-----------------------------------------------------------------------
!***  execute the last step of the dynamics.
!***  this is the run subroutine specified in
!***  the dynamics register routine called in
!***  esmf_gridcompsetservices above.
!-----------------------------------------------------------------------
!
          call esmf_logwrite("last step dynamics",esmf_log_info,rc=rc)
!
          call esmf_gridcomprun(            gc_gfs_dyn          &  
                               ,importstate=imp_gfs_dyn         &  
                               ,exportstate=exp_gfs_dyn         &  
                               ,clock      =clock               &  
                               ,rc         =rc)        
!
          call atmos_err_msg(rc,'last step dynamics',rcfinal)
!
!ifdef NMM_B
!-----------------------------------------------------------------------
!***  the forecast is finished therefore shut down the quilt servers.
!***  forecast task 0 sends -999 to task 0 in each quilt group.
!-----------------------------------------------------------------------
!
!       if(mype==0)then
!         do i=1,iquilt_group
!           call mpi_send(-999,1,mpi_integer,0,0                        &
!                        ,mpi_comm_inter_array(i),ier)
!         enddo
!       endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     endif quilt_fcst_split
!endif
!
!-----------------------------------------------------------------------
!
          call esmf_clockprint(      clock                              &
                              ,options="stoptime string"                &
                              ,rc=rc)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  the final error signal information.
!-----------------------------------------------------------------------
!
      call atmos_err_msg_final(rcfinal,'main run step',rc)
!
      rc_run=rcfinal
!
!-----------------------------------------------------------------------
!
      end subroutine atmos_run
!
!-----------------------------------------------------------------------
!***********************************************************************
!***********************************************************************
!-----------------------------------------------------------------------
!
      subroutine atmos_finalize(gc_atm,imp_state,exp_state         &
                              ,clock,rc_finalize)
!
!-----------------------------------------------------------------------
!***  this routine finalizes the main gridded component.
!-----------------------------------------------------------------------
!
      type(esmf_gridcomp),intent(inout) :: gc_atm
      type(esmf_state),   intent(inout) :: imp_state
      type(esmf_state),   intent(inout) :: exp_state
      type(esmf_clock),   intent(in)    :: clock
      type(esmf_config)                 :: cf           ! the config object
      integer,            intent(out)   :: rc_finalize
      real(kind=8)	timef
      real(kind=16)	stime, etime
      character(50)             :: mode
      logical                   :: physics_on
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer :: rc
      integer :: rcfinal      ! the final error signal variable.
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc     =esmf_success
      rcfinal=esmf_success
!
!-----------------------------------------------------------------------
!***  finalize each of the subcomponents.
!-----------------------------------------------------------------------
      call esmf_gridcompget(         gc_atm                             &  
                           ,vm      =vm                                 &
                           ,config  =cf                                 &  
                           ,rc      =rc)
!
!---------------------------------------------------------
!***  condition to run only adiabatic (dynamics only)
!---------------------------------------------------------
!
      call esmf_configgetattribute(       cf                            &
                                  ,value =mode                          &
                                  ,label ='adiabatic:'                  &
                                  ,rc    =rc)
!
      if(trim(mode)=='.true.')then
        physics_on=.false.
        print *,' finalize without physics coupling. '
      else
        physics_on=.true.
        print *,' finalize with physics coupling. '
      endif

!-----------------------------------------------------------------------
!-------------
!***  coupler
!-------------
!  
! finalize coupler component
!
      call esmf_logwrite("finalize coupler subcomponents",		&
                          esmf_log_info,rc=rc)
!
      call esmf_cplcompfinalize(            gc_atm_cpl                	&
                               ,importstate=exp_gfs_dyn           	&
                               ,exportstate=imp_gfs_phy           	&
                               ,rc         =rc)
!
      call atmos_err_msg(rc,'finalize coupler subcomponents',rcfinal)
!
! destry coupler component
!
      call esmf_logwrite("destroy coupler subcomponents",		&
                         esmf_log_info,rc=rc)
!
      call esmf_cplcompdestroy(gc_atm_cpl,rc=rc)
!
      call atmos_err_msg(rc,'destroy coupler subcomponents',rcfinal)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!--------------
!***  dynamics  
!--------------
!  
! finalize dynamics component
!
      call esmf_logwrite("finalize dynamics subcomponents",		&
                         esmf_log_info,rc=rc)
!
      call esmf_gridcompfinalize(gc_gfs_dyn  				&
                                ,importstate=imp_gfs_dyn          	&
                                ,exportstate=exp_gfs_dyn          	&
                                ,rc         =rc)
!
      call atmos_err_msg(rc,'finalize dynamics subcomponents',rcfinal)
!
! destroy dynamics states
!
      stime=timef()

      call esmf_logwrite("destroy dynamics states",esmf_log_info,rc=rc)
!
      call esmf_statedestroy(imp_gfs_dyn, rc=rc)
      call esmf_statedestroy(exp_gfs_dyn, rc=rc)
!
      call atmos_err_msg(rc,'destroy dynamics states',rcfinal)

      etime=timef()
      print *,' time spend for destroy state is ',etime-stime
!
! destroy dynamics component
!
      stime=timef()

      call esmf_logwrite("destroy dynamics subcomponents",		&
                         esmf_log_info,rc=rc)
!
! hmhj cannot destroy gc_gfs_dyn, why?
      call esmf_gridcompdestroy(gc_gfs_dyn, rc=rc)
!
      call atmos_err_msg(rc,'destroy dynamics subcomponents',rcfinal)

      etime=timef()
      print *,' time spend for destroy dynamics comp is ',etime-stime
!
!-----------------------------------------------------------------------


!--------------------------
      if( physics_on ) then
!-----------------------------------------------------------------------
!-------------
!***  physics
!-------------
!  
! finalize physics component
!
      call esmf_logwrite("finalize physics subcomponents",		&
                         esmf_log_info,rc=rc)
!
      call esmf_gridcompfinalize(gc_gfs_phy           	                &
                                ,importstate=imp_gfs_phy          	&
                                ,exportstate=exp_gfs_phy          	&
                                ,rc         =rc)
!
      call atmos_err_msg(rc,'finalize dynamics subcomponents',rcfinal)
!
! destry physics states
!
      call esmf_logwrite("destroy physics states",			&
                          esmf_log_info,rc=rc)
!
      call esmf_statedestroy(imp_gfs_phy, rc=rc)
      call esmf_statedestroy(exp_gfs_phy, rc=rc)
!
      call atmos_err_msg(rc,'destroy physics states',rcfinal)
!
! destroy physics component
!
      call esmf_logwrite("destroy physics subcomponents",		&
                         esmf_log_info,rc=rc)
!
      call esmf_gridcompdestroy(gc_gfs_phy, rc=rc)
!
      call atmos_err_msg(rc,'destroy dynamics subcomponents',rcfinal)
!
!-----------------------------------------------------------------------
      endif 	! physics_on
!---------------

!-----------------------------------------------------------------------
!***  destroy the esmf clock.
!-----------------------------------------------------------------------
!
      stime=timef()

      call esmf_logwrite("destroy the esmf clock",esmf_log_info,rc=rc)
!
      call esmf_clockdestroy(clock,rc=rc)
!
      call atmos_err_msg(rc,'destroy the esmf clock',rcfinal)

      etime=timef()
      print *,' time spend for destroy clock is ',etime-stime
!
!-----------------------------------------------------------------------
!***  the final error signal information.
!-----------------------------------------------------------------------
!
      call atmos_err_msg_final(rcfinal,'main finalize step',rc)
!
      rc_finalize=rcfinal
!
!-----------------------------------------------------------------------
!
      end subroutine atmos_finalize
!
!-----------------------------------------------------------------------
!
      end module atmos_grid_comp_mod
!
!-----------------------------------------------------------------------
