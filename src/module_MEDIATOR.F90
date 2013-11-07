#include "./ESMFVersionDefine.h"
#ifdef WITH_NUOPC

module module_MEDIATOR

  !-----------------------------------------------------------------------------
  ! NEMS Mediator Component.
  !
  ! The Mediator has two Run() phases:
  !
  !   * Run(phase=1) covers the more frequent interaction with the ATM
  !     component. The ATM exports some fields as time averages over the 
  !     integration period. Other fields are exported as instantaneous. For 
  !     the time averaged fields the Mediator is responsible to continue
  !     the averaging.
  !
  !   * Run(phase=2) is invoked for the less frequent interaction with the OCN
  !     component. Here the averaged and instantanous ATM fields are passed on
  !     to the OCN component, and OCN export fields are received by the
  !     Mediator and forwarded to the ATM component.
  !
  ! The two phases are operating on different time scales, and hence require
  ! two separate internal Component Clocks. The NUOPC layer accesses a
  ! Component's Clock through the ESMF CompGet() interface, regardless of the
  ! phase. Phase specific Clocks are implemented by swapping Clocks during the
  ! phase specific "label_SetRunClock" specialization method. These Clock
  ! objects are stored in the Component instance's own internal state.
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Mediator, only: &
    mediator_routine_SS             => routine_SetServices, &
    mediator_routine_Run            => routine_Run, &
    mediator_type_IS                => type_InternalState, &
    mediator_label_IS               => label_InternalState, &
    mediator_label_DataInitialize   => label_DataInitialize, &
    mediator_label_Advance          => label_Advance, &
    mediator_label_CheckImport      => label_CheckImport, &
    mediator_label_TimestampExport  => label_TimestampExport, &
    mediator_label_SetRunClock      => label_SetRunClock
  
  implicit none
  
  private
  
  ! private internal state to keep instance data
  type InternalStateStruct
    type(ESMF_Clock)      :: fastClock
    type(ESMF_Clock)      :: slowClock
    integer   :: slice  ! slice counter for writing to NetCDF file
    type(ESMF_Field)      :: rsns_accumulator
  end type

  type InternalState
    type(InternalStateStruct), pointer :: wrap
  end type

  public SetServices
  
  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! the NUOPC mediator component will register the generic methods
    call mediator_routine_SS(gcomp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! Provide InitP0 to overwrite the default IPD00 with IPD02
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      InitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! IPDv02 requires InitP1, where Fields should be advertised
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP1, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! IPDv02 requires InitP2, where Fields should be realized,
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_INITIALIZE, &
      userRoutine=InitializeP2, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! attach specializing method(s)
    call ESMF_MethodAdd(gcomp, label=mediator_label_DataInitialize, &
      userRoutine=DataInitialize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! Run phase 2 entry point
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_RUN, &
      userRoutine=mediator_routine_Run, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! overwrite Finalize
    call ESMF_GridCompSetEntryPoint(gcomp, ESMF_METHOD_FINALIZE, &
      userRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! attach specializing methods for Run(phase=1) "fast"
    call ESMF_MethodAdd(gcomp, label=mediator_label_SetRunClock, &
      index=1, userRoutine=SetRunClock_fast, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_MethodAdd(gcomp, label=mediator_label_CheckImport, &
      index=1, userRoutine=CheckImport_fast, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_MethodAdd(gcomp, label=mediator_label_TimestampExport, &
      index=1, userRoutine=TimestampExport_fast, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_MethodAdd(gcomp, label=mediator_label_Advance, &
      index=1, userRoutine=Advance_fast, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing methods for Run(phase=2) "slow"
    call ESMF_MethodAdd(gcomp, label=mediator_label_SetRunClock, &
      index=2, userRoutine=SetRunClock_slow, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_MethodAdd(gcomp, label=mediator_label_Advance, &
      index=2, userRoutine=Advance_slow, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! extend the NUOPC Field Dictionary to hold required entries
    if (.not.NUOPC_FieldDictionaryHasEntry("air_temperature_at_lowest_level")) &
      then
      call NUOPC_FieldDictionaryAddEntry( &
        standardName="air_temperature_at_lowest_level", &
        canonicalUnits="K", &
        defaultLongName="Air Temperature at Lowest Level", &
        defaultShortName="atll", rc=rc);
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
      
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP0(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)   :: gcomp
    type(ESMF_State)      :: importState, exportState
    type(ESMF_Clock)      :: clock
    integer, intent(out)  :: rc
    
    ! local variables    
    character(len=NUOPC_PhaseMapStringLength) :: initPhases(4)
    
    rc = ESMF_SUCCESS

    initPhases(1) = "IPDv02p1=1"
    initPhases(2) = "IPDv02p3=2"
    initPhases(3) = "IPDv02p4=3"
    initPhases(4) = "IPDv02p5=5"
    
    call ESMF_AttributeSet(gcomp, &
      name="InitializePhaseMap", valueList=initPhases, &
      convention="NUOPC", purpose="General", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------

  subroutine InitializeP1(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    rc = ESMF_SUCCESS
    
    ! importable fields:
    call NUOPC_StateAdvertiseFields(importState, StandardNames=(/ &
      "sea_surface_temperature", &
      "surface_net_downward_shortwave_flux", &
      "air_temperature_at_lowest_level" &
       /), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! exportable fields:
    call NUOPC_StateAdvertiseFields(exportState, StandardNames=(/ &
      "sea_surface_temperature", &
      "surface_net_downward_shortwave_flux", &
      "air_temperature_at_lowest_level" &
       /), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine InitializeP2(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    ! local variables    
    type(ESMF_Field)            :: field_sst, field
    type(ESMF_Grid)             :: gridIn, gridOut
    integer                     :: i, j
    real(kind=ESMF_KIND_R8),pointer :: lonPtr(:,:), latPtr(:,:)
    type(InternalState)         :: is
    integer                     :: stat
    type(ESMF_Config)           :: config
    real(ESMF_KIND_R8)          :: intervalSec
    type(ESMF_TimeInterval)     :: timeStep
    
    rc = ESMF_SUCCESS
    
    ! Allocate memory for the internal state and set it in the Component.
    allocate(is%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of the internal state memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridCompSetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! Initialize the internal state members
    is%wrap%slice = 1

    ! create a DUMMY Grid object for import and export Fields
    gridIn = ESMF_GridCreate1PeriDim(minIndex=(/1,1/), maxIndex=(/500,200/),&
      indexflag=ESMF_INDEX_GLOBAL, coordSys=ESMF_COORDSYS_SPH_DEG, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridAddCoord(gridIn, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(gridIn, coordDim=1, farrayPtr=lonPtr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    call ESMF_GridGetCoord(gridIn, coordDim=2, farrayPtr=latPtr, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    do j=lbound(lonPtr,2),ubound(lonPtr,2)
    do i=lbound(latPtr,1),ubound(latPtr,1)
      lonPtr(i,j) = 360./real(500) * (i-1)
      latPtr(i,j) = 100./real(200) * (j-1) - 50.
    enddo
    enddo
      
    gridOut = gridIn ! for now out same as in

    ! conditionally realize or remove Fields
    
    ! importable/exportable field: sea_surface_temperature
    if (NUOPC_StateIsFieldConnected(importState, fieldName="sst") .or. &
      NUOPC_StateIsFieldConnected(exportState, fieldName="sst")) then
      ! create a common Field
      field_sst = ESMF_FieldCreate(name="sst", grid=gridIn, &
        typekind=ESMF_TYPEKIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    ! importable field: sea_surface_temperature
    if (NUOPC_StateIsFieldConnected(importState, fieldName="sst")) then
      ! realize a connected Field
      call NUOPC_StateRealizeField(importState, field=field_sst, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
      ! remove a not connected Field from State
      call ESMF_StateRemove(importState, (/"sst"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    ! exportable field: sea_surface_temperature
    if (NUOPC_StateIsFieldConnected(exportState, fieldName="sst")) then
      ! realize a connected Field
      call NUOPC_StateRealizeField(exportState, field=field_sst, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
      ! remove a not connected Field from State
      call ESMF_StateRemove(exportState, (/"sst"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    ! importable field: air_temperature_at_lowest_level
    if (NUOPC_StateIsFieldConnected(importState, fieldName="atll")) then
      ! realize a connected Field
      field = ESMF_FieldCreate(name="atll", grid=gridIn, &
        typekind=ESMF_TYPEKIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_StateRealizeField(importState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
      ! remove a not connected Field from State
      call ESMF_StateRemove(importState, (/"atll"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    ! exportable field: air_temperature_at_lowest_level
    if (NUOPC_StateIsFieldConnected(exportState, fieldName="atll")) then
      ! realize a connected Field
      field = ESMF_FieldCreate(name="atll", grid=gridOut, &
        typekind=ESMF_TYPEKIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_StateRealizeField(exportState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
      ! remove a not connected Field from State
      call ESMF_StateRemove(exportState, (/"atll"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    ! importable field: surface_net_downward_shortwave_flux
    if (NUOPC_StateIsFieldConnected(importState, fieldName="rsns")) then
      ! realize a connected Field
      field = ESMF_FieldCreate(name="rsns", grid=gridIn, &
        typekind=ESMF_TYPEKIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_StateRealizeField(importState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
      ! remove a not connected Field from State
      call ESMF_StateRemove(importState, (/"rsns"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    ! exportable field: surface_net_downward_shortwave_flux
    if (NUOPC_StateIsFieldConnected(exportState, fieldName="rsns")) then
      ! realize a connected Field
      field = ESMF_FieldCreate(name="rsns", grid=gridOut, &
        typekind=ESMF_TYPEKIND_R8, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call NUOPC_StateRealizeField(exportState, field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    else
      ! remove a not connected Field from State
      call ESMF_StateRemove(exportState, (/"rsns"/), rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif

    ! need an internal accumulator Field -> could be made dependent imp/exp 
    is%wrap%rsns_accumulator = ESMF_FieldCreate(name="rsns", &
      grid=gridOut, typekind=ESMF_TYPEKIND_R8, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Initialize the internal clocks
    
    ! both fast and slow clocks start out as copies of the incoming clock
    is%wrap%fastClock = ESMF_ClockCreate(clock, rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    is%wrap%slowClock = ESMF_ClockCreate(clock, rc=rc)
    ESMF_ERR_RETURN(rc,rc)

    config = ESMF_ConfigCreate(rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    call ESMF_ConfigLoadFile(config, "nems.configure", rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
    ! ATM coupling interval -> fast time step
    call ESMF_ConfigGetAttribute(config, intervalSec, &
      label="med_atm_coupling_interval_sec:", default=-1.0_ESMF_KIND_R8, &
      rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
    if (intervalSec>0._ESMF_KIND_R8) then
      ! The coupling time step was provided
      call ESMF_TimeIntervalSet(timeStep, s_r8=intervalSec, rc=rc)
      ESMF_ERR_RETURN(rc,rc)
      call ESMF_ClockSet(is%wrap%fastClock, timestep=timeStep, rc=rc)
      ESMF_ERR_RETURN(rc,rc)
    endif
    
    ! OCN coupling interval -> slow time step
    call ESMF_ConfigGetAttribute(config, intervalSec, &
      label="med_ocn_coupling_interval_sec:", default=-1.0_ESMF_KIND_R8, &
      rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
    if (intervalSec>0._ESMF_KIND_R8) then
      ! The coupling time step was provided
      call ESMF_TimeIntervalSet(timeStep, s_r8=intervalSec, rc=rc)
      ESMF_ERR_RETURN(rc,rc)
      call ESMF_ClockSet(is%wrap%slowClock, timestep=timeStep, rc=rc)
      ESMF_ERR_RETURN(rc,rc)
    endif
    
    call ESMF_ConfigDestroy(config, rc=rc)
    ESMF_ERR_RETURN(rc,rc)
    
  end subroutine
  
  !-----------------------------------------------------------------------------

  subroutine DataInitialize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(ESMF_Time)             :: time
    type(ESMF_Field)            :: field
    type(ESMF_StateItem_Flag)   :: itemType
    logical                     :: neededCurrent
    logical                     :: allDone
    type(InternalState)         :: is
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:)

    rc = ESMF_SUCCESS
    
    ! the MED needs valid ATM export Fields to initialize its internal state

    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! Get the internal state from Component.
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! get the current time out of the clock
    call ESMF_ClockGet(clock, currTime=time, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    allDone = .true.  ! flag that can be reset if anything is not found done
    
    ! check that required Fields in the importState show correct timestamp
    
    ! -> check for "atll"
    call ESMF_StateGet(importState, itemName="atll", itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      call ESMF_StateGet(importState, field=field, itemName="atll", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      neededCurrent = NUOPC_FieldIsAtTime(field, time, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      
      if (.not.neededCurrent) then
        call ESMF_LogWrite("MED - Initialize-Data-Dependency NOT YET SATISFIED!!!", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        allDone = .false.
      else
        call ESMF_LogWrite("MED - Initialize-Data-Dependency SATISFIED!!!", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      endif
    endif
      
    ! -> check for "rsns"
    call ESMF_StateGet(importState, itemName="rsns", itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      call ESMF_StateGet(importState, field=field, itemName="rsns", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      neededCurrent = NUOPC_FieldIsAtTime(field, time, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      
      if (.not.neededCurrent) then
        call ESMF_LogWrite("MED - Initialize-Data-Dependency NOT YET SATISFIED!!!", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        allDone = .false.
      else
        call ESMF_LogWrite("MED - Initialize-Data-Dependency SATISFIED!!!", &
          ESMF_LOGMSG_INFO, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
          
        ! initialize the "rsns" accumulator
        call ESMF_FieldGet(is%wrap%rsns_accumulator, farrayPtr=dataPtr, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
          
        ! For the real case this should probably use the "rsns" field from the
        ! importState and do something with it as a sensible starting point
        ! for the accumulation field so that the OCN receives a sensible
        ! "rsns" field during its first time step. However, here for testing
        ! I simply initialize to zero.
        dataPtr(:,:) = 0._ESMF_KIND_R8  
          
      endif
    endif

    if (allDone) then
      ! -> set InitializeDataComplete Component Attribute to "true", indicating
      ! to the driver that this Component has fully initialized its data
      call ESMF_AttributeSet(gcomp, &
        name="InitializeDataComplete", value="true", &
        convention="NUOPC", purpose="General", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    endif
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine SetRunClock_fast(gcomp, rc)
    type(ESMF_GridComp)   :: gcomp
    integer, intent(out)  :: rc
    
    ! local variables
    type(InternalState)     :: is_local
    type(mediator_type_IS)  :: is

    rc = ESMF_SUCCESS
    
    ! query component for its internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set fastClock to be the component clock
    call ESMF_GridCompSet(gcomp, clock=is_local%wrap%fastClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! query component for its internal state
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, mediator_label_IS, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! check and set the model clock against the driver clock
    call NUOPC_GridCompCheckSetClock(gcomp, is%wrap%driverClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="NUOPC INCOMPATIBILITY DETECTED: between model and driver clocks", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine SetRunClock_slow(gcomp, rc)
    type(ESMF_GridComp)   :: gcomp
    integer, intent(out)  :: rc
    
    ! local variables
    type(InternalState)     :: is_local
    type(mediator_type_IS)  :: is

    rc = ESMF_SUCCESS
    
    ! query component for its internal state
    nullify(is_local%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is_local, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! set fastClock to be the component clock
    call ESMF_GridCompSet(gcomp, clock=is_local%wrap%slowClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! query component for its internal state
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, mediator_label_IS, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! check and set the model clock against the driver clock
    call NUOPC_GridCompCheckSetClock(gcomp, is%wrap%driverClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, &
      msg="NUOPC INCOMPATIBILITY DETECTED: between model and driver clocks", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine CheckImport_fast(gcomp, rc)
    type(ESMF_GridComp)   :: gcomp
    integer, intent(out)  :: rc
    
    ! This is the routine that ensures that the import Fields come in with
    ! the correct time stamps during the "fast" cycle: 
    ! -> Fields from the ATM must be at stopTime.
    ! -> Fields from the OCN must not be checked.
    
    ! local variables
    type(ESMF_Clock)        :: clock
    type(ESMF_Time)         :: stopTime
    type(ESMF_State)        :: importState
    type(ESMF_Field)        :: field
    logical                 :: atCorrectTime

    rc = ESMF_SUCCESS
    
    ! query the Component for its Clock and importState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! get the current time out of the Clock
    call ESMF_ClockGet(clock, stopTime=stopTime, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! check fields from ATM to be at stopTime
    if (NUOPC_StateIsFieldConnected(importState, fieldName="atll")) then
      call ESMF_StateGet(importState, itemName="atll", field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      atCorrectTime = NUOPC_FieldIsAtTime(field, stopTime, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (.not.atCorrectTime) then
        !TODO: introduce and use INCOMPATIBILITY return codes!!!!
        call ESMF_LogSetError(ESMF_RC_ARG_BAD, &
          msg="NUOPC INCOMPATIBILITY DETECTED: Import Fields not at correct time", &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)
        return  ! bail out
      endif
    endif

  end subroutine

  !-----------------------------------------------------------------------------
  
  subroutine TimestampExport_fast(gcomp, rc)
    type(ESMF_GridComp)   :: gcomp
    integer, intent(out)  :: rc
    
    ! This is the routine that applies the time stamp on the export Fields
    ! during the "fast" cycle: 
    ! -> Fields receive the currTime stamp of the clock.

    ! local variables
    type(ESMF_Clock)      :: clock
    type(ESMF_State)      :: exportState

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(gcomp, clock=clock, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! update timestamp on export Fields
    call NUOPC_StateSetTimestamp(exportState, clock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Advance_fast(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is
    type(ESMF_Field)            :: field
    type(ESMF_StateItem_Flag)   :: itemType
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:), accDataPtr(:,:)
    integer                     :: i,j
    
    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! Get the internal state from Component.
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MEDIATOR does the mediation of Fields that come in on the
    ! importState with a timestamp consistent to the stopTime of the 
    ! mediators Clock.
    
    ! The Mediator uses the data on the import Fields to update the data
    ! held by Fields in the exportState.
    
    ! After this routine returns the generic Mediator will correctly
    ! timestamp the export Fields to the stopTime, and then update 
    ! the Mediator Clock to:
    !
    !       currTime -> currTime + timeStep
    !
    ! Where the timeStep is equal to the parent timeStep.
    
    call NUOPC_ClockPrintCurrTime(clock, &
      "-------->MED Advance() mediating for: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! accumulation of the "rsns" Field
    call ESMF_StateGet(importState, itemName="rsns", itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      ! "rsns" is available on the import side -> do accumulation
      call ESMF_StateGet(importState, itemName="rsns", field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(is%wrap%rsns_accumulator, farrayPtr=accDataPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      do j=lbound(dataPtr,2),ubound(dataPtr,2)
      do i=lbound(dataPtr,1),ubound(dataPtr,1)
        accDataPtr(i,j) = accDataPtr(i,j) + dataPtr(i,j)
      enddo
      enddo
    endif
         
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Advance_slow(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(ESMF_Clock)            :: clock
    type(ESMF_State)            :: importState, exportState
    type(InternalState)         :: is
    type(ESMF_Field)            :: field
    type(ESMF_StateItem_Flag)   :: itemType
    real(ESMF_KIND_R8), pointer :: dataPtr(:,:), accDataPtr(:,:)
    integer                     :: i,j

    rc = ESMF_SUCCESS
    
    ! query the Component for its clock, importState and exportState
    call ESMF_GridCompGet(gcomp, clock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! Get the internal state from Component.
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! HERE THE MEDIATOR does the mediation of Fields that come in on the
    ! importState with a timestamp consistent to the currTime of the 
    ! mediators Clock.
    
    ! The Mediator uses the data on the import Fields to update the data
    ! held by Fields in the exportState.
    
    ! After this routine returns the generic Mediator will correctly
    ! timestamp the export Fields to the currentTime, and then update 
    ! the Mediator Clock to:
    !
    !       currTime -> currTime + timeStep
    !
    ! Where the timeStep is equal to the parent timeStep.
    
    call NUOPC_ClockPrintCurrTime(clock, &
      "-------->MED Advance() mediating for: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    
    ! output the imported "sst" fields
    call ESMF_StateGet(importState, itemName="sst", itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      ! Write the imported SST Field into a NetCDF file as timeslice.
      call ESMF_StateGet(importState, itemName="sst", field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
#if 1
      if (ESMF_IO_PIO_PRESENT .and. ESMF_IO_NETCDF_PRESENT) then
        call ESMF_FieldWrite(field, file="field_sst.nc", &
          timeslice=is%wrap%slice, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          call ESMF_Finalize(endflag=ESMF_END_ABORT)
      endif
#endif
      is%wrap%slice = is%wrap%slice + 1
    endif
    
    ! copy the "rsns" accumulator to the export and zero out the accumulator
    call ESMF_StateGet(exportState, itemName="rsns", itemType=itemType, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
    if (itemType /= ESMF_STATEITEM_NOTFOUND) then
      call ESMF_StateGet(exportState, itemName="rsns", field=field, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        call ESMF_Finalize(endflag=ESMF_END_ABORT)
      call ESMF_FieldGet(field, farrayPtr=dataPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call ESMF_FieldGet(is%wrap%rsns_accumulator, farrayPtr=accDataPtr, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      do j=lbound(dataPtr,2),ubound(dataPtr,2)
      do i=lbound(dataPtr,1),ubound(dataPtr,1)
        dataPtr(i,j) = accDataPtr(i,j)
        accDataPtr(i,j) = 0._ESMF_KIND_R8
      enddo
      enddo
    endif
    
  end subroutine

  !-----------------------------------------------------------------------------

  subroutine Finalize(gcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gcomp
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc
    
    type(InternalState)  :: is
    integer              :: stat

    rc = ESMF_SUCCESS
  
    ! Get the internal state from Component.
    nullify(is%wrap)
    call ESMF_GridCompGetInternalState(gcomp, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
    ! Destroy objects inside of internal state.
    call ESMF_FieldDestroy(is%wrap%rsns_accumulator, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! Deallocate the internal state memory.
    deallocate(is%wrap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of internal state memory failed.", &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out
      
  end subroutine

  !-----------------------------------------------------------------------------

end module
#endif
