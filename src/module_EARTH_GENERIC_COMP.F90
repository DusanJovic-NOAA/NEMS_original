#include "./ESMFVersionDefine.h"
#ifdef WITH_NUOPC

module module_EARTH_GENERIC_COMP

  !-----------------------------------------------------------------------------
  ! Generic NEMS Earth Driver Component
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Driver, &
    Driver_routine_SS             => routine_SetServices, &
    Driver_type_IS                => type_InternalState, &
    Driver_type_ISS               => type_InternalStateStruct, &
    Driver_label_IS               => label_InternalState, &
    Driver_label_SetModelCount    => label_SetModelCount, &
    Driver_label_SetModelPetLists => label_SetModelPetLists, &
    Driver_label_SetModelServices => label_SetModelServices, &
    Driver_label_Finalize         => label_Finalize

  implicit none
  
  private
  
  public routine_SetServices
  public type_InternalState, type_InternalStateStruct
  public label_InternalState, label_SetModelPetLists
  public label_SetModelServices, label_Finalize
  
  public NUOPC_DriverAddComp, NUOPC_DriverGetComp
  
  character(*), parameter :: &
    label_InternalState = "NemsEarthGeneric_InternalState"
  character(*), parameter :: &
    label_SetModelPetLists = "NemsEarthGeneric_SetModelPetLists"
  character(*), parameter :: &
    label_SetModelServices = "NemsEarthGeneric_SetModelServices"
  character(*), parameter :: &
    label_Finalize = "NemsEarthGeneric_Finalize"
  
  type type_InternalStateStruct
    integer, pointer    :: atmPetList(:)
    integer, pointer    :: ocnPetList(:)
    integer, pointer    :: icePetList(:)
    integer, pointer    :: medPetList(:)
    real(ESMF_KIND_R8)  :: medAtmCouplingIntervalSec
    real(ESMF_KIND_R8)  :: medOcnCouplingIntervalSec
  end type

  type type_InternalState
    type(type_InternalStateStruct), pointer :: wrap
  end type
  
  integer, parameter :: medPhase_slow = 1 ! must match MED implementation
  integer, parameter :: medPhase_fast_before = 2 ! must match MED implementation
  integer, parameter :: medPhase_fast_after  = 3 ! must match MED implementation

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine routine_SetServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    character(ESMF_MAXSTR):: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! derive from generic NUOPC_Driver
    call NUOPC_CompDerive(driver, Driver_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! attach specializing method(s)
    call NUOPC_CompSpecialize(driver, specLabel=Driver_label_SetModelCount, &
      specRoutine=SetModelCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call NUOPC_CompSpecialize(driver, specLabel=Driver_label_SetModelPetLists, &
      specRoutine=SetModelPetLists, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call NUOPC_CompSpecialize(driver, specLabel=Driver_label_SetModelServices, &
      specRoutine=SetModelServices, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call NUOPC_CompSpecialize(driver, specLabel=Driver_label_Finalize, &
      specRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------------
  
  subroutine SetModelCount(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    character(ESMF_MAXSTR):: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! set the modelCount for ATM-OCN-ICE-MED coupling
    call NUOPC_DriverSet(driver, modelCount=4, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------------
  
  subroutine SetModelPetLists(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    integer                   :: localrc, stat, i
    type(type_InternalState)  :: is
    type(Driver_type_IS)      :: superIS
    character(ESMF_MAXSTR)    :: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! allocate memory for this internal state and set it in the Component
    allocate(is%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of internal state memory failed.", &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out
    call ESMF_UserCompSetInternalState(driver, label_InternalState, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! nullify the petLists
    nullify(is%wrap%atmPetList)
    nullify(is%wrap%ocnPetList)
    nullify(is%wrap%icePetList)
    nullify(is%wrap%medPetList)
    
    ! SPECIALIZE by calling into optional attached method to set modelPetLists
    call ESMF_MethodExecute(driver, label=label_SetModelPetLists, &
      userRc=localrc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out

    ! query Component for super internal State
    nullify(superIS%wrap)
    call ESMF_UserCompGetInternalState(driver, Driver_label_IS, superIS, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! set the petLists
    i = 1 ! initialize model counter
    if (associated(is%wrap%atmPetList)) then
      superIS%wrap%modelPetLists(i)%petList => is%wrap%atmPetList
      i = i+1
    endif
    if (associated(is%wrap%ocnPetList)) then
      superIS%wrap%modelPetLists(i)%petList => is%wrap%ocnPetList
      i = i+1
    endif
    if (associated(is%wrap%icePetList)) then
      superIS%wrap%modelPetLists(i)%petList => is%wrap%icePetList
      i = i+1
    endif
    if (associated(is%wrap%medPetList)) then
      superIS%wrap%modelPetLists(i)%petList => is%wrap%medPetList
    endif
    
  end subroutine
  
  !-----------------------------------------------------------------------------
  
  subroutine SetModelServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    integer                   :: localrc, stat
    type(type_InternalState)  :: is
    character(ESMF_MAXSTR)    :: name
    type(ESMF_Clock)          :: internalClock, fastClock
    type(ESMF_TimeInterval)   :: couplingStep 
    type(ESMF_TimeInterval)   :: medAtmCouplingStep, medOcnCouplingStep

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! SPECIALIZE by calling into attached method to SetModelServices
    call ESMF_MethodExecute(driver, label=label_SetModelServices, &
      userRc=localrc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out

    ! determine the coupling time steps      
    call ESMF_GridCompGet(driver, clock=internalClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_ClockGet(internalClock, timeStep=couplingStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! query Component for this internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(driver, label_InternalState, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    if (is%wrap%medOcnCouplingIntervalSec>0._ESMF_KIND_R8) then
      ! The coupling time step was provided
      call ESMF_TimeIntervalSet(medOcnCouplingStep, &
        s_r8=is%wrap%medOcnCouplingIntervalSec, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    else
      ! Keep the default timeStep, i.e. that of parent
      medOcnCouplingStep = couplingStep
    endif
    
    if (is%wrap%medAtmCouplingIntervalSec>0._ESMF_KIND_R8) then
      ! The coupling time step was provided
      call ESMF_TimeIntervalSet(medAtmCouplingStep, &
        s_r8=is%wrap%medAtmCouplingIntervalSec, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    else
      ! Use the OCN time step as the ATM default.
      medAtmCouplingStep = medOcnCouplingStep
    endif
    
    ! The NEMS Earth Driver implements a run sequence that supports different
    ! coupling intervals for MED-ATM and MED-OCN coupling. These two coupling
    ! intervals are restricted by the following to constraints:
    ! 1) The MED-ATM coupling is the faster one:
    !      medAtmCouplingStep <= medOcnCouplingStep
    ! 2) The MED-OCN coupling interval must be a multiple of the MED-ATM inverv.
    
    if (medAtmCouplingStep > medOcnCouplingStep) then
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
        msg="The MED-ATM coupling interval must not be larger than the "// &
        "MED-OCN coupling interval!", &
        line=__LINE__, &
        file=__FILE__, rcToReturn=rc)
      return  ! bail out
    
    endif
    
    if (medAtmCouplingStep * (medOcnCouplingStep/medAtmCouplingStep) /= &
      medOcnCouplingStep) then
      call ESMF_LogSetError(ESMF_RC_NOT_VALID, &
        msg="The MED-OCN coupling interval must be a multiple of "// &
        "the MED-ATM coupling interval", &
        line=__LINE__, &
        file=__FILE__, rcToReturn=rc)
      return  ! bail out
    endif
    
    ! Implement the NEMS Earth Driver run sequence, replacing the default.
    call NUOPC_DriverNewRunSequence(driver, slotCount=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! ocn2med into slot 1
    call NUOPC_DriverAddRunElement(driver, slot=1, &
      srcCompLabel="OCN", dstCompLabel="MED", phase=1, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! med slow phase into slot 1
    call NUOPC_DriverAddRunElement(driver, slot=1, &
      compLabel="MED", phase=medPhase_slow, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! med2ocn into slot 1
    call NUOPC_DriverAddRunElement(driver, slot=1, &
      srcCompLabel="MED", dstCompLabel="OCN", phase=1, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! ocn into slot 1
    call NUOPC_DriverAddRunElement(driver, slot=1, &
      compLabel="OCN", phase=1, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! LINK from slot 1 to slot 2
    call NUOPC_DriverAddRunElement(driver, slot=1, linkSlot=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! med fast_before phase into slot 2
    call NUOPC_DriverAddRunElement(driver, slot=2, &
      compLabel="MED", phase=medPhase_fast_before, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! med2atm into slot 2
    call NUOPC_DriverAddRunElement(driver, slot=2, &
      srcCompLabel="MED", dstCompLabel="ATM", phase=1, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! med2ice into slot 2
    call NUOPC_DriverAddRunElement(driver, slot=2, &
      srcCompLabel="MED", dstCompLabel="ICE", phase=1, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! atm into slot 2
    call NUOPC_DriverAddRunElement(driver, slot=2, &
      compLabel="ATM", phase=1, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! ice into slot 2
    call NUOPC_DriverAddRunElement(driver, slot=2, &
      compLabel="ICE", phase=1, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! atm2med into slot 2
    call NUOPC_DriverAddRunElement(driver, slot=2, &
      srcCompLabel="ATM", dstCompLabel="MED", phase=1, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! ice2med into slot 2
    call NUOPC_DriverAddRunElement(driver, slot=2, &
      srcCompLabel="ICE", dstCompLabel="MED", phase=1, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! med fast_after phase into slot 2
    call NUOPC_DriverAddRunElement(driver, slot=2, &
      compLabel="MED", phase=medPhase_fast_after, &
      relaxedflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! Set the slow (MED-OCN) coupling time as time step for the internal clock.
    ! The internal clock is used by the NUOPC Layer to drive slot 1.
    call ESMF_ClockSet(internalClock, timeStep=medOcnCouplingStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! Set the fast (MED-ATM) coupling time step for slot 2
    fastClock = ESMF_ClockCreate(internalClock, rc=rc)  ! make a copy first
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_ClockSet(fastClock, timeStep=medAtmCouplingStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call NUOPC_DriverSetRunSequence(driver, slot=2, clock=fastClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out

    ! Diagnostic output
    call NUOPC_DriverPrint(driver, orderflag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
  end subroutine
    
  !-----------------------------------------------------------------------------
  
  subroutine Finalize(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    
    ! local variables
    integer                   :: localrc, stat
    type(type_InternalState)  :: is
    logical                   :: existflag
    character(ESMF_MAXSTR)    :: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(driver, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! SPECIALIZE by calling into optional attached method
    call ESMF_MethodExecute(driver, label=label_Finalize, existflag=existflag, &
      userRc=localrc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out

    ! query Component for this internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(driver, label_InternalState, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! deallocate internal state memory
    deallocate(is%wrap, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Deallocation of internal state memory failed.", &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out
      
  end subroutine
      
  !-----------------------------------------------------------------------------
  
end module
#endif
