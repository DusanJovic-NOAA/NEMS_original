#include "./ESMFVersionDefine.h"
#ifdef WITH_NUOPC

module module_EARTH_GENERIC_COMP

  !-----------------------------------------------------------------------------
  ! Generic NEMS Earth Driver Component
  !-----------------------------------------------------------------------------

  use ESMF
  use NUOPC
  use NUOPC_Driver, only: &
    Driver_routine_SS             => routine_SetServices, &
    Driver_type_IS                => type_InternalState, &
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
    integer, pointer    :: medPetList(:)
    type(ESMF_GridComp) :: atm
    type(ESMF_GridComp) :: ocn
    type(ESMF_GridComp) :: med
    type(ESMF_State)    :: atmIS, atmES
    type(ESMF_State)    :: ocnIS, ocnES
    type(ESMF_State)    :: medIS, medES
    integer, pointer    :: atm2medPetList(:)
    integer, pointer    :: ocn2medPetList(:)
    integer, pointer    :: med2atmPetList(:)
    integer, pointer    :: med2ocnPetList(:)
    type(ESMF_CplComp)  :: atm2med, ocn2med
    type(ESMF_CplComp)  :: med2atm, med2ocn
    real(ESMF_KIND_R8)  :: medAtmCouplingIntervalSec
    real(ESMF_KIND_R8)  :: medOcnCouplingIntervalSec
  end type

  type type_InternalState
    type(type_InternalStateStruct), pointer :: wrap
  end type
  
  integer, parameter  :: atm=1, ocn=2, med=3

  !-----------------------------------------------------------------------------
  contains
  !-----------------------------------------------------------------------------
  
  subroutine routine_SetServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    character(ESMF_MAXSTR):: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(gcomp, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! NUOPC_Driver registers the generic methods
    call Driver_routine_SS(gcomp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! attach specializing method(s)
    call ESMF_MethodAdd(gcomp, label=Driver_label_SetModelCount, &
      userRoutine=SetModelCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_MethodAdd(gcomp, label=Driver_label_SetModelPetLists, &
      userRoutine=SetModelPetLists, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_MethodAdd(gcomp, label=Driver_label_SetModelServices, &
      userRoutine=SetModelServices, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_MethodAdd(gcomp, label=Driver_label_Finalize, &
      userRoutine=Finalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
  end subroutine
  
  !-----------------------------------------------------------------------------
  
  subroutine SetModelCount(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    type(Driver_type_IS)  :: superIS
    character(ESMF_MAXSTR):: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(gcomp, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! query Component for super internal State
    nullify(superIS%wrap)
    call ESMF_UserCompGetInternalState(gcomp, Driver_label_IS, superIS, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! set the modelCount for ATM-OCN-MED coupling
    superIS%wrap%modelCount = 3
    
  end subroutine
  
  !-----------------------------------------------------------------------------
  
  subroutine SetModelPetLists(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    integer                   :: localrc, stat
    type(type_InternalState)  :: is
    type(Driver_type_IS)      :: superIS
    logical                   :: existflag
    character(ESMF_MAXSTR)    :: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(gcomp, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! allocate memory for this internal state and set it in the Component
    allocate(is%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of internal state memory failed.", &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out
    call ESMF_UserCompSetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! nullify the petLists
    nullify(is%wrap%atmPetList)
    nullify(is%wrap%ocnPetList)
    nullify(is%wrap%medPetList)
    nullify(is%wrap%atm2medPetList)
    nullify(is%wrap%ocn2medPetList)
    nullify(is%wrap%med2atmPetList)
    nullify(is%wrap%med2ocnPetList)
    
    ! SPECIALIZE by calling into optional attached method to set modelPetLists
    call ESMF_MethodExecute(gcomp, label=label_SetModelPetLists, &
      existflag=existflag, userRc=localrc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out

    if (existflag) then
      ! query Component for super internal State
      nullify(superIS%wrap)
      call ESMF_UserCompGetInternalState(gcomp, Driver_label_IS, superIS, rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
      ! set the petLists
      superIS%wrap%modelPetLists(atm)%petList => is%wrap%atmPetList
      superIS%wrap%modelPetLists(ocn)%petList => is%wrap%ocnPetList
      superIS%wrap%modelPetLists(med)%petList => is%wrap%medPetList
      superIS%wrap%connectorPetLists(atm,ocn)%petList => is%wrap%atm2medPetList
      superIS%wrap%connectorPetLists(ocn,med)%petList => is%wrap%ocn2medPetList
      superIS%wrap%connectorPetLists(med,atm)%petList => is%wrap%med2atmPetList
      superIS%wrap%connectorPetLists(med,ocn)%petList => is%wrap%med2ocnPetList
    endif
    
  end subroutine
  
  !-----------------------------------------------------------------------------
  
  subroutine SetModelServices(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    integer                   :: localrc, stat
    type(Driver_type_IS)      :: superIS
    type(type_InternalState)  :: is
    character(ESMF_MAXSTR)    :: name
    type(ESMF_Clock)          :: internalClock, fastClock
    type(ESMF_TimeInterval)   :: couplingStep 
    type(ESMF_TimeInterval)   :: medAtmCouplingStep, medOcnCouplingStep

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(gcomp, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! query Component for super internal State
    nullify(superIS%wrap)
    call ESMF_UserCompGetInternalState(gcomp, Driver_label_IS, superIS, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! allocate memory for this internal state and set it in the Component
    allocate(is%wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Allocation of internal state memory failed.", &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out
    call ESMF_UserCompSetInternalState(gcomp, label_InternalState, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! map components and states for ATM-OCN pair coupling
    is%wrap%atm = superIS%wrap%modelComp(atm)
    is%wrap%atmIS = superIS%wrap%modelIS(atm)
    is%wrap%atmES = superIS%wrap%modelES(atm)
    is%wrap%ocn = superIS%wrap%modelComp(ocn)
    is%wrap%ocnIS = superIS%wrap%modelIS(ocn)
    is%wrap%ocnES = superIS%wrap%modelES(ocn)
    is%wrap%med = superIS%wrap%modelComp(med)
    is%wrap%medIS = superIS%wrap%modelIS(med)
    is%wrap%medES = superIS%wrap%modelES(med)
    is%wrap%atm2med = superIS%wrap%connectorComp(atm,med)
    is%wrap%ocn2med = superIS%wrap%connectorComp(ocn,med)
    is%wrap%med2atm = superIS%wrap%connectorComp(med,atm)
    is%wrap%med2ocn = superIS%wrap%connectorComp(med,ocn)
    
    ! have the component names specified -> makes a difference in the Log
    call ESMF_GridCompSet(is%wrap%atm, name="ATM", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_GridCompSet(is%wrap%ocn, name="OCN", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_GridCompSet(is%wrap%med, name="MED", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_CplCompSet(is%wrap%atm2med, name="ATM2MED", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_CplCompSet(is%wrap%ocn2med, name="OCN2MED", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_CplCompSet(is%wrap%med2atm, name="MED2ATM", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_CplCompSet(is%wrap%med2ocn, name="MED2OCN", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
      
    ! SPECIALIZE by calling into attached method to SetModelServices
    call ESMF_MethodExecute(gcomp, label=label_SetModelServices, &
      userRc=localrc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out

    ! determine the coupling time steps      
    call ESMF_GridCompGet(gcomp, clock=internalClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_ClockGet(internalClock, timeStep=couplingStep, rc=rc)
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
    call NUOPC_RunSequenceDeallocate(superIS%wrap%runSeq, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ! add two run sequence slots: runSeq(1) and runSeq(2)
    call NUOPC_RunSequenceAdd(superIS%wrap%runSeq, 2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out    
    ! ocn2med into slot runSeq(1)
    call NUOPC_RunElementAddComp(superIS%wrap%runSeq(1), i=ocn, j=med, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ! med (phase=2) into slot runSeq(1)
    call NUOPC_RunElementAddComp(superIS%wrap%runSeq(1), i=med, phase=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ! med2ocn into slot runSeq(1)
    call NUOPC_RunElementAddComp(superIS%wrap%runSeq(1), i=med, j=ocn, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ! ocn into slot runSeq(1)
    call NUOPC_RunElementAddComp(superIS%wrap%runSeq(1), i=ocn, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ! LINK slot runSeq(2) into slot runSeq(1)
    call NUOPC_RunElementAddLink(superIS%wrap%runSeq(1), slot=2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ! med2atm into slot runSeq(2)
    call NUOPC_RunElementAddComp(superIS%wrap%runSeq(2), i=med, j=atm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ! atm into slot runSeq(2)
    call NUOPC_RunElementAddComp(superIS%wrap%runSeq(2), i=atm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ! atm2med into slot runSeq(2)
    call NUOPC_RunElementAddComp(superIS%wrap%runSeq(2), i=atm, j=med, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    ! med (phase=1) into slot runSeq(2)
    call NUOPC_RunElementAddComp(superIS%wrap%runSeq(2), i=med, phase=1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! Set the slow (MED-OCN) coupling time as time step for the internal clock.
    ! The internal clock is used by the NUOPC Layer to drive the runSeq(1) slot.
    call ESMF_ClockSet(internalClock, timeStep=medOcnCouplingStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! Set the fast (MED-ATM) coupling time step for the runSeq(2) slot.
    fastClock = ESMF_ClockCreate(internalClock, rc=rc)  ! make a copy first
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call ESMF_ClockSet(fastClock, timeStep=medAtmCouplingStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    call NUOPC_RunSequenceSet(superIS%wrap%runSeq(2), fastClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
  end subroutine
    
  !-----------------------------------------------------------------------------
  
  subroutine Finalize(gcomp, rc)
    type(ESMF_GridComp)  :: gcomp
    integer, intent(out) :: rc
    
    ! local variables
    integer                   :: localrc, stat
    type(type_InternalState)  :: is
    logical                   :: existflag
    character(ESMF_MAXSTR)    :: name

    rc = ESMF_SUCCESS

    ! query the Component for info
    call ESMF_GridCompGet(gcomp, name=name, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    
    ! SPECIALIZE by calling into optional attached method
    call ESMF_MethodExecute(gcomp, label=label_Finalize, existflag=existflag, &
      userRc=localrc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__)) return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=trim(name)//":"//__FILE__, rcToReturn=rc)) &
      return  ! bail out

    ! query Component for this internal State
    nullify(is%wrap)
    call ESMF_UserCompGetInternalState(gcomp, label_InternalState, is, rc)
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
