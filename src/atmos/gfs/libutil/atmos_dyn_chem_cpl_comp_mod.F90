!
      module atmos_dyn_chem_cpl_comp_mod

!-----------------------------------------------------------------------
!**
!** This module holds the dyn-to-chem coupler's register and run routines
!** setservices (only registers run step) is called by GFS_ATM_INIT
!** run transfer/convert data from dyn export state to chem import state
!** chem-to-dyn coupler is not needed (chem gridded component will update
!** 3D aerosol mixing ratios which are owned by dyn gridded component)
!
!
!! Code Revision:
!! 18Nov 2009     Sarah Lu, First Crack
!! 29Dec 2009     Sarah Lu, Comments added for clarification
!-----------------------------------------------------------------------

      use ESMF_MOD

      USE MODULE_ERR_MSG, ONLY: ERR_MSG, MESSAGE_CHECK

!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      private

      public :: setservices

      contains

!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine setservices(GC, RC_REG)
!
!-----------------------------------------------------------------------
!!
!! This routine register the coupler component's run routine
!!
!! Code Revision:
!! 18Nov 2009     Sarah Lu, First Crack
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      type(ESMF_cplcomp),intent(inout) :: gc         ! coupler component
!
      integer,intent(out) :: rc_reg                  ! return code for register
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer :: rc=ESMF_success                     ! the error signal variable

!-----------------------------------------------------------------------
!***  register the coupler component's run routine
!-----------------------------------------------------------------------
!
      MESSAGE_CHECK="Set Entry Point for dyn2chem coupler run"

      call ESMF_CplCompSetEntryPoint(GC                        & !<-- The gridded component
                                    ,ESMF_SETRUN               & !<-- Predefined subroutine type
                                    ,RUN                       & !<-- User's subroutineName
                                    ,ESMF_SINGLEPHASE          &
                                    ,rc)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_REG)

!-----------------------------------------------------------------------
!***  Check the error signal variable and print out the result.
!-----------------------------------------------------------------------
!
      IF(RC_REG==ESMF_SUCCESS)THEN
        WRITE(0,*)'DYN2CHEM CPL SET_SERVICES SUCCEEDED'
      ELSE
        WRITE(0,*)'DYN2CHEM CPL SET_SERVICES FAILED RC_REG=',RC_REG
      ENDIF

!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      end subroutine setservices


!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      subroutine run(GC, DYN_EXP_STATE, CHEM_IMP_STATE, CLOCK, RC_CPL)
!
!-----------------------------------------------------------------------
!!
!! This routine transfer/convert data from dyn_exp state to chem_imp state
!!
!! Code Revision:
!! 18Nov 2009     Sarah Lu, First Crack
!! 29Dec 2009     Sarah Lu, Comments added for clarification
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      type(ESMF_cplcomp),intent(inout) :: GC
      type(ESMF_state),  intent(inout) :: DYN_EXP_STATE    ! coupler import state
      type(ESMF_state),  intent(inout) :: CHEM_IMP_STATE   ! coupler export state
      type(ESMF_clock),  intent(in)    :: CLOCK
!
      integer,           intent(out)   :: RC_CPL

!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer                 :: rc=ESMF_success  ! the error signal variable
      integer                 :: i, j, k
      integer                 :: item_count
      character(20)           :: item_name(200)

! Fortran array for dyn export state
      real (ESMF_KIND_R8), pointer, dimension(:,:) :: d_ps
      real (ESMF_KIND_R8), pointer, dimension(:,:,:) :: d_u, d_v, d_t

! Fortran array for chem import state
      real (ESMF_KIND_R8), pointer, dimension(:,:) :: c_ps
      real (ESMF_KIND_R8), pointer, dimension(:,:,:) :: c_u, c_v, c_t

!---------------------------------------------
!* Get Fortran array from dyn export state
!---------------------------------------------

      MESSAGE_CHECK="Get ItemCount from dyn export state"

      call esmf_stateget(DYN_EXP_STATE                    &
                        ,itemcount = item_count           &
                        ,itemnamelist = item_name         &
                        ,rc   =rc)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      IF ( RC == ESMF_SUCCESS ) THEN
  
       WRITE(0,*)'Get Fortran array from dyn export state'
       call GetPointer_2D_(DYN_EXP_STATE,'ps', d_ps, rc)
 
      ELSE

       WRITE(0,*)'Empty dyn export state; Fortran array not filled'

      ENDIF

!---------------------------------------------
!* Get Fortran array from gocart import state
!---------------------------------------------

      MESSAGE_CHECK="Get ItemCount from chem import state"

      call esmf_stateget(CHEM_IMP_STATE                   &
                        ,itemcount = item_count           &
                        ,itemnamelist = item_name         &
                        ,rc   =rc)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      IF ( RC == ESMF_SUCCESS ) THEN
  
       WRITE(0,*)'Get Fortran array from chem import state'
       call GetPointer_2D_(CHEM_IMP_STATE,'PS',  c_ps,   rc)
 
      ELSE

       WRITE(0,*)'Empty chem import state; Fortran array not filled'

      ENDIF

!---------------------------------------------
!* Do the actual coupling (data copy)
!---------------------------------------------
! 
       IF ( RC_CPL == ESMF_SUCCESS ) THEN

        WRITE(0,*)'Coupling dyn_exp with chem_imp succeeded'
        c_ps = d_ps

       ELSE

        WRITE(0,*)'Coupling dyn_exp with chem_imp failed RC_CPL=',RC_CPL

       ENDIF

!! 
      contains 

        subroutine GetPointer_2D_ (STATE, NAME, ARRAY, RC)

! --- input/output arguments
        type(ESMF_State), intent(IN)    :: State
        character(len=*), intent(IN)    :: Name
        real(ESMF_KIND_R8), pointer     :: Array(:,:)
        integer, intent (OUT)           :: rc

! --- locals
        type(ESMF_Field)              :: Field
        integer                       :: i, j, dim1, dim2, rc1
        character(esmf_maxstr)        :: statename
!
        logical, save                 :: first

        data first /.true./
!===>  ...  begin here
!
        call ESMF_StateGet(state=State, name=statename, rc=rc ) 
!
        MESSAGE_CHECK = 'Extract '//NAME//' from '//statename
        call ESMF_StateGet(state      = STATE               &
                          ,itemName   = NAME                &
                          ,field      = Field               &
                          ,rc   =rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
        nullify(Array)
        MESSAGE_CHECK = 'Get Fortran data pointer from '//NAME
        CALL ESMF_FieldGet(field=Field, localDe=0, &
                           farray=Array, rc = rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

	print *, 'LU_DEBUG(DYN):', Array(1,1),Array(1,2),Array(2,1)
!
! 
       end subroutine GetPointer_2D_
!
      END subroutine run


      END module atmos_dyn_chem_cpl_comp_mod



