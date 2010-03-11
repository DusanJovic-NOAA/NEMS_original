!
      module atmos_chem_phy_cpl_comp_mod

!-----------------------------------------------------------------------
!
!** This module holds the chem-to-phy coupler's register and run routines
!** setservices (only registers run step) is called by GFS_ATM_INIT
!** run transfer/convert data from chem export state to phy export state
!
!! Code Revision:
!! 24Feb 2009     Sarah Lu, First Crack
!-----------------------------------------------------------------------

      use ESMF_MOD

      USE MODULE_ERR_MSG, ONLY: ERR_MSG, MESSAGE_CHECK
      use MODULE_gfs_machine,  ONLY: kind_phys

!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      private

!     Gaussian grid 
      integer(ESMF_KIND_I4),allocatable,save :: lonsperlar_r(:)
      integer, save                :: lonr, lats_node_r, lats_node_r_max
      integer, save                :: im, jm, km

!     Tracer specification
      integer, save   :: ntrac
      logical, save   :: run_DU, run_SU, run_SS, run_OC, run_BC
      real(kind_phys), allocatable :: ri(:), cpi(:)
      
   
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
!! 24Feb 2009     Sarah Lu, First Crack
!! register init
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
      MESSAGE_CHECK="Set Entry Point for chem2phy coupler run"

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
        WRITE(0,*)'CHEM2PHY CPL SET_SERVICES SUCCEEDED'
      ELSE
        WRITE(0,*)'CHEM2PHY CPL SET_SERVICES FAILED RC_REG=',RC_REG
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
      subroutine run(GC, CHEM_EXP_STATE, PHY_EXP_STATE, CLOCK, RC_CPL)
!
!-----------------------------------------------------------------------
!!
!! This routine transfer/convert data from chem_exp state to phy_exp state
!! ===> flip vertical profile index for chemical tracer arrays
!!
!! Code Revision:
!! 24Feb 2010     Sarah Lu, First Crack
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------

      type(ESMF_cplcomp),intent(inout) :: GC
      type(ESMF_state),  intent(inout) :: CHEM_EXP_STATE   ! coupler import state
      type(ESMF_state),  intent(inout) :: PHY_EXP_STATE    ! coupler export state
      type(ESMF_clock),  intent(in)    :: CLOCK
!
      integer,           intent(out)   :: RC_CPL

!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer                 :: rc=ESMF_success  ! the error signal variable
      integer                 :: i, j, k, item_count_phys, item_count_chem, kp
      character(20)           :: item_name(200)
      logical, save           :: first =  .true.
      real(ESMF_KIND_R8), pointer     :: Array(:,:,:)
      type(ESMF_Field)                :: Field

! Fortran array for phy export state
      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::                  &
!              p_spfh, p_o3mr,                                           & ! met tracers
               p_du001, p_du002, p_du003, p_du004, p_du005,              & ! DU
               p_ss001, p_ss002, p_ss003, p_ss004, p_ss005,              & ! SS
               p_msa,   p_so4,   p_so2,   p_dms,                         & ! SU
               p_ocphobic, p_ocphilic,                                   & ! OC
               p_bcphobic, p_bcphilic                                      ! BC
           

! Fortran array for chem export state
      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::                  &
!              c_o3, c_rh,                                               & ! met tracers
               c_du001, c_du002, c_du003, c_du004, c_du005,              & ! DU
               c_ss001, c_ss002, c_ss003, c_ss004, c_ss005,              & ! SS
               c_msa,   c_so4,   c_so2,   c_dms,                         & ! SU
               c_ocphobic, c_ocphilic,                                   & ! OC
               c_bcphobic, c_bcphilic                                      ! BC

!---------------------------------------------
!* Determine dimension and allocate local array
!---------------------------------------------
!
      IF ( first ) THEN
!
!  --- Retrieve attributes (lat/lon and tracer specification) 
!  --- Fill in tracer_const local arrays
        call get_attribute
!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)
!
        call ESMF_StateGet(state      = PHY_EXP_STATE       &
                          ,itemName   = 't'                &
                          ,field      = Field               &
                          ,rc   =rc)
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_CPL)
!
        nullify(Array)
        MESSAGE_CHECK = 'Get Fortran data pointer from t'
        CALL ESMF_FieldGet(field=Field, localDe=0, &
                           farray=Array, rc = rc)
        CALL ERR_MSG(RC, MESSAGE_CHECK, RC_CPL)

!       determine dimension in x-, y-, and z-direction
        im = size(Array, dim=1)
        jm = size(Array, dim=2)
        km = size(Array, dim=3)

! debug print
        print *, 'LU_TSTx: cpi =', cpi(0:ntrac)
        print *, 'LU_TSTx: ri  =', ri(0:ntrac)
        print *, 'LU_TSTx: ntrac =', ntrac
        print *, 'LU_TSTx: doing_DU =', run_DU
        print *, 'LU_TSTx: doing_SU =', run_SU
        print *, 'LU_TSTx: doing_SS =', run_SS
        print *, 'LU_TSTx: doing_OC =', run_OC
        print *, 'LU_TSTx: doing_BC =', run_BC
	print *, 'LU_TSTx: lonr =', lonr
	print *, 'LU_TSTx: lats_node_r =', lats_node_r
	print *, 'LU_TSTx: lats_node_r_max =', lats_node_r_max
	print *, 'LU_TSTx: lonsperlar_r =', lonsperlar_r(:)
        print *, 'LU_TSTx: im, jm, km=', im, jm, km

!       reset first flag
        first = .false.

      ENDIF

!---------------------------------------------
!* Get Fortran array from phy export state
!---------------------------------------------

      MESSAGE_CHECK="Get ItemCount from phy export state"

      call esmf_stateget(PHY_EXP_STATE                    &
                        ,itemcount = item_count_phys      &
                        ,itemnamelist = item_name         &
                        ,rc   =rc)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      print *,'LU_TST: CPL: phy_exp item count:',item_count_phys
      print *,'LU_TST: CPL: phy_exp item name :',      &
                       (item_name(i),i=1,item_count_phys)

      IF ( RC == ESMF_SUCCESS ) THEN

      if (item_count_phys == 0 ) then

       WRITE(0,*)'Empty phy export state; Fortran array not filled'

      else

       WRITE(0,*)'Get Fortran array from phy export state'

! get met tracers
!      call GetPointer_tracer_(PHY_EXP_STATE,'spfh', p_spfh, rc)
!      call GetPointer_tracer_(PHY_EXP_STATE,'o3mr', p_o3mr, rc)

! get chem tracers
       if ( run_DU ) then
        call GetPointer_tracer_(PHY_EXP_STATE,'du001', p_du001, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'du002', p_du002, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'du003', p_du003, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'du004', p_du004, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'du005', p_du005, rc)
       endif 

       if ( run_SS ) then
        call GetPointer_tracer_(PHY_EXP_STATE,'ss001', p_ss001, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'ss002', p_ss002, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'ss003', p_ss003, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'ss004', p_ss004, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'ss005', p_ss005, rc)
       endif

       if ( run_SU ) then
        call GetPointer_tracer_(PHY_EXP_STATE,'msa', p_msa, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'so4', p_so4, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'so2', p_so2, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'dms', p_dms, rc)
       endif

       if ( run_OC ) then
        call GetPointer_tracer_(PHY_EXP_STATE,'ocphobic', p_ocphobic, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'ocphilic', p_ocphilic, rc)
       endif

       if ( run_BC ) then
        call GetPointer_tracer_(PHY_EXP_STATE,'bcphobic', p_bcphobic, rc)
        call GetPointer_tracer_(PHY_EXP_STATE,'bcphilic', p_bcphilic, rc)
       endif

      endif

      ELSE

       WRITE(0,*)'Fail to get phy export state'

      ENDIF


!---------------------------------------------
!* Get Fortran array from gocart export state
!---------------------------------------------

      MESSAGE_CHECK="Get ItemCount from chem export state"

      call esmf_stateget(CHEM_EXP_STATE                   &
                        ,itemcount = item_count_chem      &
                        ,itemnamelist = item_name         &
                        ,rc   =rc)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      print *,'LU_TST: CPL: chem_exp item count:',item_count_chem
      print *,'LU_TST: CPL: chem_exp item name :',      &
                       (item_name(i),i=1,item_count_chem)

      IF ( RC == ESMF_SUCCESS ) THEN

      if (item_count_chem == 0 ) then

       WRITE(0,*)'Empty chem export state; Fortran array not filled'

      else
       WRITE(0,*)'Get Fortran array from chem export state'
!      call GetPointer_3D_(CHEM_IMP_STATE,'O3' , c_o3  , rc)

      endif

      ELSE

       WRITE(0,*)'Fail to get chem export state'

      ENDIF


!---------------------------------------------
!* Do the actual coupling (data copy)
!---------------------------------------------
! 
      IF ( RC_CPL == ESMF_SUCCESS ) THEN

       if (item_count_chem == 0 ) then

        WRITE(0,*)'Empty chem exp; Coupling chem_exp with phy_exp skipped'

       else

        WRITE(0,*)'Chem exp not empty; Coupling chem_exp with phy_exp'
!
! ---  data copy between phy export state and chem export state
!
        do j = 1, lats_node_r
        do i = 1, lonsperlar_r(j)

          do k = 1, km
            kp = km - k + 1
            if ( run_DU ) then
              p_du001(i,j,kp) = c_du001(i,j,k) 
              p_du002(i,j,kp) = c_du002(i,j,k) 
              p_du003(i,j,kp) = c_du003(i,j,k) 
              p_du004(i,j,kp) = c_du004(i,j,k) 
              p_du005(i,j,kp) = c_du005(i,j,k) 
            endif

            if ( run_SS ) then
              p_ss001(i,j,kp) = c_ss001(i,j,k) 
              p_ss002(i,j,kp) = c_ss002(i,j,k) 
              p_ss003(i,j,kp) = c_ss003(i,j,k) 
              p_ss004(i,j,kp) = c_ss004(i,j,k) 
              p_ss005(i,j,kp) = c_ss005(i,j,k) 
            endif

            if ( run_SU ) then
              p_msa(i,j,kp) = c_msa(i,j,k) 
              p_so2(i,j,kp) = c_so2(i,j,k) 
              p_so4(i,j,kp) = c_so4(i,j,k) 
              p_dms(i,j,kp) = c_dms(i,j,k) 
            endif

            if ( run_OC ) then
              p_ocphobic(i,j,kp) = c_ocphobic(i,j,k) 
              p_ocphilic(i,j,kp) = c_ocphilic(i,j,k) 
            endif

            if ( run_BC ) then
              p_bcphobic(i,j,kp) = c_bcphobic(i,j,k) 
              p_bcphilic(i,j,kp) = c_bcphilic(i,j,k) 
            endif

          enddo

        enddo
        enddo

       endif

      ELSE

       WRITE(0,*)'Coupling phy_exp with chem_exp skipped, RC_CPL=',RC_CPL

      ENDIF

!
!! 
      contains 


!=========================
      subroutine get_attribute

!  ---  inputs:  (in scope variables)
!  ---  outputs: (in scope variables)

!  ---  locals:
        type(ESMF_FieldBundle)   :: Bundle
        TYPE(ESMF_Logical)       :: doing_SU, doing_SS, &
                                    doing_DU, doing_OC, doing_BC
        integer                  :: status

! ---  Retrieve lat/lon info
        MESSAGE_CHECK = 'Extract lonr attribute from phy_exp'
        CALL ESMF_AttributeGet(state=PHY_EXP_STATE, name='lonr'  &
                              ,value=lonr, rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract lats_node_r attribute from phy_exp'
        CALL ESMF_AttributeGet(state=PHY_EXP_STATE, name='lats_node_r'  &  
                             ,value=lats_node_r, rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract lats_node_r_max attribute from phy_exp'
        CALL ESMF_AttributeGet(state=PHY_EXP_STATE, name='lats_node_r_max'&
                             ,value=lats_node_r_max, rc=RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        if ( .not. allocated (lonsperlar_r)) then
          allocate ( lonsperlar_r(lats_node_r_max))
        endif

        MESSAGE_CHECK = 'Extract lonsperlar_r attribute from phy_exp'
        CALL ESMF_AttributeGet(state=PHY_EXP_STATE     &  !<-- The physics export state
                           ,name ='lonsperlar_r'       &  !<-- Name of the attribute to retrieve
                           ,count = lats_node_r_max    &  !<-- Number of values in the attribute
                           ,valueList =lonsperlar_r    &  !<-- Value of the attribute
                           ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

! ---  Retrieve tracer specification
        MESSAGE_CHECK = 'Extract tracer bundle from phy_exp'
        call ESMF_StateGet(state=PHY_EXP_STATE, ItemName='tracers',    &
                          fieldbundle=Bundle, rc = rc )
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract ntrac attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                        &  !<-- Tracer bundle
               ,name  ='ntrac'                               &  !<-- Name of the attribute to retrieve
               ,value = ntrac                                &  !<-- Value of the attribute
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        allocate(ri(0:ntrac), stat=status)
           if( status .ne. 0 ) print *, 'ERROR: Fail to allocate ri'
        allocate(cpi(0:ntrac), stat=status)
           if( status .ne. 0 ) print *, 'ERROR: Fail to allocate cpi'

        MESSAGE_CHECK = 'Extract cpi_dryair attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                        &  !<-- Tracer bundle
               ,name ='cpi_dryair'                           &  !<-- Name of the attribute to retrieve
               ,value = cpi(0)                               &  !<-- Value of the attribute
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract ri_dryair attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                        &  !<-- Tracer bundle
               ,name  ='ri_dryair'                           &  !<-- Name of the attribute to retrieve
               ,value = ri(0)                                &  !<-- Value of the attribute
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract cpi attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                        &  !<-- Tracer bundle
               ,name  ='cpi'                                 &  !<-- Name of the attribute to retrieve
               ,count= ntrac                                 &  !<-- Number of values in the attribute
               ,valueList = cpi(1:ntrac)                     &  !<-- Value of the attribute
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        MESSAGE_CHECK = 'Extract ri attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                        &  !<-- Tracer bundle
               ,name  ='ri'                                  &  !<-- Name of the attribute to retrieve
               ,count= ntrac                                 &  !<-- Number of values in the attribute
               ,valueList = ri(1:ntrac)                      &  !<-- Value of the attribute
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

        run_DU = .False.
        MESSAGE_CHECK = 'Extract doing_DU attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name  ='doing_DU'                             &  !<-- Name of the logical
               ,value = doing_DU                              &  !<-- The logical being retrieved
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)
        if ( doing_DU == ESMF_True ) run_DU = .True.

        run_SU = .False.
        MESSAGE_CHECK = 'Extract doing_SU attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name  ='doing_SU'                             &  !<-- Name of the logical
               ,value = doing_SU                              &  !<-- The logical being retrieved
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)
        if ( doing_SU == ESMF_True ) run_SU = .True.

        run_SS = .False.
        MESSAGE_CHECK = 'Extract doing_SS attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name  ='doing_SS'                             &  !<-- Name of the logical
               ,value = doing_SS                              &  !<-- The logical being retrieved
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)
        if ( doing_SS == ESMF_True ) run_SS = .True.

        run_OC = .False.
        MESSAGE_CHECK = 'Extract doing_OC attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name  ='doing_OC'                             &  !<-- Name of the logical
               ,value = doing_OC                              &  !<-- The logical being retrieved
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)
        if ( doing_OC == ESMF_True ) run_OC = .True.

        run_BC = .False.
        MESSAGE_CHECK = 'Extract doing_BC attribute from phy_exp tracer bundle'
        CALL ESMF_AttributeGet(Bundle                         &  !<-- Phy export state tracer bundle
               ,name  ='doing_BC'                             &  !<-- Name of the logical
               ,value = doing_BC                              &  !<-- The logical being retrieved
               ,rc   =RC)
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)
        if ( doing_BC == ESMF_True ) run_BC = .True.

       RETURN
      end subroutine get_attribute

        subroutine GetPointer_ (STATE, NAME, Array, RC)

! --- input/output arguments
        type(ESMF_State), intent(IN)    :: State
        character(len=*), intent(IN)    :: Name
        real(ESMF_KIND_R8), pointer, intent(OUT)  :: Array(:,:)
        integer, intent (OUT)           :: rc

! --- locals
        type(ESMF_Field)                :: Field
        integer                         :: i, j, rc1
        character(esmf_maxstr)          :: statename
!
!===>  ...  begin here
!
        call ESMF_StateGet(state=State, name=statename, rc=rc )

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
!
        print *,'LU_TST:',trim(statename),Name,Array(1,1),Array(1,2),Array(2,1)
! 
!       check im and jm for physics export
        if (statename == 'Physics Export') then
          if ( lonr .ne. size ( Array, dim=1) )  print *, 'ERROR !',       &
             'Invalid lonr:',  lonr, size(Array,dim=1)
          if (lats_node_r_max .ne. size ( Array, dim=2) ) print *, 'ERROR !', &
             'Invalid lats_node_r_max',  lats_node_r_max, size(Array,dim=2)
        endif

!       patch the array with valid values for chemistry import
        if (statename == 'Chemistry Import') then
          do j = 1, lats_node_r_max
            if (lonsperlar_r(j) < lonr ) then
              Array(lonsperlar_r(j)+1:lonr,j) = Array(1,1)
              if(NAME=='smc1') print *,'LU_TST: fill smc1:',j,lonr,Array(1,1)
             endif
          enddo 
        endif

       end subroutine GetPointer_
!!!
        subroutine GetPointer_3D_ (STATE, NAME, Array, RC)

! --- input/output arguments
        type(ESMF_State), intent(IN)    :: State
        character(len=*), intent(IN)    :: Name
        real(ESMF_KIND_R8), pointer, intent(OUT)  :: Array(:,:,:)
        integer, intent (OUT)           :: rc

! --- locals
        type(ESMF_Field)                :: Field
        integer                         :: i, j, k, rc1
        character(esmf_maxstr)          :: statename
!
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

        print *, 'LU_TST:',trim(statename),Name,': ',Array(1,1,1),   &
                           Array(1,2,1), Array(im,jm,km)

!       check im, jm, km for physics export
        if (statename == 'Physics Export') then
          if ( lonr .ne. size ( Array, dim=1) )  print *, 'ERROR !',       &
             'Invalid lonr:',  lonr, size(Array,dim=1)
          if (lats_node_r_max .ne. size ( Array, dim=2) ) print *, 'ERROR !', &
             'Invalid lats_node_r_max',  lats_node_r_max, size(Array,dim=2)
          if (km .ne. size ( Array, dim=3) ) print *, 'ERROR !', &
             'Invalid km',  km, size(Array,dim=3)
        endif

!       patch the array with valid values
        if (statename == 'Chemistry Import') then
          do k = 1, km
          do j = 1, lats_node_r_max
            if (lonsperlar_r(j) < lonr ) then
              Array(lonsperlar_r(j)+1:lonr,j,k) = Array(1,1,k)
              if(NAME=='t') print *,'LU_TST: fill t:',i,lonr,Array(1,1,1)
             endif
          enddo 
          enddo 
        endif

       end subroutine GetPointer_3D_

!!!!
        subroutine GetPointer_tracer_ (STATE, NAME, Array, RC)

! --- input/output arguments
        type(ESMF_State), intent(IN)    :: State
        character(len=*), intent(IN)    :: Name
        real(ESMF_KIND_R8), pointer, intent(OUT)  :: Array(:,:,:)
        integer, intent (OUT)           :: rc

! --- locals
        type(ESMF_Field)                :: Field
        type(ESMF_FieldBundle)          :: Bundle
        integer                         :: i, j, k, rc1
        character(esmf_maxstr)          :: statename
!
!===>  ...  begin here
!
        call ESMF_StateGet(state=State, name=statename, rc=rc )

        MESSAGE_CHECK = 'Extract tracer bundle from '//statename
        call ESMF_StateGet(state=State, ItemName='tracers', &
                           fieldbundle=Bundle, rc = rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
        MESSAGE_CHECK = 'Extract '//NAME//' from tracer bundle'
        CALL ESMF_FieldBundleGet(bundle=Bundle, &
                      name=NAME, field=Field, rc = rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)
!
        nullify(Array)
        MESSAGE_CHECK = 'Get Fortran data pointer from '//NAME
        CALL ESMF_FieldGet(field=Field, localDe=0, &
                           farray=Array, rc = rc1)
        CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

        print *, 'LU_TST',trim(statename),Name,Array(1,1,1), &
                          Array(1,2,1), Array(im,jm,km)
!
! check im and jm
        if ( lonr .ne. size ( Array, dim=1) )  print *, 'ERROR !',       &
             'Invalid lonr:',  lonr, size(Array,dim=1)
        if (lats_node_r_max .ne. size ( Array, dim=2) ) print *, 'ERROR !', &
             'Invalid lats_node_r_max',  lats_node_r_max, size(Array,dim=2)
        if (km .ne. size ( Array, dim=3) ) print *, 'ERROR !', &
             'Invalid km',  km, size(Array,dim=3)

       end subroutine GetPointer_tracer_

!
      END subroutine run

      END module atmos_chem_phy_cpl_comp_mod



