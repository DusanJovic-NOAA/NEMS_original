!
      module atmos_phy_chem_cpl_comp_mod

!-----------------------------------------------------------------------
!
!** This module holds the phy-to-chem coupler's register and run routines
!** setservices (only registers run step) is called by GFS_ATM_INIT
!** run transfer/convert data from phy export state to chem import state
!** chem-to-phy coupler is not needed (chem gridded component will update
!** 3D aerosol mixing ratios which are owned by dyn gridded component)
!
!! Code Revision:
!! 11Nov 2009     Sarah Lu, First Crack
!! 18Nov 2009     Sarah Lu, Revise coupler run to do data copy
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
!! 11Nov 2009     Sarah Lu, First Crack
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
      MESSAGE_CHECK="Set Entry Point for phy2chem coupler run"

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
        WRITE(0,*)'PHY2CHEM CPL SET_SERVICES SUCCEEDED'
      ELSE
        WRITE(0,*)'PHY2CHEM CPL SET_SERVICES FAILED RC_REG=',RC_REG
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
      subroutine run(GC, PHY_EXP_STATE, CHEM_IMP_STATE, CLOCK, RC_CPL)
!
!-----------------------------------------------------------------------
!!
!! This routine transfer/convert data from phy_exp state to chem_imp state
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

      type(ESMF_cplcomp),intent(inout) :: GC
      type(ESMF_state),  intent(inout) :: PHY_EXP_STATE    ! coupler import state
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

! Fortran array for phy export state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::          &
                p_slmsk, p_hpbl, p_smc1,  p_stype,            &
                p_vtype, p_vfrac,p_rain, p_rainc, p_dtsfci,   &
                p_tsea, p_stc1, p_u10m,  p_v10m,   p_ustar, p_zorl

! Fortran array for chem import state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::          &
                c_lwi,  c_zpbl, c_frlake,  c_fraci,  c_wet1,  &
                c_lai,  c_grn,  c_cn_prcp, c_ncn_prcp, c_sh,  & 
                c_ta,   c_tsoil1, c_u10m,  c_v10m,  c_ustar,  c_z0h

! SOILTYPE and POROSITY
      integer, parameter :: DEFINED_SOIL=9
      real               :: SMCMAX, MAXSMC(DEFINED_SOIL)

      data MAXSMC /0.421, 0.464, 0.468, 0.434, 0.406, 0.465,     &
                   0.404, 0.439, 0.421 /

!---------------------------------------------
!* Get Fortran array from phy export state
!---------------------------------------------

      MESSAGE_CHECK="Get ItemCount from phy export state"

      call esmf_stateget(PHY_EXP_STATE                    &
                        ,itemcount = item_count           &
                        ,itemnamelist = item_name         &
                        ,rc   =rc)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      IF ( RC == ESMF_SUCCESS ) THEN

       WRITE(0,*)'Get Fortran array from phy export state'

       call GetPointer_(PHY_EXP_STATE,'slmsk', p_slmsk, rc)
       call GetPointer_(PHY_EXP_STATE,'hpbl',  p_hpbl,  rc)
       call GetPointer_(PHY_EXP_STATE,'smc1',  p_smc1,  rc)
       call GetPointer_(PHY_EXP_STATE,'stype', p_stype, rc)
       call GetPointer_(PHY_EXP_STATE,'vtype', p_vtype, rc)
       call GetPointer_(PHY_EXP_STATE,'vfrac', p_vfrac, rc)
       call GetPointer_(PHY_EXP_STATE,'rain',  p_rain,  rc)
       call GetPointer_(PHY_EXP_STATE,'rainc', p_rainc, rc)
       call GetPointer_(PHY_EXP_STATE,'dtsfci',p_dtsfci,rc)
       call GetPointer_(PHY_EXP_STATE,'tsea',  p_tsea,  rc)
       call GetPointer_(PHY_EXP_STATE,'stc1',  p_stc1,  rc)
       call GetPointer_(PHY_EXP_STATE,'u10m',  p_u10m,  rc)
       call GetPointer_(PHY_EXP_STATE,'v10m',  p_v10m,  rc)
       call GetPointer_(PHY_EXP_STATE,'ustar', p_ustar, rc)
       call GetPointer_(PHY_EXP_STATE,'zorl',  p_zorl,  rc)

      ELSE

       WRITE(0,*)'Empty phy export state; Fortran array not filled'

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
       call GetPointer_(CHEM_IMP_STATE,'lwi',   c_lwi,   rc)
       call GetPointer_(CHEM_IMP_STATE,'zpbl',  c_zpbl,  rc)
       call GetPointer_(CHEM_IMP_STATE,'frlake',c_frlake,rc)
       call GetPointer_(CHEM_IMP_STATE,'fraci', c_fraci, rc)
       call GetPointer_(CHEM_IMP_STATE,'wet1',  c_wet1,  rc)
       call GetPointer_(CHEM_IMP_STATE,'lai',   c_lai,   rc)
       call GetPointer_(CHEM_IMP_STATE,'grn',   c_grn,   rc)
       call GetPointer_(CHEM_IMP_STATE,'cn_prcp', c_cn_prcp,rc)
       call GetPointer_(CHEM_IMP_STATE,'ncn_prcp',c_ncn_prcp,rc)
       call GetPointer_(CHEM_IMP_STATE,'sh',    c_sh,    rc)
       call GetPointer_(CHEM_IMP_STATE,'ta',    c_ta,    rc)
       call GetPointer_(CHEM_IMP_STATE,'tsoil1',c_tsoil1,rc)
       call GetPointer_(CHEM_IMP_STATE,'u10m',  c_u10m,  rc)
       call GetPointer_(CHEM_IMP_STATE,'v10m',  c_v10m,  rc)
       call GetPointer_(CHEM_IMP_STATE,'ustar', c_ustar, rc)
       call GetPointer_(CHEM_IMP_STATE,'z0h',   c_z0h,   rc)

      ELSE

       WRITE(0,*)'Empty chem import state; Fortran array not filled'

      ENDIF

!---------------------------------------------
!* Do the actual coupling (data copy)
!---------------------------------------------
! 

      IF ( RC_CPL == ESMF_SUCCESS ) THEN

       WRITE(0,*)'Coupling phy_exp with chem_imp succeeded'

       c_frlake = 0.                  ! fraction_of_lake (1)
       c_fraci  = 0.                  ! ice_covered_fraction_of_tile (1)
       c_lai    = 3.                  ! leaf_area_index (1)
!
       c_zpbl   = p_hpbl              ! boundary layer height (m)
       c_grn    = p_vfrac             ! greeness_fraction (1)
       c_sh     = p_dtsfci            ! sensible heat flux (W/m^2)
       c_ta     = p_tsea              ! surface air Temperature (K)
       c_tsoil1 = p_stc1              ! soil temperatures layer_1 (k)
       c_u10m   = p_u10m              ! 10-meter eastward_wind (m s-1)
       c_v10m   = p_v10m              ! 10-meter northward_wind (m s-1)
       c_ustar  = p_ustar             ! surface velocity scale (m s-1)
!
       c_cn_prcp  = 1.E3*p_rainc      ! surface conv. rain flux (kg/m^2/s)
       c_ncn_prcp = 1.E3*(p_rain - p_rainc) ! Non-conv. precip rate (kg/m^2/s)
       c_z0h      = p_zorl / 1.E2     ! surface roughness (m)
! 
       c_lwi = p_slmsk                ! land-ocean-ice mask  (1)
!
       do i = 1, size(p_slmsk, dim=1)
       do j = 1, size(p_slmsk, dim=2)
         if ( p_slmsk(i,j) == 1. ) then      ! <--- over ocean
           smcmax = maxsmc(p_vtype(i,j))     ! porosity
           c_wet1(i,j) = p_smc1(i,j)/smcmax
         else
           c_wet1(i,j) = 0.
         endif
       enddo
       enddo
      ELSE

        WRITE(0,*)'Coupling phy_exp with chem_imp failed RC_CPL=',RC_CPL

      ENDIF

!
!! 
      contains 

        subroutine GetPointer_ (STATE, NAME, ARRAY, RC)

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
        integer(ESMF_KIND_I4),allocatable,save :: lonsperlar_r(:)
        integer, save                 :: lonr, lats_node_r, lats_node_r_max
        logical, save                 :: first

        data first /.true./
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
        print *, 'LU_DEBUG(PHY):', Array(1,1),Array(1,2),Array(2,1)
! 
        if ( first ) then

         MESSAGE_CHECK = 'Extract lonr Attribute'
         CALL ESMF_AttributeGet(state=STATE, name='lonr'  &
                              ,value=lonr, rc=RC1)
         CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

         MESSAGE_CHECK = 'Extract lats_node_r Attribute'
         CALL ESMF_AttributeGet(state=STATE, name='lats_node_r'  &  
                              ,value=lats_node_r, rc=RC1)
         CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

         MESSAGE_CHECK = 'Extract lats_node_r_max Attribute'
         CALL ESMF_AttributeGet(state=STATE, name='lats_node_r_max'&
                              ,value=lats_node_r_max, rc=RC1)
         CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

         if ( .not. allocated (lonsperlar_r)) then
           allocate ( lonsperlar_r(lats_node_r_max))
         endif

         MESSAGE_CHECK = 'Extract lonsperlar_r Attribute'
         CALL ESMF_AttributeGet(state=STATE             &  !<-- The physics export state
                            ,name ='lonsperlar_r'       &  !<-- Name of the attribute to insert
                            ,count = lats_node_r_max    &  !<-- Number of values in the attribute
                            ,valueList =lonsperlar_r    &  !<-- Value of the attribute
                            ,rc   =RC1)
         CALL ERR_MSG(rc1, MESSAGE_CHECK, rc)

         if ( lonr .ne. size ( Array, dim=1) )  print *, 'ERROR !',       &
             'Invalid lonr:',  lonr, size(Array,dim=1)
         if (lats_node_r_max .ne. size ( Array, dim=2) ) print *, 'ERROR !', &
             'Invalid lats_node_r_max',  lats_node_r_max, size(Array,dim=2)

         first = .False.      ! reset first                         
        endif

! patch the array with valid values
        do i = 1, lats_node_r_max
           if (lonsperlar_r(i) < lonr ) then
             Array(i,lonsperlar_r(i)+1:lonr) = Array(1,1)
           endif
        enddo 

       end subroutine GetPointer_
!
      END subroutine run


      END module atmos_phy_chem_cpl_comp_mod



