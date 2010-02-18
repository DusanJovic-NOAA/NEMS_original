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
!! 01Feb 2010     Sarah Lu, Extend to include all 2d/3d fields needed by
!!                          GOCART
!! 07Feb 2010     Sarah Lu, Add getrh subroutine
!! 11Feb 2010     Sarah Lu, Add get_attribute subroutine to retrieve tracer
!!                          specification 
!! 12Feb 2010     Sarah Lu, Include chemical tracers
!-----------------------------------------------------------------------

      use ESMF_MOD

      USE MODULE_ERR_MSG, ONLY: ERR_MSG, MESSAGE_CHECK
      use MODULE_gfs_machine,  ONLY: kind_phys
      use MODULE_gfs_funcphys, ONLY: fpvs, krealfp
      use MODULE_gfs_physcons, ONLY: con_rd,  con_fvirt, con_g, &
                                     con_eps, con_epsm1
      use MODULE_gfs_tropp,    ONLY: tpause


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
!! 01Feb 2010     Sarah Lu, Transfer/convert fields from phy export state 
!!                          to GOCART import state
!! 11Feb 2010     Sarah Lu, Add get_attribute
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
      integer                 :: i, j, k, item_count_phys, item_count_chem
      character(20)           :: item_name(200)
      logical, save           :: first =  .true.
      real(ESMF_KIND_R8), pointer     :: Array(:,:,:)
      type(ESMF_Field)                :: Field

! Fortran array for phy export state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::                    &
               p_slmsk, p_hpbl,  p_smc1,  p_stype,  p_vtype, p_vfrac,   &
               p_rain,  p_rainc, p_dtsfci,p_tsea,   p_stc1,  p_u10m,    &
               p_v10m,  p_ustar, p_zorl,  p_hs,     p_ps

      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::                  &
               p_t, p_u, p_v, p_p, p_dp, p_fcld, p_dqdt

      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::                  &
               p_spfh, p_o3mr,                                           & ! met tracers
               p_du001, p_du002, p_du003, p_du004, p_du005,              & ! DU
               p_ss001, p_ss002, p_ss003, p_ss004, p_ss005,              & ! SS
               p_msa,   p_so4,   p_so2,   p_dms,                         & ! SU
               p_ocphobic, p_ocphilic,                                   & ! OC
               p_bcphobic, p_bcphilic                                      ! BC
           

! Fortran array for chem import state
      real(ESMF_KIND_R8), pointer, dimension(:,:) ::          &
               c_lwi,   c_zpbl, c_frlake,  c_fraci,  c_wet1, c_lai,     &
               c_grn,   c_cn_prcp, c_ncn_prcp, c_sh, c_ta,   c_tsoil1,  & 
               c_u10m,  c_v10m,  c_ustar,  c_z0h,  c_tropp, c_ps

      real (ESMF_KIND_R8), pointer, dimension(:,:,:) :: c_ple, c_zle,   &
               c_airdens, c_t, c_u, c_v, c_fcld, c_dqdt

      real (ESMF_KIND_R8), pointer, dimension(:,:,:) ::                  &
               c_o3, c_rh,                                               & ! met tracers
               c_du001, c_du002, c_du003, c_du004, c_du005,              & ! DU
               c_ss001, c_ss002, c_ss003, c_ss004, c_ss005,              & ! SS
               c_msa,   c_so4,   c_so2,   c_dms,                         & ! SU
               c_ocphobic, c_ocphilic,                                   & ! OC
               c_bcphobic, c_bcphilic                                      ! BC

! SOILTYPE and POROSITY
      integer, parameter :: DEFINED_SOIL=9
      real               :: SMCMAX, MAXSMC(DEFINED_SOIL)
!
! local variables for conversion
      real                    :: hs, ps, ptp, utp, vtp, ttp, htp, shrtp, &
                                 tv1, dz, tem
      real,  allocatable, dimension(:)  :: prsln, prslk, prsik, &
                                 pm, pd, sh, t, u, v, rh, shs, rho, pi, h


      real(kind=kind_phys), parameter :: rovg = con_rd / con_g
      real(kind=kind_phys), parameter :: qmin = 1.0e-10
      real(ESMF_KIND_R8),   parameter :: f_one = 1.0
      real(ESMF_KIND_R8),   parameter :: f_zero = 0.0

      data MAXSMC /0.421, 0.464, 0.468, 0.434, 0.406, 0.465,     &
                   0.404, 0.439, 0.421 /

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

!       allocate local working arrays
        allocate (                              &
                     pm   (km),                  &
                     pd   (km),                  &
                     sh   (km),                  &
                     t    (km),                  &
                     u    (km),                  &
                     v    (km),                  &
                     rh   (km),                  &
                     rho  (km),                  &
                     shs  (km),                  &
                     prslk(km),                  &
                     prsik(km+1),                &
                     prsln(km+1),                &
                     pi   (km+1),                &
                     h    (km+1)                 &
                    )

! debug print
        print *, 'LU_TST: cpi =', cpi(0:ntrac)
        print *, 'LU_TST: ri  =', ri(0:ntrac)
        print *, 'LU_TST: ntrac =', ntrac
        print *, 'LU_TST: doing_DU =', run_DU
        print *, 'LU_TST: doing_SU =', run_SU
        print *, 'LU_TST: doing_SS =', run_SS
        print *, 'LU_TST: doing_OC =', run_OC
        print *, 'LU_TST: doing_BC =', run_BC
	print *, 'LU_TST: lonr =', lonr
	print *, 'LU_TST: lats_node_r =', lats_node_r
	print *, 'LU_TST: lats_node_r_max =', lats_node_r_max
	print *, 'LU_TST: lonsperlar_r =', lonsperlar_r(:)
        print *, 'LU_TST: im, jm, km=', im, jm, km


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

       call GetPointer_(PHY_EXP_STATE,'slmsk', p_slmsk, rc)
       call GetPointer_(PHY_EXP_STATE,'hpbl',  p_hpbl , rc)
       call GetPointer_(PHY_EXP_STATE,'smc1',  p_smc1 , rc)
       call GetPointer_(PHY_EXP_STATE,'stype', p_stype, rc)
       call GetPointer_(PHY_EXP_STATE,'vtype', p_vtype, rc)
       call GetPointer_(PHY_EXP_STATE,'vfrac', p_vfrac, rc)
       call GetPointer_(PHY_EXP_STATE,'rain',  p_rain , rc)
       call GetPointer_(PHY_EXP_STATE,'rainc', p_rainc, rc)
       call GetPointer_(PHY_EXP_STATE,'dtsfci',p_dtsfci,rc)
       call GetPointer_(PHY_EXP_STATE,'tsea',  p_tsea , rc)
       call GetPointer_(PHY_EXP_STATE,'stc1',  p_stc1 , rc)
       call GetPointer_(PHY_EXP_STATE,'u10m',  p_u10m , rc)
       call GetPointer_(PHY_EXP_STATE,'v10m',  p_v10m , rc)
       call GetPointer_(PHY_EXP_STATE,'ustar', p_ustar, rc)
       call GetPointer_(PHY_EXP_STATE,'zorl',  p_zorl , rc)
       call GetPointer_(PHY_EXP_STATE,'hs'  ,  p_hs   , rc)
       call GetPointer_(PHY_EXP_STATE,'ps'  ,  p_ps   , rc)

! for GFS, vertical index is from surface to top of atmosphere
       call GetPointer_3D_(PHY_EXP_STATE,'t' ,  p_t   , rc)
       call GetPointer_3D_(PHY_EXP_STATE,'u' ,  p_u   , rc)
       call GetPointer_3D_(PHY_EXP_STATE,'v' ,  p_v   , rc)
       call GetPointer_3D_(PHY_EXP_STATE,'p' ,  p_p   , rc)
       call GetPointer_3D_(PHY_EXP_STATE,'dp',  p_dp  , rc)
       call GetPointer_3D_(PHY_EXP_STATE,'fcld', p_fcld , rc)
       call GetPointer_3D_(PHY_EXP_STATE,'dqdt', p_dqdt , rc)

! get met tracers
       call GetPointer_tracer_(PHY_EXP_STATE,'spfh', p_spfh, rc)
       call GetPointer_tracer_(PHY_EXP_STATE,'o3mr', p_o3mr, rc)

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
!* Get Fortran array from gocart import state
!---------------------------------------------

      MESSAGE_CHECK="Get ItemCount from chem import state"

      call esmf_stateget(CHEM_IMP_STATE                   &
                        ,itemcount = item_count_chem      &
                        ,itemnamelist = item_name         &
                        ,rc   =rc)

      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_CPL)

      print *,'LU_TST: CPL: chem_imp item count:',item_count_chem
      print *,'LU_TST: CPL: chem_imp item name :',      &
                       (item_name(i),i=1,item_count_chem)

      IF ( RC == ESMF_SUCCESS ) THEN

      if (item_count_chem == 0 ) then

       WRITE(0,*)'Empty chem import state; Fortran array not filled'

      else
       WRITE(0,*)'Get Fortran array from chem import state'
       call GetPointer_(CHEM_IMP_STATE,'lwi'   ,c_lwi  , rc)
       call GetPointer_(CHEM_IMP_STATE,'zpbl'  ,c_zpbl , rc)
       call GetPointer_(CHEM_IMP_STATE,'frlake',c_frlake,rc)
       call GetPointer_(CHEM_IMP_STATE,'fraci' ,c_fraci, rc)
       call GetPointer_(CHEM_IMP_STATE,'wet1'  ,c_wet1 , rc)
       call GetPointer_(CHEM_IMP_STATE,'lai'   ,c_lai  , rc)
       call GetPointer_(CHEM_IMP_STATE,'grn'   ,c_grn  , rc)
       call GetPointer_(CHEM_IMP_STATE,'cn_prcp',c_cn_prcp,rc)
       call GetPointer_(CHEM_IMP_STATE,'ncn_prcp',c_ncn_prcp,rc)
       call GetPointer_(CHEM_IMP_STATE,'sh'    ,c_sh   , rc)
       call GetPointer_(CHEM_IMP_STATE,'ta'    ,c_ta   , rc)
       call GetPointer_(CHEM_IMP_STATE,'tsoil1',c_tsoil1,rc)
       call GetPointer_(CHEM_IMP_STATE,'u10m'  ,c_u10m , rc)
       call GetPointer_(CHEM_IMP_STATE,'v10m'  ,c_v10m , rc)
       call GetPointer_(CHEM_IMP_STATE,'ustar' ,c_ustar, rc)
       call GetPointer_(CHEM_IMP_STATE,'z0h'   ,c_z0h  , rc)
       call GetPointer_(CHEM_IMP_STATE,'TROPP' ,c_tropp, rc)
       call GetPointer_(CHEM_IMP_STATE,'PS'    ,c_ps   , rc)

       call GetPointer_3D_(CHEM_IMP_STATE,'PLE', c_ple , rc)
       call GetPointer_3D_(CHEM_IMP_STATE,'ZLE', c_zle , rc)
       call GetPointer_3D_(CHEM_IMP_STATE,'AIRDENS',  c_airdens, rc)
       call GetPointer_3D_(CHEM_IMP_STATE,'T'  , c_t   , rc)
       call GetPointer_3D_(CHEM_IMP_STATE,'U'  , c_u   , rc)
       call GetPointer_3D_(CHEM_IMP_STATE,'V'  , c_v   , rc)
       call GetPointer_3D_(CHEM_IMP_STATE,'O3' , c_o3  , rc)
       call GetPointer_3D_(CHEM_IMP_STATE,'RH' , c_rh  , rc)
       call GetPointer_3D_(CHEM_IMP_STATE,'FCLD', c_fcld , rc)
       call GetPointer_3D_(CHEM_IMP_STATE,'DQDT', c_dqdt , rc)

      endif

      ELSE

       WRITE(0,*)'Fail to get chem import state'

      ENDIF


!---------------------------------------------
!* Do the actual coupling (data copy)
!---------------------------------------------
! 
      IF ( RC_CPL == ESMF_SUCCESS ) THEN

       if (item_count_chem == 0 ) then

        WRITE(0,*)'Empty chem imp; Coupling phy_exp with chem_imp skipped'

       else

        WRITE(0,*)'Chem imp not empty; Coupling phy_exp with chem_imp'
!
! ---  data copy between phy export state and chem import state
!
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
        c_lwi    = p_slmsk             ! land-ocean-ice mask  (1)
        c_ps     = p_ps                ! surface pressure (Pa)
        c_t      = p_t                 ! air temp at mid-layer (K)
        c_u      = p_u                 ! zonal wind at mid-layer (m/s)
        c_v      = p_v                 ! meridian wind at mid-layer (m/s)
        c_o3     = p_o3mr              ! ozone mixing ratio at mid-layer (kg/kg)
        c_fcld   = p_fcld              ! cloud cover  (1)
        c_dqdt   = p_dqdt              ! total moisture tendency (kg/kg/s)
!
!
! ---  derive chem import from phy export
!
!       do i = 1, im         ! <-- only loop through valid points
!       do j = 1, jm
        do j = 1, lats_node_r
        do i = 1, lonsperlar_r(j)

!         compute soil wetness
          if ( p_slmsk(i,j) == 1. ) then      ! <--- over ocean
            smcmax = maxsmc(p_vtype(i,j))     ! porosity
            c_wet1(i,j) = p_smc1(i,j)/smcmax
          else
            c_wet1(i,j) = 0.
          endif

!         setup local surface variables from dyn export state
          hs    = p_hs(i,j)                    ! surface height (m)
          ps    = p_ps(i,j)                    ! surface pressure (Pa)

!         setup local mid-layer variables from dyn export state
          pm(:) = p_p(i,j,:)                   ! pressure (Pa)
          pd(:) = p_dp(i,j,:)                  ! pressure thickness (Pa)
          sh(:) = max( p_spfh(i,j,:), qmin)    ! specific humidity (kg/kg)
          t(:)  = p_t(i,j,:)                   ! temp (K)
          u(:)  = p_u(i,j,:)                   ! zonal wind (m s-1)
          v(:)  = p_v(i,j,:)                   ! meridian wind (m s-1)

!         computes pi (interface pressure in Pa)
          pi(1)=ps
          do k=1,km-1
            pi(k+1) = pi(k) - pd(k)
          enddo
          pi(km+1) = f_zero

! ===>  rho and h can be calculated using real temp, following
!       get_r and get_phi routines (as in phy grid component) -- Sarah Lu

!         compute rho (layer air density in kg/m^3)
!         compute h (interface height in m)
!     
!      
          h(1) =  hs
          do k = 1, km                           ! from SFC to TOA
            tv1 = t(k) * (f_one + con_fvirt * sh(k))   ! virtual temp (k)
            rho(k) = pm(k) /(con_rd * tv1)       ! air density (kg/m3)

            prsln(k) = log(pi(k))
            if ( k == km ) prsln(k+1)=log(pm(k))
            dz = rovg * (prsln(k+1)-prsln(k)) * tv1  ! thickness (m)
            if ( k == km ) dz = 2.0 * dz
            h(k+1) = h(k) + dz        ! geopotential height at interface (m)
          enddo

!         call getrh to compute saturation humidity and relative humidity
!         rh : relative humidity (percent)
          call getrh(km,pm,sh,t,shs,rh)

!         call tpause to compute tropopause level fields
!         ptp: tropopause pressure (Pa)
          call tpause(km,pm,u,v,t,h,ptp,utp,vtp,ttp,htp,shrtp)

!         pass local variables to chem import state
          c_airdens(i,j,:) =  rho(:)   ! air density at mid-layer (kg/m3)
          c_tropp(i,j) = ptp           ! tropopause (Pa)
          c_ple(i,j,:) = pi(:)         ! air pressure at interface (Pa)
          c_zle(i,j,:) = h(:)          ! geopotential height at interface (m)
          c_rh(i,j,:)  = 0.01* rh(:)   ! relative humidity at mid-layer (0-1)

        enddo
        enddo

       endif

      ELSE

       WRITE(0,*)'Coupling phy_exp with chem_imp skipped, RC_CPL=',RC_CPL

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

        print *, 'LU_TST:',trim(statename),Name,Array(1,1,1),Array(1,2,1), &
                                          Array(im,jm,km)

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

!! 
!! adopt getrh routine from /nwprod/sorc/global_postgp.fd/postgp.f

      subroutine getrh(km,p,sh,t,shs,rh)
!
! Subprogram: getrh      Compute saturation humidity and relative humidity
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes the saturation specific humidity and the
!           relative humidity.  The relative humidity is constrained to be
!           between 0 and 100.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call getrh(km,p,sh,t,shs,rh)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     sh       real (km) specific humidity (kg/kg)
!     t        real (km) temperature (K)
!   Output argument list:
!     shs      real (km) saturation specific humidity (kg/kg)
!     rh       real (km) relative humidity (percent)
!
! Modules used:
!   funcphys       Physical functions
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   fpvs           compute saturation vapor pressure
!
! Attributes:
!   Language: Fortran 90
!
!$$$
       implicit none
       integer,intent(in):: km
       real,intent(in):: p(km),sh(km),t(km)
       real,intent(out):: shs(km),rh(km)
       real(krealfp) pr,tr,es
       integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       do k=1,km
         pr=p(k)
         tr=t(k)
         es=fpvs(tr)
         es=min(es,pr)
         shs(k)=con_eps*es/(pr+con_epsm1*es)
         rh(k)=1.e2*min(max(sh(k)/shs(k),0.),1.)
       enddo
       end subroutine


      END module atmos_phy_chem_cpl_comp_mod



