    module GFS_Phy_States_Mod

!BOP

! !MODULE: GFS_Phy_States_Mod --- Define Physics Import/Export states

! !USES:
!
  use ESMF_Mod
  use machine,                         only: kind_evod
  use gfs_physics_err_msg_mod                         

  use gfs_physics_internal_state_mod    ! Physics internal state 
  use gfs_physics_namelist_mod          ! Physics configuration
  use gfs_physics_grid_create_mod, only: mgrid   

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:

  public gfs_physics_import2internal_mgrid,  &
         gfs_physics_internal2export_mgrid


! !REVISION HISTORY:
!  Sarah Lu  2009-08-04  First version
!  Sarah Lu  2009-10-12  Port to the latest trunk
!  Sarah Lu  2009-10-16  Tracer bundle added; (shum, oz, cld) removed
!  Sarah Lu  2009-11-13  2D diag fields added to export state (for GOCART)
!
!EOP

   contains

!BOP

! ========================================================================= 
! Update internal state (grid_fld) based on esmf import state                     
!
      subroutine gfs_physics_import2internal_mgrid(state, internal, rc)

! every possible import state has its own turn-on/turn-off switch flag
! which can be used to fit different interface requirement from different
! outside grid component systems.
!------------------------------------------------------------------------

      type(esmf_state),                           intent(in)    :: state 
      type(gfs_physics_internal_state), pointer,  intent(inout) :: internal
      integer, optional,                          intent(out)   :: rc     

      type(GFS_Phy_State_Namelist) :: cf

      integer                  :: rc1, rcfinal

      type(ESMF_Field)         :: Field
      type(ESMF_FieldBundle)   :: Bundle

      real, pointer     :: fArr2D(:,:)
      real, pointer     :: fArr3D(:,:,:)
!      real(ESMF_KIND_R8), pointer     :: fArr2D(:,:)
!      real(ESMF_KIND_R8), pointer     :: fArr3D(:,:,:)

      integer                  :: i, j, k

      cf = internal%esmf_sta_list

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

      call esmf_logwrite(						&
           " update internal state with the esmf import state", 	&
            esmf_log_info, rc = rc1)

! get the surface orography array from the esmf import state.
!------------------------------------------------------------
      if(cf%z_import == 1) then
          if(associated(fArr2D)) nullify(fArr2D)
          call ESMF_StateGet(state=State, ItemName='hs', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -z',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &
                             farray=fArr2D, rc = rc)
          if(internal%grid_aldata) then
           internal%grid_fld%z =  fArr2D                   
          else
           internal%grid_fld%z => fArr2D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -z',rcfinal)
      end if

! get the surface pressure array from the esmf import state.
!-----------------------------------------------------------
      if(cf%ps_import == 1) then
          if(associated(fArr2D)) nullify(fArr2D)
          call ESMF_StateGet(state=State, ItemName='ps', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -ps',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr2D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%ps =  fArr2D                   
          else
           internal%grid_fld%ps => fArr2D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -ps',rcfinal)
      end if

! get the temperature array from the esmf import state.
!------------------------------------------------------
      if(cf%temp_import == 1) then
          if(associated(fArr3D)) nullify(fArr3D)
          call ESMF_StateGet(state=State, ItemName='t', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -t',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%t =  fArr3D                   
          else
           internal%grid_fld%t => fArr3D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -t',rcfinal)
      end if

! get the zonal-wind array from the esmf import state.
!-----------------------------------------------------
      if(cf%u_import == 1) then
          if(associated(fArr3D)) nullify(fArr3D)
          call ESMF_StateGet(state=State, ItemName='u', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -u',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%u =  fArr3D                   
          else
           internal%grid_fld%u => fArr3D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -u',rcfinal)
      end if

! get the meridian-wind array from the esmf import state.
!-----------------------------------------------------
      if(cf%v_import == 1) then
          if(associated(fArr3D)) nullify(fArr3D)
          call ESMF_StateGet(state=State, ItemName='v', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -v',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%v =  fArr3D                   
          else
           internal%grid_fld%v => fArr3D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -v',rcfinal)
      end if

! get the tracer array from the esmf import state.                        
!-------------------------------------------------------------
      if(cf%tracer_import == 1) then
          call ESMF_StateGet(state=State, ItemName='tracers',    &   
                             fieldbundle=Bundle, rc = rc )          
          call gfs_physics_err_msg(rc,'retrieve bundle from state',rcfinal)
          do i = 1, internal%ntrac                                 
             if(associated(fArr3D)) nullify(fArr3D)           
             CALL ESMF_FieldBundleGet(bundle=Bundle, &                    
                                name=internal%gfs_phy_tracer%vname(i),&  
                                field=field, rc = rc)                
             call gfs_physics_err_msg(rc,'retrieve Efield from bundle',rcfinal)
             CALL ESMF_FieldGet(field=field, localDe=0, &                  
                                farray=fArr3D, rc = rc)                 
             call gfs_physics_err_msg(rc,'retrieve Farray from field',rcfinal)
             if(internal%grid_aldata) then                                 
	      internal%grid_fld%tracers(i)%flds = fArr3D
             else                                             
              internal%grid_fld%tracers(i)%flds => fArr3D 
             endif                                  
          end do                            
      end if                          


! get the pressure array from the esmf import state.
!-------------------------------------------------------------
      if(cf%p_import == 1) then
          if(associated(fArr3D)) nullify(fArr3D)
          call ESMF_StateGet(state=State, ItemName='p', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -p',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%p =  fArr3D                   
          else
           internal%grid_fld%p => fArr3D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -p',rcfinal)
      end if

! get the pressure layer depth (dp) array from the esmf import state.
!-------------------------------------------------------------
      if(cf%dp_import == 1) then
          if(associated(fArr3D)) nullify(fArr3D)
          call ESMF_StateGet(state=State, ItemName='dp', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -dp',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%dp =  fArr3D                   
          else
           internal%grid_fld%dp => fArr3D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -dp',rcfinal)
      end if

!
! get the omega (dpdt) array from the esmf import state.
!-------------------------------------------------------------
      if(cf%dpdt_import == 1) then
          if(associated(fArr3D)) nullify(fArr3D)
          call ESMF_StateGet(state=State, ItemName='dpdt', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -dpdt',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%dpdt =  fArr3D                   
          else
           internal%grid_fld%dpdt => fArr3D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -dpdt',rcfinal)
      end if
!
!
! print out the final error signal message and put it to rc.
!-----------------------------------------------------------
      call gfs_physics_err_msg_final(rcfinal,				&
                            "gfs_physics_import2internal_mgrid",rc)

      end subroutine gfs_physics_import2internal_mgrid

! ========================================================================= 
! Update esmf export state based on internal state (grid_fld)                          
!
      subroutine gfs_physics_internal2export_mgrid(internal, state, rc)

! every possible export state has its own turn-on/turn-off switch flag
! which can be used to fit different interface requirement from different
! outside grid component systems.
!------------------------------------------------------------------------

      type(esmf_state),                          intent(inout)  :: state 
      type(gfs_physics_internal_state), pointer, intent(in)     :: internal
      integer, optional,                         intent(out)    :: rc     

      type(GFS_Phy_State_Namelist) :: cf

      integer                  :: rc1, rcfinal

      type(ESMF_Field)         :: Field
      type(ESMF_FieldBundle)   :: Bundle

      real, pointer            :: fArr2D(:,:)
      real, pointer            :: fArr3D(:,:,:)
    
      integer                  :: i, j
!
!*  2-D/3-D quantities for GOCART gridded component
!*  nfld_2d and nfld_3d define the number of 2D and 3D fields
!*  the field names are specified in vname_2d and vname_3d
!
! TROPP: tropopause_pressure_based_on_blended_estimate --    
! LWI  : land-ocean-ice_mask             ==> slmsk, fice
! ZPBL : Planetary boundary layer height ==> hpbl   
! FRLAKE: fraction_of_lake  
! FRACI: ice_covered_fraction_of_tile   
! WET1 : surface_soil_wetness            ==> smc(1), stype
! LAI  : leaf_area_index                 ==> vtype  
! GRN  : greeness_fraction               ==> vfrac
! CN_PRCP: Surface Conv. rain flux needed by land ==> rainc
! NCN_PRCP: Non-convective precipitation ==> rain
! PS   : surface_pressure                ==>  ps   
! SH   : sensible_heat_flux_from_turbulence ==> dtsfci  
! TA   : surface_temperature_from_surface ==> tsea   
! TSOIL1: soil_temperatures_layer_1       ==> stc(1)  
! U10M : 10-meter_eastward_wind           ==> u10m   
! V10M : 10-meter_northward_wind          ==> v10m 
! USTAR: surface_velocity_scale           ==> uustar
! Z0H: surface_roughness_for_heat         ==> zorl 
!
      integer, parameter      :: nfld_2d = 16
      integer, parameter      :: nfld_3d = 2
      integer(ESMF_KIND_I4), allocatable, save :: lonsperlar_r(:)
      integer, save           :: lonr, lats_node_r, lats_node_r_max
      integer                 :: ilat
      character*50            :: msg
      character*8             :: vname, vname_2d(nfld_2d), &
                                 vname_3d(nfld_3d)

      data vname_2d /'slmsk', 'fice', 'hpbl', 'smc1',     &
                     'stype', 'vtype', 'vfrac', 'rainc', &
                     'rain', 'dtsfci', 'tsea', 'stc1',   &
                     'u10m', 'v10m',  'ustar','zorl'/

      data vname_3d /'cldcov','dqdt'/
 
      cf = internal%esmf_sta_list

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

      call esmf_logwrite(						&
           " point esmf export state to internal state", 	        &
            esmf_log_info, rc = rc1)

! put the surface orography array to esmf export state.
!------------------------------------------------------------
      if(cf%z_export == 1) then
       if(associated(fArr2D)) nullify(fArr2D)                              
       fArr2D => internal%grid_fld%z
       field = ESMF_FieldCreate(name='hs', grid=mgrid, fArray=fArr2D,&
               gridToFieldMap=(/1,2,0/),                              &
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -z",rcfinal)
      end if

! put the surface pressure array to the esmf export state.
!-----------------------------------------------------------
      if(cf%ps_export == 1) then
       if(associated(fArr2D)) nullify(fArr2D)                               
       fArr2D => internal%grid_fld%ps
       field = ESMF_FieldCreate(name='ps', grid=mgrid, fArray=fArr2D,&
               gridToFieldMap=(/1,2,0/),                              &
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -ps",rcfinal)
      end if

! put the temperature array to the esmf export state.
!------------------------------------------------------
      if(cf%temp_export == 1) then
       if(associated(fArr3D)) nullify(fArr3D)                              
       fArr3D => internal%grid_fld%t
       Field = ESMF_FieldCreate(name='t', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -t",rcfinal)
      end if

! put the zonal-wind array to the esmf export state.
!-----------------------------------------------------
      if(cf%u_export == 1) then
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_fld%u
       field = ESMF_FieldCreate(name='u', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -u",rcfinal)
      end if

! put the meridian-wind array to the esmf export state.
!-----------------------------------------------------
      if(cf%v_export == 1) then
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_fld%v
       field = ESMF_FieldCreate(name='v', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -v",rcfinal)
      end if

! put the tracer array to the esmf export state.                      
!-------------------------------------------------------------
      if(cf%tracer_export == 1) then
       bundle = ESMF_FieldBundleCreate(name='tracers', &       
                grid=mgrid, rc=rc)                            
       call gfs_physics_err_msg(rc,"create empty fieldbundle",rcfinal)   
       do i = 1, internal%ntrac                                     
          if(associated(fArr3D)) nullify(fArr3D)                  
          fArr3D => internal%grid_fld%tracers(i)%flds                
          field = ESMF_FieldCreate(name=internal%gfs_phy_tracer%vname(i),&
                  grid=mgrid, fArray=fArr3D, &                           
                  indexFlag=ESMF_INDEX_DELOCAL, rc=rc)                
          call ESMF_FieldBundleAdd(bundle,field,rc=rc)                
          call gfs_physics_err_msg(rc,"add Efield to bundle",rcfinal)    
       end do                                                         
       call ESMF_StateAdd(state,Bundle,rc=rc)                            
       call gfs_physics_err_msg(rc,"add to esmf state -tracer",rcfinal) 
      end if                                                       

! put the pressure array to the esmf export state.
!-------------------------------------------------------------
      if(cf%p_export == 1) then
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_fld%p
       field = ESMF_FieldCreate(name='p', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -p",rcfinal)
      end if

! put the pressure layer depth (dp) array to the esmf export state.
!-------------------------------------------------------------
      if(cf%dp_export == 1) then
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_fld%dp
       field = ESMF_FieldCreate(name='dp', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -dp",rcfinal)
      end if

!
! put the omega (dpdt) array to the esmf export state.
!-------------------------------------------------------------
      if(cf%dpdt_export == 1) then
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_fld%dpdt
       field = ESMF_FieldCreate(name='dpdt', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -dpdt",rcfinal)
      end if

!!****************************************************************************
!!****************************************************************************
!! add 2D diag fields (for GOCART) to esmf export state
!!------------------------------------------------------------

      lonr            = internal%lonr
      lats_node_r     = internal%lats_node_r
      lats_node_r_max = internal%lats_node_r_max

      if ( .not. allocated (lonsperlar_r)) then
        allocate ( lonsperlar_r(lats_node_r_max))
        lonsperlar_r(:) = 0
        do i = 1, lats_node_r         
          ilat = internal%global_lats_r(internal%ipt_lats_node_r-1+i) 
          lonsperlar_r(i) =  internal%lonsperlar(ilat)                        
        enddo
      endif

      print *, 'LU_TST: lonr            = ', lonr
      print *, 'LU_TST: lats_node_r     = ', lats_node_r
      print *, 'LU_TST: lats_node_r_max = ', lats_node_r_max
      print *, 'LU_TST: lonsperlar_r    = ', lonsperlar_r(:)

      CALL ESMF_AttributeSet(state=state                   &  !<-- The physics export state
                            ,name ='lats_node_r'           &  !<-- Name of the attribute to insert
                            ,value= lats_node_r            &  !<-- Value of the attribute
                            ,rc   =RC)

      CALL ESMF_AttributeSet(state=state                   &  
                            ,name ='lats_node_r_max'       & 
                            ,value= lats_node_r_max        & 
                            ,rc   =RC)

      CALL ESMF_AttributeSet(state=state                   &  
                            ,name ='lonr'                  & 
                            ,value= lonr                   &  
                            ,rc   =RC)

      CALL ESMF_AttributeSet(state=state                   &  !<-- The physics export state
                            ,name ='lonsperlar_r'          &  !<-- Name of the attribute to insert
                            ,count= lats_node_r_max        &  !<-- Number of values in the attribute
                            ,valueList =lonsperlar_r       &  !<-- Value of the attribute
                            ,rc   =RC)

! loop through the 2D diag fields 

      lab_do_2d : DO i = 1, nfld_2d 

        vname  = trim(vname_2d(i))
        if(associated(fArr2D)) nullify(fArr2D)                              

        SELECT CASE (vname)

!!        LWI: land-ocean-ice_mask
          CASE ('slmsk')   ! Land-sea mask (1=land; 0=sea)   
            fArr2D => internal%sfc_fld%slmsk

          CASE ('fice')    ! Ice concentration (ice>0; no ice=0)     
            fArr2D => internal%sfc_fld%fice

!!        Planetary boundary layer height (m)
          CASE ('hpbl')    ! Boundary layer height (m)  
            fArr2D => internal%flx_fld%hpbl

!!        surface_soil_wetness
          CASE ('smc1')      
            fArr2D => internal%sfc_fld%smc(1,:,:)
          CASE ('stype')      
            fArr2D => internal%sfc_fld%stype

!!        LAI: leaf_area_index
          CASE ('vtype')      
            fArr2D => internal%sfc_fld%vtype

!!        GRN: greeness_fraction
          CASE ('vfrac')      
            fArr2D => internal%sfc_fld%vfrac

!!        CN_PRCP: Surface Conv. rain flux (kg/m^2/s)
          CASE ('rainc')      
            fArr2D => internal%flx_fld%rainc

!!        NCN_PRCP: Non-convective precipitation rate (kg/m^2/s)
          CASE ('rain')      
            fArr2D => internal%flx_fld%rain

!!        SHFX: sensible_heat_flux_from_turbulence (W m-2)
          CASE ('dtsfci')      
            fArr2D => internal%flx_fld%dtsfci

!!        TA: Surface Air Temperature (K)
          CASE ('tsea')      
            fArr2D => internal%sfc_fld%tsea

!!        TSOIL1: soil_temperatures_layer_1 (k)
          CASE ('stc1')      
            fArr2D => internal%sfc_fld%stc(1,:,:)

!!        U10M: 10-meter_eastward_wind (m s-1)
          CASE ('u10m')      
            fArr2D => internal%flx_fld%u10m

!!        V10M: 10-meter_northward_wind (m s-1)
          CASE ('v10m')      
            fArr2D => internal%flx_fld%v10m

!!        USTAR: surface_velocity_scale (m s-1)
          CASE ('ustar')      
            fArr2D => internal%sfc_fld%uustar

!!        Z0H: surface_roughness_for_heat (m)
          CASE ('zorl')      
            fArr2D => internal%sfc_fld%zorl

        END SELECT 
       
        msg   = "Create ESMF Field from "//vname
        field = ESMF_FieldCreate(name=vname, grid=mgrid,       &
                fArray=fArr2D, gridToFieldMap=(/1,2,0/),       &
                indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
        call gfs_physics_err_msg(rc,msg,rcfinal)    

        msg   = "Add to Physics Export State"
        call ESMF_StateAdd(state, field, rc=rc)
        call gfs_physics_err_msg(rc, msg ,rcfinal)

      ENDDO   lab_do_2D

!
! print out the final error signal message and put it to rc.
!-----------------------------------------------------------
      call gfs_physics_err_msg_final(rcfinal,				&
                            "gfs_physics_internal2export_mgrid",rc)

      end subroutine gfs_physics_internal2export_mgrid

! =========================================================================

    end module GFS_Phy_States_Mod
