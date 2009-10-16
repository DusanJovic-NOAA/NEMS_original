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

      real, pointer            :: fArr2D(:,:)
      real, pointer            :: fArr3D(:,:,:)

      integer                  :: i, j                  

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

! get the moisture array from the esmf import state.
!---------------------------------------------------
      if(cf%q_import == 1) then
          if(associated(fArr3D)) nullify(fArr3D)
          call ESMF_StateGet(state=State, ItemName='shum', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -q',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%q =  fArr3D                   
          else
           internal%grid_fld%q => fArr3D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -q',rcfinal)
      end if

! get the ozone array from the esmf import state.
!------------------------------------------------
      if(cf%oz_import == 1) then
          if(associated(fArr3D)) nullify(fArr3D)
          call ESMF_StateGet(state=State, ItemName='oz', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -oz',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%oz =  fArr3D                   
          else
           internal%grid_fld%oz => fArr3D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -oz',rcfinal)
      end if

! get the cloud liquid water array from the esmf import state.
!-------------------------------------------------------------
      if(cf%cld_import == 1) then
          if(associated(fArr3D)) nullify(fArr3D)
          call ESMF_StateGet(state=State, ItemName='cld', &
                             field=Field, rc=rc)
          call gfs_physics_err_msg(rc,'retrieve Efield from state -cld',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          if(internal%grid_aldata) then
           internal%grid_fld%cld =  fArr3D                   
          else
           internal%grid_fld%cld => fArr3D                   
          endif
          call gfs_physics_err_msg(rc,'retrieve Farray from field -cld',rcfinal)
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

      real, pointer            :: fArr2D(:,:)
      real, pointer            :: fArr3D(:,:,:)


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

! put the moisture array to the esmf export state.
!---------------------------------------------------
      if(cf%q_export == 1) then
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_fld%q
       field = ESMF_FieldCreate(name='shum', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -q",rcfinal)
      end if

! put the ozone array to the esmf export state.
!------------------------------------------------
      if(cf%oz_export == 1) then
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_fld%oz
       field = ESMF_FieldCreate(name='oz', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -oz",rcfinal)
      end if

! put the cloud liquid water array to the esmf export state.
!-------------------------------------------------------------
      if(cf%cld_export == 1) then
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_fld%cld
       field = ESMF_FieldCreate(name='cld', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_physics_err_msg(rc,"add to esmf export state -cld",rcfinal)
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
!
!
! print out the final error signal message and put it to rc.
!-----------------------------------------------------------
      call gfs_physics_err_msg_final(rcfinal,				&
                            "gfs_physics_internal2export_mgrid",rc)

      end subroutine gfs_physics_internal2export_mgrid

! =========================================================================

    end module GFS_Phy_States_Mod
