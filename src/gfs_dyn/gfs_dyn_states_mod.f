   module GFS_Dyn_States_Mod

!BOP

! !MODULE: GFS_Dyn_States_Mod --- Define Dynamics Import/Export states

! !USES:
!
  use ESMF_Mod
  use gfs_dyn_machine,              only: kind_evod
  use gfs_dynamics_err_msg_mod                     

  use gfs_dynamics_internal_state_mod   ! Dynamics internal state 
  use gfs_dynamics_namelist_mod         ! Dynamics configuration
  use gfs_dynamics_grid_create_mod, only: mgrid   

  implicit none

  private

! !PUBLIC MEMBER FUNCTIONS:

  public GFS_Dyn_StatesDefine,               &             
         gfs_dynamics_import2internal_mgrid, &            
         gfs_dynamics_internal2export_mgrid              


! !REVISION HISTORY:
!  da Silva  2009-01-22  First version
!  Sarah Lu  2009-01-26  Revised to include all 3d atmos fields
!  Sarah Lu  2009-08-06  Add gfs_dynamics_import2internal_mgrid and
!                        gfs_dynamics_internal2export_mgrid
!  Sarah Lu  2009-10-12  Port to the latest trunk
!  Sarah Lu  2009-10-17  Tracer bundle added; (shum, oz, cld) removed
!  Jun Wang  2009-11-09  set difital filter variables to export state
!                 
!EOP

   contains

!BOP

! !IROUTINE: GFS_Dyn_StatesDefine ---

! !INTERFACE:

      subroutine GFS_Dyn_StatesDefine ( import, export, internal, rc ) 

! !ARGUMENTS:

      type(GFS_Dynamics_Internal_State), intent(in)    :: Internal
      integer,                           intent(out)   :: rc           

      type(ESMF_STATE),                  intent(inout) :: Import     
      type(ESMF_STATE),                  intent(inout) :: Export
!
! !DESCRIPTION: Given previously created but otherwise empty Import/
!               Export ESMF states, fill those states with the contents
!               of a previously allocated (non-ESMF) Internal state.  
!
! !BOP

    type(GFS_Dyn_State_Namelist) :: cf
    integer                      :: rc1, rcfinal                  

    rc1     = esmf_success                                        
    rcfinal = esmf_success                                      

    cf = internal%esmf_sta_list

!    if ( .not. associated(Internal%grid_gr) ) &
!         call die_gfs ('Dynamics internal state not initialized')

!    print *, 'LU_DYN: call StateDefine for import state'        
    call StateDefine_ ( import, internal, mgrid,                       &
                        cf%z_import, cf%ps_import, cf%u_import,        &
                        cf%v_import, cf%temp_import, cf%tracer_import, &
                        cf%p_import, cf%dp_import, cf%dpdt_import, rc1)
    call gfs_dynamics_err_msg(rc1,"StateDefine for import state",rcfinal) 

!    print *, 'LU_DYN: call StateDefine for export state'        
    call StateDefine_ ( export, internal, mgrid,                       &
                        cf%z_export, cf%ps_export, cf%u_export,        &
                        cf%v_export, cf%temp_export, cf%tracer_export, &
                        cf%p_export, cf%dp_export, cf%dpdt_export, rc1)
    call gfs_dynamics_err_msg(rc1,"StateDefine for export state",rcfinal) 

! print out the final error signal message and put it to rc.          
!-----------------------------------------------------------          
    call gfs_dynamics_err_msg_final(rcfinal,                          &
                            "done gfs_dyn_statesdefine",rc)            

    end subroutine GFS_Dyn_StatesDefine

! ---
    subroutine StateDefine_ ( state, internal, mgrid,        &
                              z, ps, u, v, temp, tracer,     &
                              p, dp, dpdt, rc1) 

    type(ESMF_Grid), intent(in) :: mGrid

    integer,         intent(in) :: z, ps, u, v, temp, tracer, &
                                   p, dp, dpdt 

    type(GFS_Dynamics_Internal_State), intent(in)    :: Internal
    type(ESMF_STATE),                  intent(inout) :: State
    integer,                           intent(out)   :: rc1             
    
!                          ---

    type(ESMF_Field)       :: Field

    integer                :: rc, rcfinal
    real, pointer          :: fArr2D(:,:)
    real, pointer          :: fArr3D(:,:,:)

    integer                :: i, j

! initialize the error signal variables.                                 
!---------------------------------------                                 
      rc     = esmf_success                                              
      rcfinal = esmf_success                                             

!   Surface orography 
    if ( z == 1 ) then
       i = internal%g_gz
       if(associated(fArr2D)) nullify(fArr2D)        
       fArr2D => internal%grid_gr(:,:,i)
       field = ESMF_FieldCreate(name='hs', grid=mgrid, fArray=fArr2D,&
               gridToFieldMap=(/1,2,0/),                              &
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf state - z",rcfinal)      
    end if

!   Surface pressure
    if ( ps == 1 ) then
       i = internal%g_zq
       if(associated(fArr2D)) nullify(fArr2D)        
       fArr2D => internal%grid_gr(:,:,i)
       field = ESMF_FieldCreate(name='ps', grid=mgrid, fArray=fArr2D,&
               gridToFieldMap=(/1,2,0/),                              &
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf state - ps",rcfinal)    
    end if

!   Temperature
    if ( temp == 1 ) then
       i = internal%g_t
       j = internal%g_t + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)        
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='t', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf state - temp",rcfinal)  
    end if

!   Zonal-wind
    if ( u == 1 ) then
       i = internal%g_u
       j = internal%g_u + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)        
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='u', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf state - u",rcfinal)     
    end if

!   Meridian-wind
    if ( v == 1 ) then
       i = internal%g_v
       j = internal%g_v + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)        
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='v', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf state - v",rcfinal)     
    end if

!   Pressure
    if ( p == 1 ) then
       i = internal%g_p
       j = internal%g_p + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)        
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='p', grid=mgrid, fArray=fArr3D,& 
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf state - p",rcfinal)     
    end if

!   Pressure layer depth 
    if ( dp == 1 ) then
       i = internal%g_dp
       j = internal%g_dp + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)        
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='dp', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf state - dp",rcfinal)     
    end if

!   Omega
    if ( dpdt == 1 ) then
       i = internal%g_dpdt
       j = internal%g_dpdt + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)        
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='dpdt', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf state - dpdt",rcfinal)     
    end if


! print out the final error signal message and put it to rc1.                
!-----------------------------------------------------------                  
      call gfs_dynamics_err_msg_final(rcfinal,                          &     
                            "StateDefine_",rc1)                               

  end subroutine StateDefine_

! ========================================================================= 
! Update internal state (grid_gr) based on esmf import state                     
!
      subroutine gfs_dynamics_import2internal_mgrid(state, internal, rc,exp_gfs_dyn)

      type(esmf_state),                           intent(in)    :: state 
      type(gfs_dynamics_internal_state), pointer, intent(inout) :: internal
      integer, optional,                          intent(out)   :: rc     
      type(esmf_state),optional,                  intent(in)    :: exp_gfs_dyn

      type(GFS_Dyn_State_Namelist) :: cf

      integer                  :: rc1, rcfinal

      type(ESMF_Field)         :: Field
      type(ESMF_FieldBundle)   :: Bundle

      real, pointer            :: fArr2D(:,:)
      real, pointer            :: fArr3D(:,:,:)
      integer                  :: i, j, k 

      cf = internal%esmf_sta_list

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

      call esmf_logwrite(						&
           " update internal state with the esmf state", 	        &
            esmf_log_info, rc = rc1)

! get the surface orography array from the esmf import state.
!------------------------------------------------------------
      if(cf%z_import == 1) then
        i = internal%g_gz
        if(.not. present(exp_gfs_dyn)) then  
          call ESMF_StateGet(state=State, ItemName='hs',              & 
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from state -z',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &           
                             farray=fArr2D, rc = rc)             
          internal%grid_gr(:,:,i) = fArr2D                   
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -z',rcfinal)
        else
          call ESMF_StateGet(state=exp_gfs_dyn, ItemName='hs_dfi',   &
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from expstate -z',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &
                             farray=fArr2D, rc = rc)
          internal%grid_gr(:,:,i) = fArr2D
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -z',rcfinal)
        endif
      end if

! get the surface pressure array from the esmf import state.
!-----------------------------------------------------------
      if(cf%ps_import == 1) then
        i = internal%g_zq
        if(.not. present(exp_gfs_dyn)) then  
          call ESMF_StateGet(state=State, ItemName='ps',              & 
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from state -ps',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &           
                             farray=fArr2D, rc = rc)             
          internal%grid_gr(:,:,i) = fArr2D                   
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -ps',rcfinal)
        else
          call ESMF_StateGet(state=exp_gfs_dyn, ItemName='ps_dfi',   &
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from expstate -ps',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &
                             farray=fArr2D, rc = rc)
          internal%grid_gr(:,:,i) = fArr2D
!          print *,'from imp2int,ps=',maxval(internal%grid_gr(:,:,i)),   &
!            minval(internal%grid_gr(:,:,i))
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -ps',rcfinal)
        endif
      end if

! get the temperature array from the esmf import state.
!------------------------------------------------------
      if(cf%temp_import == 1) then
        i = internal%g_t
        j = i + internal%levs - 1
        if(.not. present(exp_gfs_dyn)) then  
          call ESMF_StateGet(state=State, ItemName='t',              & 
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from state -t',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          internal%grid_gr(:,:,i:j) = fArr3D                   
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -t',rcfinal)
        else
          call ESMF_StateGet(state=exp_gfs_dyn, ItemName='t_dfi',    &
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from expstate -t',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &
                             farray=fArr3D, rc = rc)
          internal%grid_gr(:,:,i:j) = fArr3D
!          print *,'from imp2int,t=',maxval(internal%grid_gr(:,:,i:j)),   &
!            minval(internal%grid_gr(:,:,i:j))
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -t',rcfinal)
        endif
      end if

! get the zonal-wind array from the esmf import state.
!-----------------------------------------------------
      if(cf%u_import == 1) then
        i = internal%g_u
        j = i + internal%levs - 1
        if(.not. present(exp_gfs_dyn)) then
          call ESMF_StateGet(state=State, ItemName='u',              & 
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from state -u',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          internal%grid_gr(:,:,i:j) = fArr3D                   
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -u',rcfinal)
        else
          call ESMF_StateGet(state=exp_gfs_dyn, ItemName='u_dfi',   &
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from expstate -u',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &
                             farray=fArr3D, rc = rc)
          internal%grid_gr(:,:,i:j) = fArr3D
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -u',rcfinal)
        endif
      end if

! get the meridian-wind array from the esmf import state.
!-----------------------------------------------------
      if(cf%v_import == 1) then
        i = internal%g_v
        j = i + internal%levs - 1
        if(.not. present(exp_gfs_dyn)) then
          call ESMF_StateGet(state=State, ItemName='v',              & 
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from state -v',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          internal%grid_gr(:,:,i:j) = fArr3D                   
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -v',rcfinal)
        else
          call ESMF_StateGet(state=exp_gfs_dyn, ItemName='v_dfi',   &
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from expstate -v',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &
                             farray=fArr3D, rc = rc)
          internal%grid_gr(:,:,i:j) = fArr3D
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -v',rcfinal)
        endif
      end if

! get the tracer array from the esmf import state.                           
!-------------------------------------------------------------           
      if(cf%tracer_import == 1) then
        if(.not. present(exp_gfs_dyn)) then
          call ESMF_StateGet(state=State, ItemName='tracers',    &      
                             fieldbundle=Bundle, rc = rc )                    
          call gfs_dynamics_err_msg(rc,'retrieve Ebundle from state',rcfinal) 
          do k = 1, internal%ntrac                   
             i = internal%g_rt+ (k-1) * internal%levs                     
             j = i + internal%levs - 1                                    
             if(associated(fArr3D)) nullify(fArr3D)                         
             CALL ESMF_FieldBundleGet(bundle=Bundle, &                       
                             name=internal%gfs_dyn_tracer%vname(k),&     
                             field=Field, rc = rc)                     
             call gfs_dynamics_err_msg(rc,'retrieve Efield from bundle',rcfinal)
             CALL ESMF_FieldGet(field=field, localDe=0, &                    
                                farray=fArr3D, rc = rc)                  
             call gfs_dynamics_err_msg(rc,'retrieve Farray from field',rcfinal)
             internal%grid_gr(:,:,i:j) = fArr3D                         
          end do                                                   
        else
          do k = 1, internal%ntrac
             i = internal%g_rt+ (k-1) * internal%levs
             j = i + internal%levs - 1
             if(associated(fArr3D)) nullify(fArr3D)
             CALL ESMF_StateGet(state=exp_gfs_dyn,                 &
                      ItemName=trim(internal%gfs_dyn_tracer%vname(k))//'_dfi',  &
                      field=Field, rc = rc )
             call gfs_dynamics_err_msg(rc,'retrieve Efield from bundle',rcfinal)
             CALL ESMF_FieldGet(field=field, localDe=0, &
                                farray=fArr3D, rc = rc)
             call gfs_dynamics_err_msg(rc,'retrieve Farray from field',rcfinal)
             internal%grid_gr(:,:,i:j) = fArr3D
!          print *,'from imp2int,tracert=',trim(internal%gfs_dyn_tracer%vname(k)), &
!            maxval(internal%grid_gr(:,:,i:j)),   &
!            minval(internal%grid_gr(:,:,i:j))
          end do
        endif 
      end if                                                

! get the pressure array from the esmf import state.
!-------------------------------------------------------------
      if(cf%p_import == 1) then
        i = internal%g_p
        j = i + internal%levs - 1
        if(.not. present(exp_gfs_dyn)) then
          call ESMF_StateGet(state=State, ItemName='p',            & 
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from state -p',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          internal%grid_gr(:,:,i:j) = fArr3D                   
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -p',rcfinal)
        else
          call ESMF_StateGet(state=exp_gfs_dyn, ItemName='p_dfi',            &
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from expstate -p',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &
                             farray=fArr3D, rc = rc)
          internal%grid_gr(:,:,i:j) = fArr3D
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -p',rcfinal)
        endif
      end if

! get the pressure layer depth (dp) array from the esmf import state.
!-------------------------------------------------------------
      if(cf%dp_import == 1) then
        i = internal%g_dp
        j = i + internal%levs - 1
        if(.not. present(exp_gfs_dyn)) then
          call ESMF_StateGet(state=State, ItemName='dp',            & 
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from state -dp',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          internal%grid_gr(:,:,i:j) = fArr3D                   
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -dp',rcfinal)
        else
          call ESMF_StateGet(state=exp_gfs_dyn, ItemName='dp_dfi',            &
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from expstate -dp',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &
                             farray=fArr3D, rc = rc)
          internal%grid_gr(:,:,i:j) = fArr3D
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -dp',rcfinal)
        endif
      end if

!
! get the omega (dpdt) array from the esmf import state.
!-------------------------------------------------------------
      if(cf%dpdt_import == 1) then
        i = internal%g_dpdt
        j = i + internal%levs - 1
        if(.not. present(exp_gfs_dyn)) then
        
          call ESMF_StateGet(state=State, ItemName='dpdt',            & 
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from state -dpdt',rcfinal)
          CALL ESMF_FieldGet(field=Field, localDe=0, &           
                             farray=fArr3D, rc = rc)             
          internal%grid_gr(:,:,i:j) = fArr3D                   
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -dpdt',rcfinal)
        else
          call ESMF_StateGet(state=exp_gfs_dyn, ItemName='dpdt_dfi',            &
                             field=Field, rc = rc )
          call gfs_dynamics_err_msg(rc,'retrieve Efield from expstate -dpdt',rcfinal)
          CALL ESMF_FieldGet(field=field, localDe=0, &
                             farray=fArr3D, rc = rc)
          internal%grid_gr(:,:,i:j) = fArr3D
          call gfs_dynamics_err_msg(rc,'retrieve Farray from field -dpdt',rcfinal)
        endif
      end if
!
!
! print out the final error signal message and put it to rc.
!-----------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,				&
                            'gfs_dynamics_import2internal_mgrid',rc)

      end subroutine gfs_dynamics_import2internal_mgrid

! ========================================================================= 
! Update esmf export state based on internal state (grid_gr)                          
!
      subroutine gfs_dynamics_internal2export_mgrid(internal, state, rc)

! every possible export state has its own turn-on/turn-off switch flag
! which can be used to fit different interface requirement from different
! outside grid component systems.
!------------------------------------------------------------------------

      type(esmf_state),                           intent(inout) :: state 
      type(gfs_dynamics_internal_state), pointer, intent(in)    :: internal
      integer, optional,                          intent(out)   :: rc     

      type(GFS_Dyn_State_Namelist) :: cf

      integer                      :: rc1, rcfinal

      type(ESMF_Field)             :: Field
      type(ESMF_FieldBundle)       :: Bundle

      real, pointer                :: fArr2D(:,:)
      real, pointer                :: fArr3D(:,:,:)

      integer                      :: i, j, k
!
!jw test:
       type(ESMF_Field)             :: Field1
       type(ESMF_GRID)              :: mgrid1
       integer :: gridrank,gridrank1
!jw test end
!
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
       i = internal%g_gz
       if(associated(fArr2D)) nullify(fArr2D)                              
       fArr2D => internal%grid_gr(:,:,i)
       field = ESMF_FieldCreate(name='hs', grid=mgrid, fArray=fArr2D, &
               gridToFieldMap=(/1,2,0/),                              &
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf export state -z",rcfinal)
!
       if(internal%ndfi>0 .and. cf%z_import==1) then
         if(associated(fArr2D)) nullify(fArr2D)                              
         i=1
         fArr2D => internal%grid_gr_dfi%hs(:,:,i)
         field = ESMF_FieldCreate(name='hs_dfi', grid=mgrid,          &
               fArray=fArr2D,gridToFieldMap=(/1,2,0/),                &
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
         call ESMF_StateAdd(state,field,rc=rc)
         call gfs_dynamics_err_msg(rc,"add to esmf export state -z-dfi",rcfinal)
       endif
      end if

! put the surface pressure array to the esmf export state.
!-----------------------------------------------------------
      if(cf%ps_export == 1) then
       i = internal%g_zq
       if(associated(fArr2D)) nullify(fArr2D)                               
       fArr2D => internal%grid_gr(:,:,i)
       field = ESMF_FieldCreate(name='ps', grid=mgrid, fArray=fArr2D,&
               gridToFieldMap=(/1,2,0/),                              &
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf export state -ps",rcfinal)
!      
       if(internal%ndfi>0 .and. cf%ps_import==1) then
         if(associated(fArr2D)) nullify(fArr2D)
         i=1
         fArr2D => internal%grid_gr_dfi%ps(:,:,i)
         field = ESMF_FieldCreate(name='ps_dfi', grid=mgrid,          &
               fArray=fArr2D,gridToFieldMap=(/1,2,0/),                &
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
         call ESMF_StateAdd(state,field,rc=rc)
         call gfs_dynamics_err_msg(rc,"add to esmf export state -ps",rcfinal)
       endif

      end if

! put the temperature array to the esmf export state.
!------------------------------------------------------
      if(cf%temp_export == 1) then
       i = internal%g_t
       j = internal%g_t + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)                              
       fArr3D => internal%grid_gr(:,:,i:j)
       Field = ESMF_FieldCreate(name='t', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf export state -t",rcfinal)
!      
       if(internal%ndfi>0 .and. cf%temp_import==1) then
         if(associated(fArr3D)) nullify(fArr3D)
         fArr3D => internal%grid_gr_dfi%t(:,:,:)
         field = ESMF_FieldCreate(name='t_dfi', grid=mgrid,          &
               fArray=fArr3D,indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
         call ESMF_StateAdd(state,field,rc=rc)
         call gfs_dynamics_err_msg(rc,"add to esmf export state -ps",rcfinal)
       endif

      end if

! put the zonal-wind array to the esmf export state.
!-----------------------------------------------------
      if(cf%u_export == 1) then
       i = internal%g_u
       j = internal%g_u + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='u', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf export state -u",rcfinal)
!      
       if(internal%ndfi>0 .and. cf%u_import==1) then
         if(associated(fArr3D)) nullify(fArr3D)
         fArr3D => internal%grid_gr_dfi%u(:,:,:)
         field = ESMF_FieldCreate(name='u_dfi', grid=mgrid,          &
               fArray=fArr3D,indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
         call ESMF_StateAdd(state,field,rc=rc)
         call gfs_dynamics_err_msg(rc,"add to esmf export state -ps",rcfinal)
       endif
!
      end if

! put the meridian-wind array to the esmf export state.
!-----------------------------------------------------
      if(cf%v_export == 1) then
       i = internal%g_v
       j = internal%g_v + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='v', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf export state -v",rcfinal)
!      
       if(internal%ndfi>0 .and. cf%v_import==1) then
         if(associated(fArr3D)) nullify(fArr3D)
         fArr3D => internal%grid_gr_dfi%v(:,:,:)
         field = ESMF_FieldCreate(name='v_dfi', grid=mgrid,          &
               fArray=fArr3D,indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
         call ESMF_StateAdd(state,field,rc=rc)
         call gfs_dynamics_err_msg(rc,"add to esmf export state -ps",rcfinal)
       endif
!
      end if

! put the tracer array to the esmf export state.                            
!-------------------------------------------------------------           
      if(cf%tracer_export == 1) then
       bundle = ESMF_FieldBundleCreate(name='tracers', grid=mgrid, rc=rc)
       call gfs_dynamics_err_msg(rc,"create empty fieldbundle",rcfinal)     
       do k = 1, internal%ntrac               
           i = internal%g_rt+ (k-1) * internal%levs           
           j = i + internal%levs - 1                        
           if(associated(fArr3D)) nullify(fArr3D)                           
           fArr3D => internal%grid_gr(:,:,i:j)                               
           field = ESMF_FieldCreate(name=internal%gfs_dyn_tracer%vname(k),&  
                   grid=mgrid, fArray=fArr3D, &                            
                   indexFlag=ESMF_INDEX_DELOCAL, rc=rc)                  
           call ESMF_FieldBundleAdd(bundle,field,rc=rc)                     
           call gfs_dynamics_err_msg(rc,"add field to bundle",rcfinal)       
       end do                                                                
       call ESMF_StateAdd(state,Bundle,rc=rc)                               
       call gfs_dynamics_err_msg(rc,"add to esmf state - tracer",rcfinal)   
!
       if(internal%ndfi>0 .and. cf%tracer_import==1) then
         do k = 1, internal%ntrac
           i = (k-1) * internal%levs+1
           j = i + internal%levs - 1    
           if(associated(fArr3D)) nullify(fArr3D)
           fArr3D => internal%grid_gr_dfi%tracer(:,:,i:j)
           field = ESMF_FieldCreate(name=trim(internal%gfs_dyn_tracer%vname(k))//'_dfi',&
                   grid=mgrid, fArray=fArr3D, &
                   indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
!           print *,'int2exp,k=',k,i,j,'name=',trim(internal%gfs_dyn_tracer%vname(k))//'_dfi', 'rc=',rc
           call ESMF_StateAdd(state,field,rc=rc)
!           print *,'add tracer field into state'
           call gfs_dynamics_err_msg(rc,"add to esmf export state -ps",rcfinal)
         enddo
       endif
!
      end if                                                               

! put pressure array to the esmf export state.
!-------------------------------------------------------------
      if(cf%p_export == 1) then
       i = internal%g_p
       j = internal%g_p + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='p', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf export state -p",rcfinal)
!
       if(internal%ndfi>0 .and. cf%p_import==1) then
         if(associated(fArr3D)) nullify(fArr3D)
         fArr3D => internal%grid_gr_dfi%p(:,:,:)
         field = ESMF_FieldCreate(name='p_dfi', grid=mgrid,          &
               fArray=fArr3D,indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
         call ESMF_StateAdd(state,field,rc=rc)
         call gfs_dynamics_err_msg(rc,"add to esmf export state -p",rcfinal)
       endif
!
      end if

! put the pressure layer depth (dp) array to the esmf export state.
!-------------------------------------------------------------
      if(cf%dp_export == 1) then
       i = internal%g_dp
       j = internal%g_dp + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='dp', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf export state -dp",rcfinal)
!
       if(internal%ndfi>0 .and. cf%dp_import==1) then
         if(associated(fArr3D)) nullify(fArr3D)
         fArr3D => internal%grid_gr_dfi%dp(:,:,:)
         field = ESMF_FieldCreate(name='dp_dfi', grid=mgrid,          &
               fArray=fArr3D,indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
         call ESMF_StateAdd(state,field,rc=rc)
         call gfs_dynamics_err_msg(rc,"add to esmf export state -dp",rcfinal)
       endif
!
      end if

!
! put the omega (dpdt) array to the esmf export state.
!-------------------------------------------------------------
      if(cf%dpdt_export == 1) then
       i = internal%g_dpdt
       j = internal%g_dpdt + internal%levs - 1
       if(associated(fArr3D)) nullify(fArr3D)                               
       fArr3D => internal%grid_gr(:,:,i:j)
       field = ESMF_FieldCreate(name='dpdt', grid=mgrid, fArray=fArr3D,&
               indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
       call ESMF_StateAdd(state,field,rc=rc)
       call gfs_dynamics_err_msg(rc,"add to esmf export state -dpdt",rcfinal)
!
       if(internal%ndfi>0 .and. cf%dpdt_import==1) then
         if(associated(fArr3D)) nullify(fArr3D)
         fArr3D => internal%grid_gr_dfi%dpdt(:,:,:)
         field = ESMF_FieldCreate(name='dpdt_dfi', grid=mgrid,          &
               fArray=fArr3D,indexFlag=ESMF_INDEX_DELOCAL, rc=rc)
         call ESMF_StateAdd(state,field,rc=rc)
         call gfs_dynamics_err_msg(rc,"add to esmf export state -dpdt",rcfinal)
       endif

      end if
!
!
! print out the final error signal message and put it to rc.
!-----------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,				&
                            "gfs_dynamics_internal2export_mgrid",rc)

      end subroutine gfs_dynamics_internal2export_mgrid
 

  end module GFS_Dyn_States_Mod
