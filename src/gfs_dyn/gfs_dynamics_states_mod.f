       module gfs_dynamics_states_mod

! 
!  june 2005 		weiyu yang             initial code.
!  february 2007 	hann-ming henry juang  
!			for gfs dynamics and gaussian grid data.
!  2009/10/05           Sarah Lu, grid_gr unfolded to 3D
!  nov 2009             Jun wang 
!                       set digital filter variables to export state
!
!!uses:
!
      use esmf_mod                 ! the esmf library.

! the derived type of the internal state.
!----------------------------------------
      use gfs_dynamics_internal_state_mod

! routines which can be used to add a fortran array to 
! an esmf state and to get a fortran array from an 
! esmf state
!-----------------------------------------------------
      use gfs_dynamics_grid_create_mod
      use gfs_dynamics_add_get_state_mod
      use gfs_dynamics_err_msg_mod

      use gfs_dyn_mpi_def
      use gfs_dyn_date_def

      implicit none

      contains

! =======================================================================


      subroutine gfs_dynamics_import2internal(gc_gfs_dyn, 		&
                              imp_gfs_dyn, int_state, rc)

! this subroutine can be used to update the initial condition 
! fields of the internal state from the esmf inport state.
!------------------------------------------------------------

! every possible import state has its own turn-on/turn-off switch flag
! which can be used to fit different interface requirement from different
! outside grid component systems.
!------------------------------------------------------------------------

      type(esmf_gridcomp),                        intent(inout) :: gc_gfs_dyn     
      type(esmf_state)                                          :: imp_gfs_dyn  
      type(gfs_dynamics_internal_state), pointer, intent(inout) :: int_state 
      integer,                                    intent(out)   :: rc     

      type(esmf_vm)            :: vm  
      integer                  :: rc1, rcfinal
!     integer                  :: llgg_s, lev_s
      integer                  :: k, krq, kstr, kend
      real, dimension(:,:), pointer :: hold0
      real, dimension(:,:), pointer :: hold1
      real, dimension(:,:), pointer :: hold2


! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

      call esmf_logwrite(						&
           " update internal state with the esmf import state", 	&
            esmf_log_info, rc = rc1)

! check the internal state, if it is ready to run in the start_step
! then we don't need to do import to internal state
!-----------------------------------------------------------------------
      if( int_state% start_step ) then
        print *,' It is starting, so no need for import_state2internal '
        return
      else
        print *,' do import state to internal state '
      endif

! create the gathering parameter arrays which will be used
! to gather the distributed esmf state local array to the global array.
!----------------------------------------------------------------------
!     llgg_s = int_state%llgg_s
!     lev_s = int_state%levs
!     print *,' llgg_s lev_s ',llgg_s,lev_s

! idate1_im:  (1) --- fhour (integer), (2) - (5) --- idate.
!-----------------------------------------------------------

      if(int_state%esmf_sta_list%idate1_import == 1) then

          call getf90arrayfromstate(imp_gfs_dyn,'date',hold0,0, rc = rc1)
          int_state%fhour_idate(:,:)=hold0
          call gfs_dynamics_err_msg(rc1,"done idate1_import.",rcfinal)

      end if

! get the surface orography array from the esmf import state.
!------------------------------------------------------------
      if(int_state%esmf_sta_list%z_import == 1) then

          kstr=int_state%g_gz
          kend=int_state%g_gz
          call getf90arrayfromstate(imp_gfs_dyn, 'hs', hold1, 0, rc = rc1)
!*        int_state%grid_gr(:,kstr:kend)=hold1
          int_state%grid_gr(:,:,kstr:kstr)=reshape(hold1,      &
             (/int_state%lonf,int_state%lats_node_a_max,1/)) 
          call gfs_dynamics_err_msg(rc1,"gete esmf state - hs_im",rcfinal)

      end if

! get the surface pressure array from the esmf import state.
! for the detailed comments for every computational steps
! please refer to the surface orography array code.
!-----------------------------------------------------------
      if(int_state%esmf_sta_list%ps_import == 1) then

          kstr=int_state%g_zq
          kend=int_state%g_zq
          call getf90arrayfromstate(imp_gfs_dyn, 'ps', hold1, 0, rc = rc1)
!*        int_state%grid_gr(:,kstr:kend)=hold1
          int_state%grid_gr(:,:,kstr:kstr)=reshape(hold1,      &
             (/int_state%lonf,int_state%lats_node_a_max,1/)) 

          call gfs_dynamics_err_msg(rc1,"gete esmf state - ps_im",rcfinal)

      end if

! get the temperature array from the esmf import state.
!------------------------------------------------------
      if(int_state%esmf_sta_list%temp_import == 1) then

          call getf90arrayfromstate(imp_gfs_dyn, 't', hold2, 0, rc = rc1)
          kstr = int_state%g_t
          kend = kstr + int_state%levs - 1
!*        int_state%grid_gr(:,kstr:kend)=hold2
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          & 
               int_state%levs/))                                   

          call gfs_dynamics_err_msg(rc1,"gete esmf state - t_im",rcfinal)

      end if

! get the zonal-wind array from the esmf import state.
! for detailed line by line comments please refer to 
! the temperature code.
!-----------------------------------------------------
      if(int_state%esmf_sta_list%u_import == 1) then

          call getf90arrayfromstate(imp_gfs_dyn, 'u', hold2, 0, rc = rc1)
          kstr = int_state%g_u
          kend = kstr + int_state%levs - 1
!*        int_state%grid_gr(:,kstr:kend)=hold2
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          & 
               int_state%levs/))                                   

          call gfs_dynamics_err_msg(rc1,"gete esmf state - u_im",rcfinal)


      end if

! get the meridian-wind array from the esmf import state.
!-----------------------------------------------------
      if(int_state%esmf_sta_list%v_import == 1) then

          call getf90arrayfromstate(imp_gfs_dyn, 'v', hold2, 0, rc = rc1)
          kstr = int_state%g_v
          kend = kstr + int_state%levs - 1
!*        int_state%grid_gr(:,kstr:kend)=hold2
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          & 
               int_state%levs/))                                   

          call gfs_dynamics_err_msg(rc1,"gete esmf state - v_im",rcfinal)

      end if

! get the moisture array from the esmf import state.
!---------------------------------------------------
      if(int_state%esmf_sta_list%q_import == 1) then

          krq=int_state%g_rt
          call getf90arrayfromstate(imp_gfs_dyn, 'shum', hold2,0, rc = rc1)
          kstr = krq
          kend = kstr + int_state%levs - 1
!*        int_state%grid_gr(:,kstr:kend)=hold2
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          & 
               int_state%levs/))                                   

          call gfs_dynamics_err_msg(rc1,"gete esmf state - shum_im",rcfinal)

      end if

! get the ozone array from the esmf import state.
!------------------------------------------------
      if(int_state%esmf_sta_list%oz_import == 1) then

          krq = krq + int_state%levs
          call getf90arrayfromstate(imp_gfs_dyn, 'oz', hold2, 0, rc = rc1)
          kstr = krq
          kend = kstr + int_state%levs - 1
!*        int_state%grid_gr(:,kstr:kend)=hold2
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          & 
               int_state%levs/))                                   

          call gfs_dynamics_err_msg(rc1,"gete esmf state - oz_im",rcfinal)

      end if

! get the cloud liquid water array from the esmf import state.
!-------------------------------------------------------------
      if(int_state%esmf_sta_list%cld_import == 1) then

          krq = krq + int_state%levs
          call getf90arrayfromstate(imp_gfs_dyn, 'cld', hold2, 0, rc = rc1)
          kstr = krq
          kend = kstr + int_state%levs - 1
!*        int_state%grid_gr(:,kstr:kend)=hold2
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          & 
               int_state%levs/))                                   

          call gfs_dynamics_err_msg(rc1,"gete esmf state - cld_im",rcfinal)

      end if

! get the pressure array from the esmf import state.
!-------------------------------------------------------------
      if(int_state%esmf_sta_list%p_import == 1) then

          call getf90arrayfromstate(imp_gfs_dyn, 'p', hold2, 0, rc = rc1)
          kstr = int_state%g_p
          kend = kstr + int_state%levs - 1
!*        int_state%grid_gr(:,kstr:kend)=hold2
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          & 
               int_state%levs/))                                   

          call gfs_dynamics_err_msg(rc1,"gete esmf state - p_im",rcfinal)

      end if

!

! get the pressure layer depth (dp) array from the esmf import state.
!-------------------------------------------------------------
      if(int_state%esmf_sta_list%dp_import == 1) then

          call getf90arrayfromstate(imp_gfs_dyn, 'dp', hold2, 0, rc = rc1)
          kstr = int_state%g_dp
          kend = kstr + int_state%levs - 1
!*        int_state%grid_gr(:,kstr:kend)=hold2
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          & 
               int_state%levs/))                                   

          call gfs_dynamics_err_msg(rc1,"gete esmf state - dp_im",rcfinal)

      end if

!
! get the omega (dpdt) array from the esmf import state.
!-------------------------------------------------------------
      if(int_state%esmf_sta_list%dpdt_import == 1) then

          call getf90arrayfromstate(imp_gfs_dyn, 'dpdt', hold2, 0, rc = rc1)
          kstr = int_state%g_dpdt
          kend = kstr + int_state%levs - 1
!*        int_state%grid_gr(:,kstr:kend)=hold2 
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          & 
               int_state%levs/))                                   

          call gfs_dynamics_err_msg(rc1,"gete esmf state - dpdt_im",rcfinal)

      end if

!
!
! print out the final error signal message and put it to rc.
!-----------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,				&
                            "gfs_dynamics_import2internal",rc)

      end subroutine gfs_dynamics_import2internal


! =========================================================================


      subroutine gfs_dynamics_internal2export(gc_gfs_dyn,               &
                              int_state, exp_gfs_dyn, rc)

! 
! this subroutine will be changed that all export
! esmf states will get data directly from the gfs internal structure data arrays
! 

!
!!uses:
!
!
! !input/output variables and parameters:
!----------------------------------------

      type(esmf_gridcomp),                        intent(inout) :: gc_gfs_dyn     
      type(esmf_state),                           intent(inout) :: exp_gfs_dyn 
      type(gfs_dynamics_internal_state), pointer, intent(inout) :: int_state 
      integer,                          intent(out)   :: rc       

! local array size parameter of the esmf export state arrays.
!------------------------------------------------------------
!     integer               :: llgg_s, lev_s

      type(esmf_vm)         :: vm     ! esmf virtual machine.

      integer               :: rc1     ! error signal work variable.
      integer               :: rcfinal ! the dinal error signal variable.
      integer               :: k, krq, kstr, kend
      real, dimension(:,:), pointer :: hold0
      real, dimension(:,:), pointer :: hold_z, hold_ps, hold_temp, hold_u,   &
                                       hold_v, hold_q,  hold_oz,   hold_cld, &
                                       hold_p, hold_dp, hold_dpdt
      character(20) :: imp_item_name(20), exp_item_name(20)
      integer :: imp_item, exp_item, n
      integer :: ii1
      logical, save :: first
      data first/.true./
      save hold0,  hold_z,  hold_ps,   hold_temp, hold_u,   &
           hold_v, hold_q,  hold_oz,   hold_cld,            &
           hold_p, hold_dp, hold_dpdt

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! Allocate all working arrays.
!-----------------------------
!! grid_gr unfolded from 2D to 3D (Sarah Lu)
      IF(first) THEN
	  ii1 =size(int_state%fhour_idate, dim=1)
          allocate(hold0(ii1, 5))
!*        ii1 =size(int_state%grid_gr, dim=1)
          ii1 =size(int_state%grid_gr, dim=1) *     &
               size(int_state%grid_gr, dim=2)
          allocate(hold_z(ii1, 1))
          allocate(hold_ps(ii1, 1))
          allocate(hold_temp(ii1, int_state%levs))
          allocate(hold_u(ii1, int_state%levs))
          allocate(hold_v(ii1, int_state%levs))
          allocate(hold_q(ii1, int_state%levs))
          allocate(hold_oz(ii1, int_state%levs))
          allocate(hold_cld(ii1, int_state%levs))
          allocate(hold_p(ii1, int_state%levs))
          allocate(hold_dp(ii1, int_state%levs))
          allocate(hold_dpdt(ii1, int_state%levs))
      END IF
          
! print out the information.
!-----------------------------------------------
      print*, 'do int_state to exp_gfs_dyn'

      call esmf_logwrite("begining to put the esmf export state.", &
                    esmf_log_info, rc = rc1)

! getting the global vm for the gathering data purpose
! and the grid3, which is for the gaussian grid arrays.
!------------------------------------------------------
      call esmf_gridcompget(gc_gfs_dyn, vm = vm, rc = rc1)

      call gfs_dynamics_err_msg(rc1,                                   &
              "get the global vm in int2exp",rcfinal)


! get the local array size, which are for the parallel esmf interface states,
! from the internal state.  lnt2_s is for the spectral arrays.  llgg_s
! is for the longitude and latitude of the gaussian grid arrays, respectively.
!----------------------------------------------------------------------------- 
!     llgg_s = int_state%llgg_s
!     lev_s = int_state%levs


! idate1_im and idate1_ex:  (1) --- fhour (integer), (2) - (5) --- idate.
!-------------------------------------------------------------------------

      if(int_state%esmf_sta_list%idate1_export == 1) then
          hold0 = int_state%fhour_idate(:,:)
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn,grid0,'date',hold0,rc=rc1)
              call gfs_dynamics_err_msg(rc1,"put to esmf state - date_ex",rcfinal)
          END IF
      end if

! orography field. gaussian grid
!----------------------------------------------------------------

      if(int_state%esmf_sta_list%z_export == 1) then

          kstr=int_state%g_gz
          kend=int_state%g_gz
!*        hold_z = int_state%grid_gr(:,kstr:kend)
          hold_z = reshape(int_state%grid_gr(:,:,kstr), &
                   (/int_state%lonf*int_state%lats_node_a_max,1/))
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid3, 'hs',             & 
                                      hold_z, rc = rc1)
              call gfs_dynamics_err_msg(rc1,"done z_export.",rcfinal)
          END IF
      end if

! surface pressure
!----------------------------------------------------------------

      if(int_state%esmf_sta_list%ps_export == 1) then

          kstr=int_state%g_zq
          kend=int_state%g_zq
!*        hold_ps = int_state%grid_gr(:,kstr:kend)
          hold_ps = reshape(int_state%grid_gr(:,:,kstr),  &
                   (/int_state%lonf*int_state%lats_node_a_max,1/))
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid3, 'ps',             & 
                                      hold_ps, rc = rc1)
    
              call gfs_dynamics_err_msg(rc1,"put to esmf state - ps_ex",rcfinal)
          END IF
      end if

! add the temperature fileds into the esmf export state.
!-------------------------------------------------------

      if(int_state%esmf_sta_list%temp_export == 1) then

          kstr = int_state%g_t
          kend = kstr + int_state%levs - 1
!*        hold_temp = int_state%grid_gr(:,kstr:kend)
          hold_temp = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                   (/int_state%lonf*int_state%lats_node_a_max,   &
                     int_state%levs/))                           

          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid4, 't',             & 
                                      hold_temp, rc = rc1)
              call gfs_dynamics_err_msg(rc1,"put to esmf state - t_ex",rcfinal)
          END IF
      end if

! to add the u field into the esmf export state.  for the detailed
! description comments please refer to the temperature field.
!--------------------------------------------------------------------------

      if(int_state%esmf_sta_list%u_export == 1) then

          kstr = int_state%g_u
          kend = kstr + int_state%levs - 1
!*        hold_u = int_state%grid_gr(:,kstr:kend)
          hold_u = reshape(int_state%grid_gr(:,:,kstr:kend),     &
                   (/int_state%lonf*int_state%lats_node_a_max,   &
                     int_state%levs/))                           
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid4, 'u',             & 
                                      hold_u, rc = rc1)

              call gfs_dynamics_err_msg(rc1,"put to esmf state - u_ex",rcfinal)
          END IF

      end if

! add the v field into the esmf export state.
!----------------------------------------------------

      if(int_state%esmf_sta_list%v_export == 1) then

          kstr = int_state%g_v
          kend = kstr + int_state%levs - 1
!*        hold_v = int_state%grid_gr(:,kstr:kend)
          hold_v = reshape(int_state%grid_gr(:,:,kstr:kend),     &
                   (/int_state%lonf*int_state%lats_node_a_max,   &
                     int_state%levs/))                           
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid4, 'v',             & 
                                      hold_v, rc = rc1)

              call gfs_dynamics_err_msg(rc1,"put to esmf state - v_ex",rcfinal)
          END IF

      end if

! add the moisture field into the esmf export state.
!---------------------------------------------------

      if(int_state%esmf_sta_list%q_export == 1) then

          krq=int_state%g_rt

          kstr = krq
          kend = kstr + int_state%levs - 1
!*        hold_q = int_state%grid_gr(:,kstr:kend)
          hold_q = reshape(int_state%grid_gr(:,:,kstr:kend),     &
                   (/int_state%lonf*int_state%lats_node_a_max,   &
                     int_state%levs/))                           
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid4, 'shum',          & 
                                      hold_q, rc = rc1)

               call gfs_dynamics_err_msg(rc1,"put to esmf state - shum_ex",rcfinal)
          END IF

      end if

! add the ozone field into the esmf export state.
!------------------------------------------------

      if(int_state%esmf_sta_list%oz_export == 1) then

          krq = krq + int_state%levs

          kstr = krq
          kend = kstr + int_state%levs - 1
!*        hold_oz = int_state%grid_gr(:,kstr:kend)
          hold_oz = reshape(int_state%grid_gr(:,:,kstr:kend),    &
                   (/int_state%lonf*int_state%lats_node_a_max,   &
                     int_state%levs/))                           
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid4, 'oz',          & 
                                      hold_oz, rc = rc1)

              call gfs_dynamics_err_msg(rc1,"put to esmf state - oz_ex",rcfinal)
          END IF

      end if

! add the cloud liquid water field into the esmf export state.
!-------------------------------------------------------------

      if(int_state%esmf_sta_list%cld_export == 1) then

          krq = krq + int_state%levs

          kstr = krq
          kend = kstr + int_state%levs - 1
!*        hold_cld = int_state%grid_gr(:,kstr:kend)
          hold_cld = reshape(int_state%grid_gr(:,:,kstr:kend),   &
                   (/int_state%lonf*int_state%lats_node_a_max,   &
                     int_state%levs/))                           
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid4, 'cld',          & 
                                      hold_cld, rc = rc1)

              call gfs_dynamics_err_msg(rc1,"put to esmf state - cld_ex",rcfinal)
          END IF

      end if

! add layer pressure (pp) fileds into the esmf export state.
!-------------------------------------------------------

      if(int_state%esmf_sta_list%p_export == 1) then

          kstr = int_state%g_p
          kend = kstr + int_state%levs - 1
!*        hold_p = int_state%grid_gr(:,kstr:kend)
          hold_p = reshape(int_state%grid_gr(:,:,kstr:kend),     &
                   (/int_state%lonf*int_state%lats_node_a_max,   &
                     int_state%levs/))                           
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid4, 'p',          & 
                                      hold_p, rc = rc1)

              call gfs_dynamics_err_msg(rc1,                                &
                      "put to esmf state - p_ex",rcfinal)
          END IF

      end if

! add pressure depth (dp) fileds into the esmf export state.
!-------------------------------------------------------

      if(int_state%esmf_sta_list%dp_export == 1) then

          kstr = int_state%g_dp
          kend = kstr + int_state%levs - 1
!*        hold_dp = int_state%grid_gr(:,kstr:kend)
          hold_dp = reshape(int_state%grid_gr(:,:,kstr:kend),    &
                   (/int_state%lonf*int_state%lats_node_a_max,   &
                     int_state%levs/))                           
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid4, 'dp',          & 
                                      hold_dp, rc = rc1)

              call gfs_dynamics_err_msg(rc1,                                &
                      "put to esmf state - dp_ex",rcfinal)
          END IF

      end if


! add omega (dpdt) fileds into the esmf export state.
!-------------------------------------------------------

      if(int_state%esmf_sta_list%dpdt_export == 1) then

          kstr = int_state%g_dpdt
          kend = kstr + int_state%levs - 1
!*        hold_dpdt = int_state%grid_gr(:,kstr:kend)
          hold_dpdt = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                   (/int_state%lonf*int_state%lats_node_a_max,   &
                     int_state%levs/))                           
          IF(first) THEN
              call addf90arraytostate(exp_gfs_dyn, grid4, 'dpdt',          & 
                                      hold_dpdt, rc = rc1)

              call gfs_dynamics_err_msg(rc1,                                &
                      "put to esmf state - dpdt_ex",rcfinal)
          END IF

      end if

      call esmf_stateget(exp_gfs_dyn                                      &
                        ,itemcount = imp_item                           &
                        ,itemnamelist = imp_item_name                   &
                        ,rc   =rc)
!      print *,' dynamics state, item count is ',imp_item
!      print *,' dynamics state, item name ',(imp_item_name(n),n=1,imp_item)

!
!! print out the final error signal information and put it to the rc.
!!-------------------------------------------------------------------
      call gfs_dynamics_err_msg_final(rcfinal,                         &
                    "gfs_dynamics_internal2export",rc)

      IF(first) first = .false.

      end subroutine gfs_dynamics_internal2export


! ==========================================================================


      end module gfs_dynamics_states_mod
