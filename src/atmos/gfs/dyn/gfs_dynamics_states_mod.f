       module gfs_dynamics_states_mod

! 
!  june 2005 		weiyu yang             initial code.
!  february 2007 	hann-ming henry juang  
!			for gfs dynamics and gaussian grid data.
!  March 2009           Weiyu Yang, modified for the ensemble NEMS run.
!  2009/10/05           Sarah Lu, grid_gr unfolded to 3D
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

      implicit none

      contains

! =======================================================================


      subroutine gfs_dynamics_import2internal(imp_gfs_dyn, int_state, rc)

! this subroutine can be used to update the initial condition 
! fields of the internal state from the esmf inport state.
!------------------------------------------------------------

! every possible import state has its own turn-on/turn-off switch flag
! which can be used to fit different interface requirement from different
! outside grid component systems.
!------------------------------------------------------------------------

      type(esmf_state)                                          :: imp_gfs_dyn  
      type(gfs_dynamics_internal_state), pointer, intent(inout) :: int_state 
      integer, optional,                          intent(out)   :: rc     

      integer                  :: rc1, rcfinal
      integer                  :: k, krq, kstr, kend
      INTEGER                  :: mstr, mend
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
      if( int_state%start_step ) then
        print *,' It is starting, so no need for import_state2internal '
        return
!     else
!       print *,' do import state to internal state '
!     endif

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

          int_state%grid_gr(:,:,kstr:kstr)=reshape(hold1,      &
             (/int_state%lonf,int_state%lats_node_a_max,1/))
          call gfs_dynamics_err_msg(rc1,"gete esmf state - hs_im",rcfinal)

      end if

! get the surface pressure array from the esmf import state.
! for the detailed comments for every computational steps
! please refer to the surface orography array code.
!-----------------------------------------------------------
      if(int_state%esmf_sta_list%ps_import == 1) then

          IF(int_state%ENS .AND. int_state%Cpl_flag == ESMF_TRUE) THEN
              mstr = int_state%g_q
              mend = int_state%g_q
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'pps', hold1, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mstr)=reshape(hold1,      &
                 (/int_state%lonf,int_state%lats_node_a_max,1/))
              mstr = int_state%g_qm
              mend = int_state%g_qm
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'psm', hold1, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mstr)=reshape(hold1,      &
                 (/int_state%lonf,int_state%lats_node_a_max,1/))
          ELSE
              kstr=int_state%g_zq
              kend=int_state%g_zq
              call getf90arrayfromstate(imp_gfs_dyn, 'ps', hold1, 0, rc = rc1)
              int_state%grid_gr(:,:,kstr:kstr)=reshape(hold1,      &
                 (/int_state%lonf,int_state%lats_node_a_max,1/))
          END IF

          call gfs_dynamics_err_msg(rc1,"gete esmf state - ps_im",rcfinal)

      end if

! get the temperature array from the esmf import state.
!------------------------------------------------------
      if(int_state%esmf_sta_list%temp_import == 1) then

          IF(int_state%ENS .AND. int_state%Cpl_flag == ESMF_TRUE) THEN
              mstr = int_state%g_tt
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'tt', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
              mstr = int_state%g_ttm
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'tm', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          ELSE
              call getf90arrayfromstate(imp_gfs_dyn, 't', hold2, 0, rc = rc1)
              kstr = int_state%g_t
              kend = kstr + int_state%levs - 1
              int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          END IF

          call gfs_dynamics_err_msg(rc1,"gete esmf state - t_im",rcfinal)

      end if

! get the zonal-wind array from the esmf import state.
! for detailed line by line comments please refer to 
! the temperature code.
!-----------------------------------------------------
      if(int_state%esmf_sta_list%u_import == 1) then

          IF(int_state%ENS .AND. int_state%Cpl_flag == ESMF_TRUE) THEN
              mstr = int_state%g_uu
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'uu', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
              mstr = int_state%g_uum
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'um', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          ELSE
              call getf90arrayfromstate(imp_gfs_dyn, 'u', hold2, 0, rc = rc1)
              kstr = int_state%g_u
              kend = kstr + int_state%levs - 1
              int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          END IF

          call gfs_dynamics_err_msg(rc1,"gete esmf state - u_im",rcfinal)


      end if

! get the meridian-wind array from the esmf import state.
!-----------------------------------------------------
      if(int_state%esmf_sta_list%v_import == 1) then

          IF(int_state%ENS .AND. int_state%Cpl_flag == ESMF_TRUE) THEN
              mstr = int_state%g_vv
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'vv', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
              mstr = int_state%g_vvm
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'vm', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          ELSE
              call getf90arrayfromstate(imp_gfs_dyn, 'v', hold2, 0, rc = rc1)
              kstr = int_state%g_v
              kend = kstr + int_state%levs - 1
              int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          END IF

          call gfs_dynamics_err_msg(rc1,"gete esmf state - v_im",rcfinal)

      end if

! get the moisture array from the esmf import state.
!---------------------------------------------------
      if(int_state%esmf_sta_list%q_import == 1) then

          IF(int_state%ENS .AND. int_state%Cpl_flag == ESMF_TRUE) THEN
              mstr = int_state%g_rq
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'sshum', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
              mstr = int_state%g_rm
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'shumm', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          ELSE
              kstr = int_state%g_rt
              kend = kstr + int_state%levs - 1
              call getf90arrayfromstate(imp_gfs_dyn, 'shum', hold2,0, rc = rc1)
              int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          END IF

          call gfs_dynamics_err_msg(rc1,"gete esmf state - shum_im",rcfinal)

      end if

! get the ozone array from the esmf import state.
!------------------------------------------------
      if(int_state%esmf_sta_list%oz_import == 1) then

          IF(int_state%ENS .AND. int_state%Cpl_flag == ESMF_TRUE) THEN
              mstr = int_state%g_rq + int_state%levs
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'ooz', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
              mstr = int_state%g_rm + int_state%levs
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'ozm', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          ELSE
              kstr = int_state%g_rt + int_state%levs
              kend = kstr + int_state%levs - 1
              call getf90arrayfromstate(imp_gfs_dyn, 'oz', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          END IF

          call gfs_dynamics_err_msg(rc1,"gete esmf state - oz_im",rcfinal)

      end if

! get the cloud liquid water array from the esmf import state.
!-------------------------------------------------------------
      if(int_state%esmf_sta_list%cld_import == 1) then

          IF(int_state%ENS .AND. int_state%Cpl_flag == ESMF_TRUE) THEN
              mstr = int_state%g_rq + int_state%levs * 2
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'ccld', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
              mstr = int_state%g_rm + int_state%levs * 2
              mend = mstr + int_state%levs - 1
              CALL GetF90ArrayFromState(imp_gfs_dyn, 'cldm', hold2, 0, rc = rc1)
              int_state%grid_gr(:,:,mstr:mend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          ELSE
              kstr = int_state%g_rt + int_state%levs * 2
              kend = kstr + int_state%levs - 1
              call getf90arrayfromstate(imp_gfs_dyn, 'cld', hold2, 0, rc = rc1)

              int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
                 (/int_state%lonf,int_state%lats_node_a_max,          &
                   int_state%levs/))
          END IF

          call gfs_dynamics_err_msg(rc1,"gete esmf state - cld_im",rcfinal)

      end if

! get the pressure array from the esmf import state.
!-------------------------------------------------------------
      if(int_state%esmf_sta_list%p_import == 1) then

          call getf90arrayfromstate(imp_gfs_dyn, 'p', hold2, 0, rc = rc1)
          kstr = int_state%g_p
          kend = kstr + int_state%levs - 1
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
          int_state%grid_gr(:,:,kstr:kend)=reshape(hold2,         &
             (/int_state%lonf,int_state%lats_node_a_max,          &
               int_state%levs/))

          call gfs_dynamics_err_msg(rc1,"gete esmf state - dpdt_im",rcfinal)

      end if

!
!
! print out the final error signal message and put it to rc.
!-----------------------------------------------------------
      IF(PRESENT(rc)) call gfs_dynamics_err_msg_final(rcfinal, &
          "gfs_dynamics_import2internal",rc)

      end subroutine gfs_dynamics_import2internal


! =========================================================================


      subroutine gfs_dynamics_internal2export(int_state, exp_gfs_dyn, rc)

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
      type(esmf_state),                           intent(inout) :: exp_gfs_dyn 
      type(gfs_dynamics_internal_state), pointer, intent(inout) :: int_state 
      integer, optional,                          intent(out)   :: rc     

! local array size parameter of the esmf export state arrays.
!------------------------------------------------------------
!     integer               :: llgg_s, lev_s

      integer               :: rc1     ! error signal work variable.
      integer               :: rcfinal ! the final error signal variable.
      integer               :: k, krq, kstr, kend
      INTEGER               :: mstr, mend
      real(kind=kind_evod), dimension(:,:), pointer :: hold0
      real(kind=kind_evod), dimension(:,:), pointer :: hold_z, hold_ps, hold_temp, hold_u,      &
                                                       hold_v, hold_q,  hold_oz,   hold_cld,    &
                                                       hold_p, hold_dp, hold_dpdt
      real(kind=kind_evod), dimension(:,:), pointer :: hold_pps, hold_ttemp, hold_uu,           &
                                                       hold_vv, hold_qq, hold_ooz, hold_ccld
      real(kind=kind_evod), dimension(:,:), pointer :: hold_psm, hold_tempm, hold_um,           &
                                                       hold_vm, hold_qm, hold_ozm, hold_cldm
      real(kind=kind_evod), dimension(:,:), pointer :: hold_pps6, hold_ttemp6, hold_uu6,        &
                                                       hold_vv6, hold_qq6, hold_ooz6, hold_ccld6
      real(kind=kind_evod), dimension(:,:), pointer :: hold_psm6, hold_tempm6, hold_um6,        &
                                                       hold_vm6, hold_qm6, hold_ozm6, hold_cldm6
      character(20) :: imp_item_name(100), exp_item_name(100)
      integer :: imp_item, exp_item, n
      integer :: ii1
      logical, save :: first
      data first/.true./
      save hold0,  hold_z,  hold_ps,   hold_temp, hold_u,   &
           hold_v, hold_q,  hold_oz,   hold_cld,            &
           hold_p, hold_dp, hold_dpdt,                      &
           hold_psm, hold_tempm, hold_um,                   &
           hold_vm, hold_qm, hold_ozm, hold_cldm,           &
           hold_pps, hold_ttemp, hold_uu,                   &
           hold_vv, hold_qq, hold_ooz, hold_ccld,           &
           hold_psm6, hold_tempm6, hold_um6,                &
           hold_vm6, hold_qm6, hold_ozm6, hold_cldm6,       &
           hold_pps6, hold_ttemp6, hold_uu6,                &
           hold_vv6, hold_qq6, hold_ooz6, hold_ccld6

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
          hold0     = 0.0
          hold_z    = 0.0
          hold_ps   = 0.0
          hold_temp = 0.0
          hold_u    = 0.0
          hold_v    = 0.0
          hold_q    = 0.0
          hold_oz   = 0.0
          hold_cld  = 0.0
          hold_p    = 0.0
          hold_dp   = 0.0
          hold_dpdt = 0.0

          IF(int_state%ENS) THEN
              ALLOCATE(hold_pps   (ii1, 1))
              ALLOCATE(hold_ttemp (ii1, int_state%levs))
              ALLOCATE(hold_uu    (ii1, int_state%levs))
              ALLOCATE(hold_vv    (ii1, int_state%levs))
              ALLOCATE(hold_qq    (ii1, int_state%levs))
              ALLOCATE(hold_ooz   (ii1, int_state%levs))
              ALLOCATE(hold_ccld  (ii1, int_state%levs))
              ALLOCATE(hold_psm   (ii1, 1))
              ALLOCATE(hold_tempm (ii1, int_state%levs))
              ALLOCATE(hold_um    (ii1, int_state%levs))
              ALLOCATE(hold_vm    (ii1, int_state%levs))
              ALLOCATE(hold_qm    (ii1, int_state%levs))
              ALLOCATE(hold_ozm   (ii1, int_state%levs))
              ALLOCATE(hold_cldm  (ii1, int_state%levs))
              ALLOCATE(hold_pps6  (ii1, 1))
              ALLOCATE(hold_ttemp6(ii1, int_state%levs))
              ALLOCATE(hold_uu6   (ii1, int_state%levs))
              ALLOCATE(hold_vv6   (ii1, int_state%levs))
              ALLOCATE(hold_qq6   (ii1, int_state%levs))
              ALLOCATE(hold_ooz6  (ii1, int_state%levs))
              ALLOCATE(hold_ccld6 (ii1, int_state%levs))
              ALLOCATE(hold_psm6  (ii1, 1))
              ALLOCATE(hold_tempm6(ii1, int_state%levs))
              ALLOCATE(hold_um6   (ii1, int_state%levs))
              ALLOCATE(hold_vm6   (ii1, int_state%levs))
              ALLOCATE(hold_qm6   (ii1, int_state%levs))
              ALLOCATE(hold_ozm6  (ii1, int_state%levs))
              ALLOCATE(hold_cldm6 (ii1, int_state%levs))
              hold_pps   = 0.0
              hold_ttemp = 0.0
              hold_uu    = 0.0
              hold_vv    = 0.0
              hold_qq    = 0.0
              hold_ooz   = 0.0
              hold_ccld  = 0.0
              hold_psm   = 0.0
              hold_tempm = 0.0
              hold_um    = 0.0
              hold_vm    = 0.0
              hold_qm    = 0.0
              hold_ozm   = 0.0
              hold_cldm  = 0.0
              hold_pps6   = 0.0
              hold_ttemp6 = 0.0
              hold_uu6    = 0.0
              hold_vv6    = 0.0
              hold_qq6    = 0.0
              hold_ooz6   = 0.0
              hold_ccld6  = 0.0
              hold_psm6   = 0.0
              hold_tempm6 = 0.0
              hold_um6    = 0.0
              hold_vm6    = 0.0
              hold_qm6    = 0.0
              hold_ozm6   = 0.0
              hold_cldm6  = 0.0
          END IF

      END IF
          
! print out the information.
!-----------------------------------------------
      print*, 'do int_state to exp_gfs_dyn'

      call esmf_logwrite("begining to put the esmf export state.", &
                    esmf_log_info, rc = rc1)

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

          IF(int_state%ENS .AND. (int_state%end_step .OR. first)) THEN
              mstr            = int_state%g_q
              hold_pps        = reshape(int_state%grid_gr(:,:,mstr),  &
                               (/int_state%lonf*int_state%lats_node_a_max,1/))
              hold_pps6(:, 1) = int_state%grid_gr6(:, mstr - 1)
              mstr            = int_state%g_qm
              hold_psm        = reshape(int_state%grid_gr(:,:,mstr),  &
                               (/int_state%lonf*int_state%lats_node_a_max,1/))
              hold_psm6(:, 1) = int_state%grid_gr6(:, mstr - 1)
          ELSE
              kstr=int_state%g_zq
              hold_ps = reshape(int_state%grid_gr(:,:,kstr),  &
                       (/int_state%lonf*int_state%lats_node_a_max,1/))
          END IF

          IF(first) THEN
              IF(int_state%ENS) THEN
                  kstr=int_state%g_zq
                  hold_ps = reshape(int_state%grid_gr(:,:,kstr),  &
                           (/int_state%lonf*int_state%lats_node_a_max,1/))
              END IF
              call addf90arraytostate(exp_gfs_dyn, grid3, 'ps',             & 
                                      hold_ps, rc = rc1)
              IF(int_state%ENS) THEN
                  call addf90arraytostate(exp_gfs_dyn, grid3, 'pps',        & 
                                          hold_pps, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid3, 'psm',        & 
                                          hold_psm, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid3, 'pps6',       & 
                                          hold_pps6, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid3, 'psm6',       & 
                                          hold_psm6, rc = rc1)
              END IF
              call gfs_dynamics_err_msg(rc1,"put to esmf state - ps_ex",rcfinal)
          END IF
      end if

! add the temperature fileds into the esmf export state.
!-------------------------------------------------------

      if(int_state%esmf_sta_list%temp_export == 1) then

          IF(int_state%ENS .AND. (int_state%end_step .OR. first)) THEN
              mstr        = int_state%g_tt
              mend        = mstr + int_state%levs - 1
              hold_ttemp  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                         (/int_state%lonf*int_state%lats_node_a_max,   &
                           int_state%levs/))
              hold_ttemp6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
              mstr        = int_state%g_ttm
              mend        = mstr + int_state%levs - 1
              hold_tempm  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                         (/int_state%lonf*int_state%lats_node_a_max,   &
                           int_state%levs/))
              hold_tempm6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
          ELSE
              kstr = int_state%g_t
              kend = kstr + int_state%levs - 1
              hold_temp = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                       (/int_state%lonf*int_state%lats_node_a_max,   &
                         int_state%levs/))
          END IF

          IF(first) THEN
              IF(int_state%ENS) THEN
                  kstr = int_state%g_t
                  kend = kstr + int_state%levs - 1
                  hold_temp = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                           (/int_state%lonf*int_state%lats_node_a_max,   &
                             int_state%levs/))
              END IF
              call addf90arraytostate(exp_gfs_dyn, grid4, 't',             & 
                                      hold_temp, rc = rc1)
              IF(int_state%ENS) THEN
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'tt',        & 
                                          hold_ttemp, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'tm',        & 
                                          hold_tempm, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'tt6',       & 
                                          hold_ttemp6, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'tm6',       & 
                                          hold_tempm6, rc = rc1)
              END IF
              call gfs_dynamics_err_msg(rc1,"put to esmf state - t_ex",rcfinal)
          END IF
      end if

! to add the u field into the esmf export state.  for the detailed
! description comments please refer to the temperature field.
!--------------------------------------------------------------------------

      if(int_state%esmf_sta_list%u_export == 1) then

          IF(int_state%ENS .AND. (int_state%end_step .OR. first)) THEN
              mstr     = int_state%g_uu
              mend     = mstr + int_state%levs - 1
              hold_uu  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                      (/int_state%lonf*int_state%lats_node_a_max,   &
                        int_state%levs/))
              hold_uu6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
              mstr     = int_state%g_uum
              mend     = mstr + int_state%levs - 1
              hold_um  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                      (/int_state%lonf*int_state%lats_node_a_max,   &
                        int_state%levs/))
              hold_um6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
          ELSE
              kstr = int_state%g_u
              kend = kstr + int_state%levs - 1
              hold_u  = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                     (/int_state%lonf*int_state%lats_node_a_max,   &
                       int_state%levs/))
          END IF

          IF(first) THEN
              IF(int_state%ENS) THEN
                  kstr = int_state%g_u
                  kend = kstr + int_state%levs - 1
                  hold_u  = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                         (/int_state%lonf*int_state%lats_node_a_max,   &
                           int_state%levs/))
              END IF
              call addf90arraytostate(exp_gfs_dyn, grid4, 'u',             & 
                                      hold_u, rc = rc1)
              IF(int_state%ENS) THEN
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'uu',        & 
                                          hold_uu, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'um',        & 
                                          hold_um, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'uu6',       & 
                                          hold_uu6, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'um6',       & 
                                          hold_um6, rc = rc1)
              END IF
              call gfs_dynamics_err_msg(rc1,"put to esmf state - u_ex",rcfinal)
          END IF

      end if

! add the v field into the esmf export state.
!----------------------------------------------------

      if(int_state%esmf_sta_list%v_export == 1) then

          IF(int_state%ENS .AND. (int_state%end_step .OR. first)) THEN
              mstr     = int_state%g_vv
              mend     = mstr + int_state%levs - 1
              hold_vv  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                      (/int_state%lonf*int_state%lats_node_a_max,   &
                        int_state%levs/))
              hold_vv6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
              mstr     = int_state%g_vvm
              mend     = mstr + int_state%levs - 1
              hold_vm  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                      (/int_state%lonf*int_state%lats_node_a_max,   &
                        int_state%levs/))
              hold_vm6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
          ELSE
              kstr = int_state%g_v
              kend = kstr + int_state%levs - 1
              hold_v  = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                     (/int_state%lonf*int_state%lats_node_a_max,   &
                       int_state%levs/))
          END IF

          IF(first) THEN
              IF(int_state%ENS) THEN
                  kstr = int_state%g_v
                  kend = kstr + int_state%levs - 1
                  hold_v  = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                         (/int_state%lonf*int_state%lats_node_a_max,   &
                           int_state%levs/))
              END IF
              call addf90arraytostate(exp_gfs_dyn, grid4, 'v',             & 
                                      hold_v, rc = rc1)
              IF(int_state%ENS) THEN
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'vv',        & 
                                          hold_vv, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'vm',        & 
                                          hold_vm, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'vv6',       & 
                                          hold_vv6, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'vm6',       & 
                                          hold_vm6, rc = rc1)
              END IF
              call gfs_dynamics_err_msg(rc1,"put to esmf state - v_ex",rcfinal)
          END IF

      end if

! add the moisture field into the esmf export state.
!---------------------------------------------------

      if(int_state%esmf_sta_list%q_export == 1) then

          IF(int_state%ENS .AND. (int_state%end_step .OR. first)) THEN
              mstr     = int_state%g_rq
              mend     = mstr + int_state%levs - 1
              hold_qq  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                      (/int_state%lonf*int_state%lats_node_a_max,   &
                        int_state%levs/))
              hold_qq6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
              mstr     = int_state%g_rm
              mend     = mstr + int_state%levs - 1
              hold_qm  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                      (/int_state%lonf*int_state%lats_node_a_max,   &
                        int_state%levs/))
              hold_qm6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
          ELSE
              kstr = int_state%g_rt
              kend = kstr + int_state%levs - 1
              hold_q  = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                     (/int_state%lonf*int_state%lats_node_a_max,   &
                       int_state%levs/))
          END IF

          IF(first) THEN
              IF(int_state%ENS) THEN
                  kstr = int_state%g_rt
                  kend = kstr + int_state%levs - 1
                  hold_q  = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                         (/int_state%lonf*int_state%lats_node_a_max,   &
                           int_state%levs/))
              END IF
              call addf90arraytostate(exp_gfs_dyn, grid4, 'shum',          & 
                                      hold_q, rc = rc1)
              IF(int_state%ENS) THEN
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'sshum',     & 
                                          hold_qq, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'shumm',     & 
                                          hold_qm, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'sshum6',    & 
                                          hold_qq6, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'shumm6',    & 
                                          hold_qm6, rc = rc1)
              END IF
              call gfs_dynamics_err_msg(rc1,"put to esmf state - shum_ex",rcfinal)
          END IF

      end if

! add the ozone field into the esmf export state.
!------------------------------------------------

      if(int_state%esmf_sta_list%oz_export == 1) then

          IF(int_state%ENS .AND. (int_state%end_step .OR. first)) THEN
              mstr      = int_state%g_rq + int_state%levs
              mend      = mstr + int_state%levs - 1
              hold_ooz  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                       (/int_state%lonf*int_state%lats_node_a_max,   &
                         int_state%levs/))
              hold_ooz6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
              mstr      = int_state%g_rm + int_state%levs
              mend      = mstr + int_state%levs - 1
              hold_ozm  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                       (/int_state%lonf*int_state%lats_node_a_max,   &
                         int_state%levs/))
              hold_ozm6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
          ELSE
              kstr = int_state%g_rt + int_state%levs
              kend = kstr + int_state%levs - 1
              hold_oz = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                     (/int_state%lonf*int_state%lats_node_a_max,   &
                       int_state%levs/))
          END IF

          IF(first) THEN
              IF(int_state%ENS) THEN
                  kstr = int_state%g_rt + int_state%levs
                  kend = kstr + int_state%levs - 1
                  hold_oz = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                         (/int_state%lonf*int_state%lats_node_a_max,   &
                           int_state%levs/))
              END IF
              call addf90arraytostate(exp_gfs_dyn, grid4, 'oz',          & 
                                      hold_oz, rc = rc1)
              IF(int_state%ENS) THEN
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'ooz',     & 
                                          hold_ooz, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'ozm',     & 
                                          hold_ozm, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'ooz6',    & 
                                          hold_ooz6, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'ozm6',    & 
                                          hold_ozm6, rc = rc1)
              END IF
              call gfs_dynamics_err_msg(rc1,"put to esmf state - oz_ex",rcfinal)
          END IF

      end if

! add the cloud liquid water field into the esmf export state.
!-------------------------------------------------------------

      if(int_state%esmf_sta_list%cld_export == 1) then

          IF(int_state%ENS .AND. (int_state%end_step .OR. first)) THEN
              mstr       = int_state%g_rq + int_state%levs * 2
              mend       = mstr + int_state%levs - 1
              hold_ccld  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                        (/int_state%lonf*int_state%lats_node_a_max,   &
                          int_state%levs/))
              hold_ccld6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
              mstr       = int_state%g_rm + int_state%levs * 2
              mend       = mstr + int_state%levs - 1
              hold_cldm  = reshape(int_state%grid_gr(:,:,mstr:mend),  &
                        (/int_state%lonf*int_state%lats_node_a_max,   &
                          int_state%levs/))
              hold_cldm6 = int_state%grid_gr6(:, mstr - 1 : mend - 1)
          ELSE
              kstr = int_state%g_rt + int_state%levs * 2
              kend = kstr + int_state%levs - 1
              hold_cld = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                      (/int_state%lonf*int_state%lats_node_a_max,   &
                        int_state%levs/))
          END IF

          IF(first) THEN
              IF(int_state%ENS) THEN
                  kstr = int_state%g_rt + int_state%levs * 2
                  kend = kstr + int_state%levs - 1
                  hold_cld = reshape(int_state%grid_gr(:,:,kstr:kend),  &
                          (/int_state%lonf*int_state%lats_node_a_max,   &
                            int_state%levs/))
              END IF
              call addf90arraytostate(exp_gfs_dyn, grid4, 'cld',          & 
                                      hold_cld, rc = rc1)
              IF(int_state%ENS) THEN
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'ccld',     & 
                                          hold_ccld, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'cldm',     & 
                                          hold_cldm, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'ccld6',    & 
                                          hold_ccld6, rc = rc1)
                  call addf90arraytostate(exp_gfs_dyn, grid4, 'cldm6',    & 
                                          hold_cldm6, rc = rc1)
              END IF
              call gfs_dynamics_err_msg(rc1,"put to esmf state - cld_ex",rcfinal)
          END IF

      end if

! add layer pressure (pp) fileds into the esmf export state.
!-------------------------------------------------------

      if(int_state%esmf_sta_list%p_export == 1) then

          kstr = int_state%g_p
          kend = kstr + int_state%levs - 1
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
                        ,rc   =rc1)
      call gfs_dynamics_err_msg(rc1,                                &
          "esmf_stateget in gfs_dynamics_internal2export ",rcfinal)
!      print *,' dynamics state, item count is ',imp_item
!      print *,' dynamics state, item name ',(imp_item_name(n),n=1,imp_item)

!
!! print out the final error signal information and put it to the rc.
!!-------------------------------------------------------------------
      IF(PRESENT(rc)) call gfs_dynamics_err_msg_final(rcfinal, &
          "gfs_dynamics_internal2export",rc)

      first = .false.

      end subroutine gfs_dynamics_internal2export


! ==========================================================================


      end module gfs_dynamics_states_mod
