      module gfs_physics_grid_create_mod
!
!-------------------------------------------------------------------
! this code is used to create the esmf grids for the gfs esmf model.
! weiyu yang, 09/2005.
! updated by henry juang 04/2007
! updated by shrinivas moorthi for physics on 07/2007
! updated by henry juang 11/2007
! weiyu yang, 02/2008, updated to use the ESMF 3.1.0 library.
!-------------------------------------------------------------------
!
!!uses:
!
      use esmf_mod,                         ONLY: esmf_grid, esmf_vm,                &
                                                  ESMF_DistGrid, esmf_success,       &
                                                  ESMF_LogWrite, ESMF_LOG_INFO,      &
                                                  ESMF_DistGridCreate,               &
                                                  ESMF_LogMsgFoundError,             &
                                                  ESMF_FAILURE, ESMF_GridCreate
      use gfs_physics_internal_state_mod,   ONLY: gfs_physics_internal_state

      implicit none

      type(esmf_grid), save :: grid0   ! the esmf grid type array. for the 
                                       ! gfs start date and time information.
      type(esmf_grid), save :: grid3   ! the esmf grid type array.
                                       ! for the single gaussian grid arrays.
      type(esmf_grid), save :: grid4   ! the esmf grid type array.
                                       ! for the multiple gaussian grid arrays.


      contains

!------------------------------------------------------------------------
      subroutine gfs_physics_grid_create_gauss(vm, int_state,  &
                                               DistGrid0, DistGrid3, DistGrid4, rc)
!
! this routine create Gaussian grid type for single and multiple levels
! grid 3 (single) and grid4(multiple)
!
      type(esmf_vm),                     intent(inout) :: vm   
      type(gfs_physics_internal_state),  intent(inout) :: int_state
      integer,                           intent(out)   :: rc

      TYPE(ESMF_DistGrid),               INTENT(inout) :: DistGrid0    ! the ESMF DistGrid.
      TYPE(ESMF_DistGrid),               INTENT(inout) :: DistGrid3    ! the ESMF DistGrid.
      TYPE(ESMF_DistGrid),               INTENT(inout) :: DistGrid4    ! the ESMF DistGrid.

      integer                           :: rc1
      integer                           :: rcfinal

      integer,            dimension(2)  :: counts
      integer,            dimension(2)  :: arraystr, arrayend

      rc1     = esmf_success
      rcfinal = esmf_success

! create grid.
!=====================================================================
! set up parameter arrays for the esmf grid of the gaussian grid space.
! the first dimension is the multiple of latitudian and longitudian
! the second dimension is single level.
!----------------------------------------------------------------------
      counts(1)      = int_state%lonr*int_state%lats_node_r_max
      counts(1)      = counts(1) * int_state%nodes
      counts(2)      = 1
      arraystr(1)    = 1
      arraystr(2)    = 1
      arrayend(1)    = counts(1)
      arrayend(2)    = counts(2)

! Create the ESMF DistGrid3 using the 1-D default decomposition.
!---------------------------------------------------------------
      CALL ESMF_LogWrite("Create DistGrid3", ESMF_LOG_INFO, rc = rc1)

      DistGrid3 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create DistGrid3")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid3, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid3 based on the created ESMF DistGrid3 information.
! Grid3 is the grid for the Gaussian grid space.
!-----------------------------------------------------------------------
      CALL ESMF_LogWrite("create gfs_phy grid3", ESMF_LOG_INFO, rc = rc1)

      grid3 = ESMF_GridCreate(name = "gfs_phy grid3", distgrid = DistGrid3, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create Grid3")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid3, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

!=====================================================================
! set up parameter arrays for the esmf grid of the gaussian grid space.
! the first dimension is the multiple of latitudian and longitudian
! the second dimension is single level.
!--------------------------------------------------
      counts(1)      = int_state%lonr*int_state%lats_node_r_max
      counts(1)      = counts(1) * int_state%nodes
      counts(2)      = int_state%levs
      arraystr(1)    = 1
      arraystr(2)    = 1
      arrayend(1)    = counts(1)
      arrayend(2)    = counts(2)

! Create the ESMF DistGrid4 using the 1-D default decomposition.
!---------------------------------------------------------------
      CALL ESMF_LogWrite("Create DistGrid4", ESMF_LOG_INFO, rc = rc1)

      DistGrid4 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create DistGrid4")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid4, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid4 based on the created ESMF DistGrid4 information.
! Grid4 is the grid for the multiple level Gaussian grid space.
!-----------------------------------------------------------------------
      CALL ESMF_LogWrite("create gfs_phy grid4", ESMF_LOG_INFO, rc = rc1)

      grid4 = ESMF_GridCreate(name = "gfs_phy grid4", distgrid = DistGrid4, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create Grid4")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid4, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! set up parameter arrays for the esmf grid used for the date and time
! information to run the gfs.  all processors contains the same five date
! and time valus.
!------------------------------------------------------------------------
      counts(1)      = int_state%nodes
      counts(2)      = 5
      arraystr(1)    = 1
      arraystr(2)    = 1
      arrayend(1)    = counts(1)
      arrayend(2)    = counts(2)

! Create the ESMF DistGrid0 using the 1-D default decomposition.
!---------------------------------------------------------------
      CALL ESMF_LogWrite("Create DistGrid0", ESMF_LOG_INFO, rc = rc1)

      DistGrid0 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create DistGrid0")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid0, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! create the esmf grid for the date and time information.
!--------------------------------------------------------
      CALL ESMF_LogWrite("create create gfs_phy grid0", ESMF_LOG_INFO, rc = rc1)

      grid0 = ESMF_GridCreate(name = "gfs_phy grid0", distgrid = DistGrid0, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create Grid0")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid0, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

      if(rcfinal == esmf_success) then
          print*, "pass: gfs_physics_grid_create_gauss."
      else
          print*, "fail: gfs_physics_grid_create_gauss."
      end if

      rc = rcfinal

      end subroutine gfs_physics_grid_create_gauss

      end module gfs_physics_grid_create_mod
