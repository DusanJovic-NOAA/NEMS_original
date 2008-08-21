      module gfs_dynamics_grid_create_mod
!
!-------------------------------------------------------------------
! this code is used to create the esmf grids for the gfs esmf model.
! weiyu yang, 09/2005.
! updated by henry juang 04/2007
! weiyu yang, 02/2008, updated to use the ESMF 3.1.0 library.
!-------------------------------------------------------------------
!
!!uses:
!
      use esmf_mod                        ! the esmf library.
      use gfs_dynamics_internal_state_mod ! the contents of the esmf internal state.

      implicit none

      type(esmf_grid), save :: grid0   ! the esmf grid type array. for the 
                                       ! gfs start date and time information.
      type(esmf_grid), save :: grid1   ! the esmf grid type array.
                                       ! for the single level spectral arrays.
      type(esmf_grid), save :: grid2   ! the esmf grid type array.
                                       ! for the multiple level spectral arrays.
      type(esmf_grid), save :: grid3   ! the esmf grid type array.
                                       ! for the single gaussian grid arrays.
      type(esmf_grid), save :: grid4   ! the esmf grid type array.
                                       ! for the multiple gaussian grid arrays.


      contains


      subroutine gfs_dynamics_grid_create_spect(vm, int_state, &
                                                DistGrid0, DistGrid1, DistGrid2, rc)

!
! this routine create grid type of spectral grid, in single and multiple levels
! spectral grid types: grid1 (single) and grid2 (mutiple)
!
      implicit none

      type(esmf_vm),                     intent(inout) :: vm 
      type(gfs_dynamics_internal_state), intent(inout) :: int_state  
      integer,                           intent(out)   :: rc  

      TYPE(ESMF_DistGrid),               INTENT(inout) :: DistGrid0    ! the ESMF DistGrid.
      TYPE(ESMF_DistGrid),               INTENT(inout) :: DistGrid1    ! the ESMF DistGrid.
      TYPE(ESMF_DistGrid),               INTENT(inout) :: DistGrid2    ! the ESMF DistGrid.

      integer                 :: rc1          ! error signal work variable.
      integer                 :: rcfinal      ! the final error signal variable.

      integer,  dimension(2)  :: counts       ! parameter array to set up the 
                                              ! size of the 2-d esmf grid.
      integer,  dimension(2)  :: arraystr, arrayend     
                                              ! parameter arrays to set up the
                                              ! start number and the end number of
                                              ! the esmf grid in each dimension.

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success
      rcfinal = esmf_success

! create grid.
! use uniform grid to represent both the gaussian grid 
! and the spectral space grids, since no dx, dy is needed.
!---------------------------------------------------------

!===========================================================================
! create the single level spectral esmf grid.  the first dimension is the
! spectral coefficient, that is a 1-d array.  thus the second dimension
! size is one.  grid starts from 1 and end at the total number of the
! spectral coefficients.
!------------------------------------------------------------------------
      counts(1)        = (int_state%jcap+1)*(int_state%jcap+2)
      counts(2)        = 1
      arraystr(1)      = 1
      arraystr(2)      = 1
      arrayend(1)      = counts(1)
      arrayend(2)      = counts(2)

      print *,' gfs_dyn grid1 dimension ',counts

! Create the ESMF DistGrid1 using the 1-D default decomposition.
!---------------------------------------------------------------
      CALL ESMF_LogWrite("Create DistGrid1", ESMF_LOG_INFO, rc = rc1)

      DistGrid1 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create DistGrid1")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid1, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid1 based on the created ESMF DistGrid1 information.
!-----------------------------------------------------------------------
      CALL ESMF_LogWrite("create gfs_dyn grid1", ESMF_LOG_INFO, rc = rc1)
 
      grid1 = ESMF_GridCreate(name = "gfs_dyn grid1", distgrid = DistGrid1, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create Grid1")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid1, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

!===========================================================================
! create the multiple level spectral esmf grid.  the first dimension is the
! spectral coefficient, that is a 2-d array.  thus the second dimension
! size is levs.  grid starts from 1 and end at the total number of the
! spectral coefficients.
!------------------------------------------------------------------------
      counts(1)        = (int_state%jcap+1)*(int_state%jcap+2)
      counts(2)        = int_state%levs
      arraystr(1)      = 1
      arraystr(2)      = 1
      arrayend(1)      = counts(1)
      arrayend(2)      = counts(2)

      print *,' gfs_dyn grid2 dimension ',counts

! Create the ESMF DistGrid2 using the 1-D default decomposition.
!---------------------------------------------------------------
      CALL ESMF_LogWrite("Create DistGrid2", ESMF_LOG_INFO, rc = rc1)

      DistGrid2 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create DistGrid2")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid2, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid2 based on the created ESMF DistGrid2 information.
!-----------------------------------------------------------------------
      CALL ESMF_LogWrite("create gfs_dyn grid2", ESMF_LOG_INFO, rc = rc1)

      grid2 = ESMF_GridCreate(name = "gfs_dyn grid2", distgrid = DistGrid2, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create Grid2")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid2, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! set up parameter arrays for the esmf grid used for the date and time
! information to run the gfs.  all processors contains the same five date
! and time valus.
!------------------------------------------------------------------------
      counts(1)      = int_state%nodes
      counts(2)      = 5
      arraystr(1)      = 1
      arraystr(2)      = 1
      arrayend(1)    = counts(1)
      arrayend(2)    = counts(2)

      print *,' gfs_dyn grid0 dimension ',counts

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
      CALL ESMF_LogWrite("create gfs_dyn grid0", ESMF_LOG_INFO, rc = rc1)

      grid0 = ESMF_GridCreate(name = "gfs_dyn grid0", distgrid = DistGrid0, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create Grid0")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid0, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! finally print out the error signal information and put it to "rc".
!-------------------------------------------------------------------
      if(rcfinal == esmf_success) then
          print*, "pass: gfs_dynamics_grid_create_spect."
      else
          print*, "fail: gfs_dynamics_grid_create_spect."
      end if

      rc = rcfinal

      end subroutine gfs_dynamics_grid_create_spect






      subroutine gfs_dynamics_grid_create_gauss(vm, int_state, &
                                                DistGrid0, DistGrid3, DistGrid4, rc)
!
! this routine create Gaussian grid type for single and multiple levels
! grid 3 (single) and grid4(multiple)
!
      use esmf_mod     
      use gfs_dynamics_internal_state_mod

      type(esmf_vm),                     intent(inout) :: vm   
      type(gfs_dynamics_internal_state), intent(inout) :: int_state
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
      counts(1)      = int_state%lonf*int_state%lats_node_a_max
      counts(1)      = counts(1) * int_state%nodes
      counts(2)      = 1
      arraystr(1)    = 1
      arraystr(2)    = 1
      arrayend(1)    = counts(1)
      arrayend(2)    = counts(2)
                                                                                
      print *,' gfs_dyn grid3 dimension ',counts

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
      CALL ESMF_LogWrite("create gfs_dyn grid3", ESMF_LOG_INFO, rc = rc1)

      grid3 = ESMF_GridCreate(name = "gfs_dyn grid3", distgrid = DistGrid3, rc = rc1)

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
      counts(1)      = int_state%lonf*int_state%lats_node_a_max
      counts(1)      = counts(1) * int_state%nodes
      counts(2)      = int_state%levs
      arraystr(1)    = 1
      arraystr(2)    = 1
      arrayend(1)    = counts(1)
      arrayend(2)    = counts(2)

      print *,' gfs_dyn grid4 dimension ',counts

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
      CALL ESMF_LogWrite("create gfs_dyn grid4", ESMF_LOG_INFO, rc = rc1)

      grid4 = ESMF_GridCreate(name = "gfs_dyn grid4", distgrid = DistGrid4, rc = rc1)

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
      CALL ESMF_LogWrite("create gfs_dyn grid0", ESMF_LOG_INFO, rc = rc1)

      grid0 = ESMF_GridCreate(name = "gfs_dyn grid0", distgrid = DistGrid0, rc = rc1)

      IF(ESMF_LogMsgFoundError(rc1, "Create Grid0")) THEN
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating Grid0, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

      if(rcfinal == esmf_success) then
          print*, "pass: gfs_dynamics_grid_create_gauss."
      else
          print*, "fail: gfs_dynamics_grid_create_gauss."
      end if

      rc = rcfinal

      end subroutine gfs_dynamics_grid_create_gauss

      end module gfs_dynamics_grid_create_mod
