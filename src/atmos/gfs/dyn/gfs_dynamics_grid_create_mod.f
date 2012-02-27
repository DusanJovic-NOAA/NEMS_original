#include "../../../ESMFVersionDefine.h"

#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)
#undef ESMF_520r
#define ESMF_LogFoundError ESMF_LogMsgFoundError
#else
#define ESMF_520r
#endif

      module gfs_dynamics_grid_create_mod
!
!-------------------------------------------------------------------
! this code is used to create the esmf grids for the gfs esmf model.
! weiyu yang, 09/2005.
! updated by henry juang 04/2007
! weiyu yang, 02/2008, updated to use the ESMF 3.1.0 library.
! Sarah Lu, 10/2009, add gfs_dynamics_grid_create_Gauss3D routine
!                    that creates mGrid (3D gaussian grid)
! Weiyu Yang, 05/2011, modified for using the ESMF 5.2.0r_beta_snapshot_07.
!-------------------------------------------------------------------
!
!!uses:
!
#ifdef ESMF_520r
      USE esmf
#else
      USE esmf_mod
#endif
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
      type(ESMF_Grid), save :: mGrid   ! Mid-layer 3D gaussian grid: (im,jm,km)


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

! Create the ESMF DistGrid1 using the 1-D default decomposition.
!---------------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("Create DistGrid1", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("Create DistGrid1", ESMF_LOG_INFO, rc = rc1)
#endif

      DistGrid1 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create DistGrid1")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create DistGrid1")) THEN
#endif
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid1, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid1 based on the created ESMF DistGrid1 information.
!-----------------------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("create gfs_dyn grid1", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("create gfs_dyn grid1", ESMF_LOG_INFO, rc = rc1)
#endif
 
      grid1 = ESMF_GridCreate(name = "gfs_dyn grid1", distgrid = DistGrid1, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create Grid1")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create Grid1")) THEN
#endif
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

! Create the ESMF DistGrid2 using the 1-D default decomposition.
!---------------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("Create DistGrid2", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("Create DistGrid2", ESMF_LOG_INFO, rc = rc1)
#endif

      DistGrid2 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create DistGrid2")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create DistGrid2")) THEN
#endif
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid2, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid2 based on the created ESMF DistGrid2 information.
!-----------------------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("create gfs_dyn grid2", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("create gfs_dyn grid2", ESMF_LOG_INFO, rc = rc1)
#endif

      grid2 = ESMF_GridCreate(name = "gfs_dyn grid2", distgrid = DistGrid2, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create Grid2")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create Grid2")) THEN
#endif
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

! Create the ESMF DistGrid0 using the 1-D default decomposition.
!---------------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("Create DistGrid0", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("Create DistGrid0", ESMF_LOG_INFO, rc = rc1)
#endif

      DistGrid0 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create DistGrid0")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create DistGrid0")) THEN
#endif
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid0, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! create the esmf grid for the date and time information.
!--------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("create gfs_dyn grid0", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("create gfs_dyn grid0", ESMF_LOG_INFO, rc = rc1)
#endif

      grid0 = ESMF_GridCreate(name = "gfs_dyn grid0", distgrid = DistGrid0, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create Grid0")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create Grid0")) THEN
#endif
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
#ifdef ESMF_520r
      USE esmf
#else
      USE esmf_mod
#endif
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
                                                                                
! Create the ESMF DistGrid3 using the 1-D default decomposition.
!---------------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("Create DistGrid3", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("Create DistGrid3", ESMF_LOG_INFO, rc = rc1)
#endif

      DistGrid3 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create DistGrid3")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create DistGrid3")) THEN
#endif
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid3, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid3 based on the created ESMF DistGrid3 information.
! Grid3 is the grid for the Gaussian grid space.
!-----------------------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("create gfs_dyn grid3", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("create gfs_dyn grid3", ESMF_LOG_INFO, rc = rc1)
#endif

      grid3 = ESMF_GridCreate(name = "gfs_dyn grid3", distgrid = DistGrid3, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create Grid3")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create Grid3")) THEN
#endif
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

! Create the ESMF DistGrid4 using the 1-D default decomposition.
!---------------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("Create DistGrid4", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("Create DistGrid4", ESMF_LOG_INFO, rc = rc1)
#endif

      DistGrid4 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create DistGrid4")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create DistGrid4")) THEN
#endif
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid4, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! Create the ESMF grid4 based on the created ESMF DistGrid4 information.
! Grid4 is the grid for the multiple level Gaussian grid space.
!-----------------------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("create gfs_dyn grid4", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("create gfs_dyn grid4", ESMF_LOG_INFO, rc = rc1)
#endif

      grid4 = ESMF_GridCreate(name = "gfs_dyn grid4", distgrid = DistGrid4, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create Grid4")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create Grid4")) THEN
#endif
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
#ifdef ESMF_520r
      CALL ESMF_LogWrite("Create DistGrid0", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("Create DistGrid0", ESMF_LOG_INFO, rc = rc1)
#endif

      DistGrid0 = ESMF_DistGridCreate(arraystr, arrayend, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create DistGrid0")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create DistGrid0")) THEN
#endif
          rcfinal = ESMF_FAILURE
          PRINT*, 'Error Happened When Creating DistGrid0, rc = ', rc1
          rc1     = ESMF_SUCCESS
      END IF

! create the esmf grid for the date and time information.
!--------------------------------------------------------
#ifdef ESMF_520r
      CALL ESMF_LogWrite("create gfs_dyn grid0", ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite("create gfs_dyn grid0", ESMF_LOG_INFO, rc = rc1)
#endif

      grid0 = ESMF_GridCreate(name = "gfs_dyn grid0", distgrid = DistGrid0, rc = rc1)

#ifdef ESMF_520r
      IF(ESMF_LogFoundError(rc1, msg="Create Grid0")) THEN
#else
      IF(ESMF_LogFoundError(rc1,     "Create Grid0")) THEN
#endif
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


!..................................................................................
!BOP
!
!    !IROUTINE: gfs_dynamics_grid_create_Gauss3D
!
!    !INTERFACE:
!
      subroutine gfs_dynamics_grid_create_Gauss3D(vm, iState, distGrid, rc )

!     !ARGUMENTS:

      type(ESMF_VM),                     intent(in)  :: vm
      type(gfs_dynamics_internal_state), intent(in)  :: iState
      TYPE(ESMF_DistGrid),               intent(out) :: distGrid

      integer,                           intent(out)  :: rc

!     !DESCRIPTION:
!
!     This routine creates a 3D Gaussian grid on mid-layer vertical levels.
!     The actual distribution has already been determined and such information
!     is contained in the internal state *iState*.
!
!     Consistent with the current design, the actual grid is returned in
!     module-scoped variable *mGrid*.
!
!     !REVISION HISTORY:
!      da Silva  05Feb2009  Based of similar MAPL routine
!      Sarah Lu  28Jul2009  Remove MAPL exception handling; specify coord info
!      Sarah Lu  02Oct2009  Port from local branch to the trunk
!
!EOP
!                                      ---

      integer                 :: I1, IN, J1, JN
      integer                 :: rc1, i, j

!     Local space for coordinate information
!     --------------------------------------
      real(kind=kind_grid), pointer :: centerX(:,:)
      real(kind=kind_grid), pointer :: centerY(:,:)

! initialize the error signal variables.
!---------------------------------------
      rc1     = esmf_success

#ifdef ESMF_520r
      CALL ESMF_LogWrite('CreateGauss3D', ESMF_LOGMSG_INFO, rc = rc1)
#else
      CALL ESMF_LogWrite('CreateGauss3D', ESMF_LOG_INFO, rc = rc1)
#endif

!     Create grid with index-space information from internal state
!     ------------------------------------------------------------
#ifdef ESMF_520r
      mGrid = ESMF_GridCreate(                                &
                   countsPerDEDim1=(/iState%lonf/),           &
                   countsPerDEDim2=iState%lats_nodes_a_fix,   &
                   countsPerDEDim3=(/iState%levs/),           &
                   coordDep1 = (/1,2/),                       &
                   coordDep2 = (/1,2/),                       &
                   coordDep3 = (/3/),                         &
                   gridEdgeLWidth = (/0,0,0/),                &
                   gridEdgeUWidth = (/0,0,0/),                &
                   indexflag = ESMF_INDEX_DELOCAL,            &
                   name = 'mGrid',                            &
                   rc = rc1 )
#else
      mGrid = ESMF_GridCreateShapeTile ( name='mGrid',        &
                   countsPerDEDim1=(/iState%lonf/),           &
                   countsPerDEDim2=iState%lats_nodes_a_fix,   &
                   countsPerDEDim3=(/iState%levs/),           &
                   coordDep1 = (/1,2/),                       &
                   coordDep2 = (/1,2/),                       &
                   coordDep3 = (/3/),                         &
                   gridEdgeLWidth = (/0,0,0/),                &
                   gridEdgeUWidth = (/0,0,0/),                &
                   indexflag = ESMF_INDEX_DELOCAL,            &
                   rc = rc1 )
#endif

!  Add coordinate information
!  --------------------------
   call ESMF_GridAddCoord(mGrid, rc = rc1 )

!  Retrieve the coordinates so we can set them
!  -------------------------------------------
#ifdef ESMF_520r
   call ESMF_GridGetCoord (mGrid, 1, localDE=0,               &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                           FARRAYPTR=centerX, rc = rc1 )

   call ESMF_GridGetCoord (mGrid, 2, localDE=0,               &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                           FARRAYPTR=centerY, rc = rc1 )
#else
   call ESMF_GridGetCoord (mGrid, 0, 1,                       &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                                FPTR=centerX, rc = rc1 )

   call ESMF_GridGetCoord (mGrid, 0, 2,                       &
                           staggerloc=ESMF_STAGGERLOC_CENTER, &
                                FPTR=centerY, rc = rc1 )
#endif
   call GridGetInterior_ (mGrid,i1,in,j1,jn)
   print *, 'LU_DYN: mGrid, i1,in,j1,jn ', i1,in,j1,jn        !chlu_debug

   do i = 1, iState%lonf
   do j = 1, iState%lats_node_a
     centerX(i,j) = iState%xlon(i,j)
     centerY(i,j) = iState%xlat(i,j)
   enddo
   enddo

!  Make sure we've got it right
!  ----------------------------
   call ESMF_GridValidate(mGrid, rc=rc1 )

   return
   end subroutine gfs_dynamics_grid_create_Gauss3D

!..................................................................................

!  This routine came from MAPL...
   subroutine GridGetInterior_(GRID,I1,IN,J1,JN)
    type (ESMF_Grid), intent(IN) :: grid
    integer, intent(OUT)         :: I1, IN, J1, JN

! local vars
    integer                               :: status
    character(len=ESMF_MAXSTR)            :: IAm='MAPL_GridGetInterior'

    type (ESMF_DistGrid)                  :: distGrid
    type(ESMF_DELayout)                   :: LAYOUT
    integer,               allocatable    :: AL(:,:)
    integer,               allocatable    :: AU(:,:)
    integer                               :: nDEs
    integer                               :: deId
    integer                               :: gridRank
    integer                               :: deList(1)

    integer                               :: ierr, rc1


#ifdef ESMF_520r
    CALL ESMF_LogWrite("GridGetInterior_", ESMF_LOGMSG_INFO, rc=rc1)
#else
    CALL ESMF_LogWrite("GridGetInterior_", ESMF_LOG_INFO, rc=rc1)
#endif

    call ESMF_GridGet    (GRID, dimCount=gridRank, distGrid=distGrid, rc=rc1)
    call ESMF_DistGridGet(distGRID, delayout=layout, rc=STATUS)
    call ESMF_DELayoutGet(layout, deCount =nDEs, localDeList=deList, rc=rc1 )
    deId = deList(1)

    allocate (AL(gridRank,0:nDEs-1), stat = ierr )
    allocate (AU(gridRank,0:nDEs-1), stat = ierr )

#ifdef ESMF_520r
    call ESMF_DistGridGet(distgrid, &
         minIndexPDe=AL, maxIndexPDe=AU, rc=rc1 )
#else
    call ESMF_DistGridGet(distgrid, &
         minIndexPDimPDe=AL, maxIndexPDimPDe=AU, rc=rc1 )
#endif

    I1 = AL(1, deId)
    IN = AU(1, deId)
    J1 = AL(2, deId)
    JN = AU(2, deId)
    deallocate(AU, AL)

  end subroutine GridGetInterior_

      end module gfs_dynamics_grid_create_mod
