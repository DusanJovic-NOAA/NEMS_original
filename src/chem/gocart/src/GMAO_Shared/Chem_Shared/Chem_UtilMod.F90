#include "MAPL_Generic.h"
#if 0
#define min(x,y) amin(real(x),real(y))
#define MIN(x,y) AMIN(real(x),real(y))
#define max(x,y) amax(real(x),real(y))
#define MAX(x,y) AMAX(real(x),real(y))
#endif

!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_UtilMod --- Assorted Utilities for fvChem
!
! !INTERFACE:
!

   module  Chem_UtilMod

! !USES:

   use ESMF_Mod
   use MAPL_Mod

   use MAPL_CommsMod
   use MAPL_CFIOMod
!  use HorzBinMod
   use MAPL_HorzTransformMod

   use Chem_Mod                  ! Chemistry Base Class
   use mod_diag                  ! fvGCM diagnostics
   use m_die
   use m_StrTemplate             ! string templates
   use m_chars, only: uppercase

   implicit NONE

!
! !PUBLIIC MEMBER FUNCTIONS:
!

   PRIVATE
   PUBLIC  Chem_UtilMPread          ! Array reader
   PUBLIC  Chem_UtilNegFiller       ! Fills negative values in a column
   PUBLIC  Chem_UtilTroppFixer      ! Repairs tropopause pressure bad values
   PUBLIC  Chem_UtilGetTimeInfo     ! Time info on file
   PUBLIC  Chem_UtilExtractIntegers ! Extract integers from a character-delimited string

   PUBLIC tick      ! GEOS-4 stub
   PUBLIC mcalday   ! GEOS-4 stub
   PUBLIC pmaxmin   ! functional
   PUBLIC zenith    ! GEOS-4 stub

!
! !DESCRIPTION:
!
!  This module implements assorted odds & ends for fvChem.
!
! !REVISION HISTORY:
!
!  29oct2003  da Silva  First crack.
!  16aug2005  da Silva  Introduced scatter from MAPL_CommsMod.
!
!EOP
!-------------------------------------------------------------------------

   interface pmaxmin
     module procedure pmaxmin2d
     module procedure pmaxmin3d
   end interface

CONTAINS

#if 0

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_UtilMPread --- Reads fields from file and distribute
!
! !INTERFACE:
!
   subroutine Chem_UtilMPread_g5 ( filen, varn, nymd, nhms, &
                                i1, i2, ig, im, j1, j2, jg, jm, km, &
                                var2d, var3d, cyclic, grid )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

  character(len=*), intent(in) :: filen   ! GFIO compatible file name
  character(len=*), intent(in) :: varn    ! variable name
  integer, intent(in)          :: nymd, nhms ! date/time

                                          ! Distributed grid info:
  integer,      intent(in)     :: i1, i2  !   local  zonal indices
  integer,      intent(in)     :: ig      !   zonal ghosting
  integer,      intent(in)     :: im      !   global zonal dimension
  integer,      intent(in)     :: j1, j2  !   local meridional indices
  integer,      intent(in)     :: jg      !   meridional ghosting
  integer,      intent(in)     :: jm      !   global zonal dimension
  integer,      intent(in)     :: km      !   vertical dimension

  logical, OPTIONAL, intent(in) :: cyclic ! whether time dimension is periodic
 
                                          ! ESMF Grid; this is required
                                          !  in GEOS-5 under ESMF
  type(ESMF_Grid), OPTIONAL, intent(in) :: grid 

                                                  


! !OUTPUT PARAMETERS:

  real, OPTIONAL, intent(out), target :: var2d(i1-ig:i2+ig,j1-jg:j2+jg)
  real, OPTIONAL, intent(out), target :: var3d(i1-ig:i2+ig,j1-jg:j2+jg,km)

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  15Aug2006  da Silva  It is now a simple wrap around CFIOReadArray().
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  Iam = 'Chem_UtilMPread_g5'
    type(ESMF_TIME)             :: time
    integer                     :: yy, mm, dd, h, m, s, rc, STATUS

    real, pointer               :: ptr2(:,:), ptr3(:,:,:)

    
!   Convert time to ESMF format
!   ---------------------------
    call parseIntTime_ ( nymd, yy, mm, dd )
    call parseIntTime_ ( nhms,  h,  m,  s )
    call ESMF_TimeSet(time, yy=yy, mm=mm, dd=dd,  h=h,  m=m, s=s, rc=status )
    if ( status /= 0 ) call die(Iam,'failed to convert time')

!   Read either 2D or 3D array
!   --------------------------
    if ( .not. present(grid) ) &
       call die ( Iam,'when running under the ESMF "grid" must be specified' )

    if ( present(var2d) ) then

         ptr2 => var2d
         call MAPL_CFIORead ( varn, filen, time, grid, ptr2, rc=STATUS,  &
                              verbose = .true., time_is_cyclic=cyclic,   &
                              time_interp = .true. )

    else if ( present(var3d) ) then

         ptr3 => var3d
         call MAPL_CFIORead ( varn, filen, time, grid, ptr3, rc=STATUS,  &
                              verbose = .true., time_is_cyclic=cyclic,   &
                              time_interp = .true. )

    else

       call die ( Iam,'either "var2d" or "var3d" must be specified' )

    end if

    if ( status /= 0 ) call die(Iam,'cannot read '//trim(varn))


CONTAINS
    subroutine parseIntTime_ ( hhmmss, hour, min, sec )      
      integer, intent(in)  :: hhmmss
      integer, intent(out) :: hour, min, sec 
      hour = hhmmss / 10000
      min  = mod(hhmmss,10000)/100
      sec  = mod(hhmmss,100)
    end subroutine parseIntTime_

  end subroutine Chem_UtilMPread_g5

#endif

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_UtilMPread_g4 --- Reads fields from file and distribute
!
! !INTERFACE:
!
   subroutine Chem_UtilMPread ( filen, varn, nymd, nhms,             &
                                i1, i2, ig, im, j1, j2, jg, jm, km,  &
                                var2d, var3d, cyclic, grid,          &
                                ForceBinning,                        &
                                maskString, gridMask, instanceNumber )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

  character(len=*), intent(in) :: filen   ! GFIO compatible file name
  character(len=*), intent(in) :: varn    ! variable name
  integer, intent(in)          :: nymd, nhms ! date/time

                                          ! Distributed grid info:
  integer,      intent(in)     :: i1, i2  !   local  zonal indices
  integer,      intent(in)     :: ig      !   zonal ghosting
  integer,      intent(in)     :: im      !   global zonal dimension
  integer,      intent(in)     :: j1, j2  !   local meridional indices
  integer,      intent(in)     :: jg      !   meridional ghosting
  integer,      intent(in)     :: jm      !   global zonal dimension
  integer,      intent(in)     :: km      !   vertical dimension

  logical, OPTIONAL, intent(in) :: cyclic ! whether time dimension is periodic
                                          ! ESMF Grid; this is required
                                          !  in GEOS-5 under ESMF

                                       
       
  logical, OPTIONAL, intent(in) :: ForceBinning ! Whether remapping uses
                                                !  binning or interpolation
                                                ! Starting with ARCTAS default
                                                ! is .TRUE. so binning always
                                                ! takes place.

  type(ESMF_Grid), OPTIONAL, intent(in) :: grid

  character(len=*), OPTIONAL, intent(in) :: maskString !Delimited string of integers
  real, OPTIONAL, intent(in) :: gridMask(i1:i2,j1:j2)  !Grid mask (NOTE: No ghosting) 
  integer, OPTIONAL, intent(in) :: instanceNumber      !Instantiation, for debugging.

! !OUTPUT PARAMETERS:

  real, OPTIONAL, intent(out)   :: var2d(i1-ig:i2+ig,j1-jg:j2+jg)
  real, OPTIONAL, intent(out)   :: var3d(i1-ig:i2+ig,j1-jg:j2+jg,km)

! !DESCRIPTION: 
!
! !REVISION HISTORY:
!
!  28oct2003 da Silva  First crack.
!  03ayg2004 da Silva  Uses GetVarT for time interpolation
!  18nov2004 da Silva  Added cyclic option for climatological files.
!  30may2005 da Silva  Introduced template expansing in file names.
!  14aug2006 da Silva  Swap longitudes of input arrays (always), and
!                      made it read only one layer at a time. Notice that
!                      the lon-swap is dangerous and it is a temporary
!                      device until re retire this routine in favor of
!                      CFIOArrayRead().
!  15Aug2006 da Silva  Renamed with _g4 as this is now obsolete.
!  21Aug2006 da Silva  Back to GFIO compatible version; renamed CFIO one g5;
!                      now the GFIO one does horizontal interp/binning
!                      and swap only if necessary.
!  29Feb2008 Nielsen   Masking
!
!EOP
!-------------------------------------------------------------------------

    character(len=*), parameter ::  myname = 'Chem_UtilMPread_g4'
    character(len=*), parameter ::  Iam=myname 
    logical :: tcyclic, verb

    integer  :: READ_ONLY=1, nokm=0
    integer :: fid, rc, ios, k
    character(len=257) :: fname, vname
    real :: const
    integer :: imf, jmf, kmf                  ! dimensions on file
    real,  pointer :: latf(:), lonf(:), levf(:)
    real,  pointer :: lat(:),  lon(:)

    real, pointer :: ptr2(:,:), ptr2f(:,:)

    logical :: doingMasking, ForceBinning_

    INTEGER, ALLOCATABLE :: regionNumbers(:),flag(:)
    INTEGER, ALLOCATABLE :: mask(:,:)

!   type(HorzBinTransform)     :: Trans
    type(MAPL_HorzTransform)   :: Trans

    vname = trim(varn) ! make a copy

!   Consistency check
!   -----------------
    if ( .not. ( present(var2d) .or. present(var3d) ) ) then
       call die ( myname, 'missing var2d or var3d' )
    else if ( present(var2d) .and. present(var3d) ) then
       call die ( myname, 'either var2d or var3d, but not both' )
    end if

    if ( present(ForceBinning) ) then
         ForceBinning_ = ForceBinning
    else
         ForceBinning_ = .true.
    endif

    if ( present(cyclic) ) then
        tcyclic = cyclic
    else
        tcyclic = .false. ! by default time dimension is not periodic
    end if

    if ( .not. present(grid) ) then
       call die ( myname, 'ESMF is a must when running under the ESMF' )
    end if

!   When masking, both the mask and the string
!   of integers (region numbers) must be present
!   --------------------------------------------
    IF ( (PRESENT(maskString) .AND. .NOT. PRESENT(gridMask)) .OR.   &
         (PRESENT(gridMask) .AND. .NOT. PRESENT(maskString))      ) &
       CALL die ( myname, ": Both gridMask and maskString must be specified." )

    IF(PRESENT(gridMask)) THEN
       IF(TRIM(maskString) == "-1") THEN
        doingMasking = .FALSE.
       ELSE
        doingMasking = .TRUE.
       END IF
    ELSE
     doingMasking = .FALSE.
    END IF

!   Masking initialization
!   ----------------------
    SetMask: IF(doingMasking) THEN

     k = 32
     ALLOCATE(regionNumbers(k),flag(k),mask(i1:i2,j1:j2),STAT=ios)
     IF ( ios /= 0 ) CALL die ( myname, ": Cannot allocate for masking.")

!   Obtain region numbers from delimited list of integers
!   -----------------------------------------------------
     regionNumbers(:) = 0
     CALL Chem_UtilExtractIntegers(maskString,k,regionNumbers,RC=ios)
     IF ( ios /= 0 ) CALL die ( myname, ": Unable to extract integers for regionNumbers.")

!   How many integers were found?
!   -----------------------------
     flag(:) = 1
     WHERE(regionNumbers(:) == 0) flag(:) = 0
     k = SUM(flag)
     DEALLOCATE(flag,STAT=ios)
     IF ( ios /= 0 ) CALL die ( myname, ": Cannot dallocate flag.")

!   Set local mask to 1 where gridMask matches each integer (within precision!) 
!   ---------------------------------------------------------------------------
     mask(i1:i2,j1:j2) = 0
     DO ios=1,k
      WHERE(regionNumbers(ios)-0.01 <= gridMask .AND. &
            gridMask <= regionNumbers(ios)+0.01) mask = 1
     END DO

    END IF SetMask

!   Expand templates
!   ----------------
    if ( index(filen,'%') .gt. 0 ) then
         call StrTemplate ( fname, filen, xid='unknown', &
                            nymd=nymd, nhms=nhms )
    else
         fname = filen
    end if

!   Trick: when filename is "/dev/null" it is assumed the user wants to
!   set the veriable to a constant
!   -------------------------------------------------------------------
    if ( fname(1:9) == '/dev/null' ) then    
         ios = -1
         k = index(fname,':')
         if ( k > 9 ) then
              read(fname(k+1:),*,iostat=ios) const
         end if
         if ( ios /= 0 ) const = 0.0

         IF ( PRESENT(var2d) ) THEN
          var2d = const
          IF(doingMasking .AND. const /= 0.0) THEN
           WHERE(mask == 0) var2d = 0.0
          END IF
         END IF

         IF ( PRESENT(var3d) ) THEN
          var3d = const
          IF(doingMasking .AND. const /= 0.0) THEN
           DO k=1,km
            WHERE(mask == 0) var3d(:,:,k) = 0.0
           END DO
          END IF
         END IF

         if ( MAPL_am_I_root() ) then
            print *, myname // ': input file is /dev/null'
            print *, myname // ':    setting variable ' // trim(vname) // &
                               ' to constant ', const
         end if
         return    ! exit here
    end if

!   Get the lat/lons from the input grid
!   ------------------------------------
    call GridGetLatLons_ ( grid, lon, lat )

!   Read file
!   ---------
    if ( MAPL_am_I_root() ) then

!      Allocate work space for scatter
!      -------------------------------
       allocate(ptr2(im,jm),stat=ios)
       if ( ios /= 0 ) call die ( myname, 'cannot allocate ptr2' )


       print *, myname // ': Reading ' // trim(vname) // ' from GFIO file ' &
                // trim(fname) // ' at ', nymd, nhms

!      Open file
!      ---------
       call GFIO_Open ( fname, READ_ONLY, fid, rc )
       if ( rc .ne. 0 ) then
          call die(myname,'cannot open GFIO file '//trim(fname))
       end if

!      Query file metadata
!      -------------------
       call queryFile_ ( fid, imf, jmf, kmf, lonf, latf, levf, vname )

!      If binning, define transform
!      ----------------------------
       if (im < imf .OR. jm < jmf .OR. ForceBinning_ ) then
          call MAPL_HorzTransformCreate (Trans, imf, jmf, im, jm, rc=rc)
          if ( rc /= 0 ) call die(myname,'cannot create transform',rc)
       end if

!      Allocate work space to reada data in
!      ------------------------------------
       allocate(ptr2f(imf,jmf),stat=ios)
       if ( ios /= 0 ) call die ( myname, 'cannot allocate ptr2f' )

    end if ! masterproc

!   2D Variables
!   ------------
    verb = .true.
    if ( present(var2d) ) then

!        Read global array and swap longitudes
!        -------------------------------------
         if ( MAPL_am_I_root() ) then

!           Read the file
!           -------------
            call GFIO_GetVarT1 ( fid, trim(vname), nymd, nhms, imf, jmf, &
                                 nokm, 1, ptr2f, rc, tcyclic, fid )
            if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(vname) )

!           Interpolate/bin if necessary
!           ----------------------------
            call Regrid_ ( ptr2f, imf, jmf, ptr2, im, jm, lonf, lon, verb )

         end if ! masterproc

!        Scatter the array
!        -----------------           
         call ArrayScatter ( var2d, ptr2, grid, rc=ios )
         if ( ios /= 0 ) call die ( myname, 'cannot scatter '//trim(vname) )

!        Apply mask when present
!        -----------------------
         IF(doingMasking) THEN
          WHERE(mask(i1:i2,j1:j2) == 0) var2d(i1:i2,j1:j2) = 0.00
         END IF

!   3D Variables
!   ------------
    else

      do k = 1, km

!        Read 1 level of global array and swap longitudes
!        ------------------------------------------------
         if ( MAPL_am_I_root() ) then
    
!            Read the file
!            -------------
             call GFIO_GetVarT1 ( fid, trim(vname), nymd, nhms, imf, jmf, k, &
                                  1, ptr2f, rc, tcyclic, fid )
             if ( rc .ne. 0 ) call die(myname,'cannot read '//trim(vname) )

!            Interpolate/bin if necessary
!            ----------------------------
             call Regrid_ ( ptr2f, imf, jmf, ptr2, im, jm, lonf, lon, verb ) 

         end if ! masterproc

!        Scatter the array with 1 level
!        ------------------------------           
         call ArrayScatter ( var3d(:,:,k), ptr2, grid, rc=ios )
         if ( ios /= 0 ) call die ( myname, 'cannot scatter v3d'//trim(vname))

         verb = .false. ! true only for k=1

!        Apply mask when present
!        -----------------------
         IF(doingMasking) THEN
          WHERE(mask(i1:i2,j1:j2) == 0) var3d(i1:i2,j1:j2,k) = 0.00
         END IF

        end do ! over k

      end if ! 2D or 3D

!   Close file
!   ----------
    if ( MAPL_am_I_root() ) then

!      If binning, destroy transform
!      -----------------------------
       if (im < imf .OR. jm < jmf .or. ForceBinning_ ) then
          call MAPL_HorzTransformDestroy (Trans )
       end if

       call GFIO_Close ( fid, rc )
!       print *, myname // ': Closing GFIO file ' // trim(fname)

       deallocate(ptr2, ptr2f, lonf, latf, levf, stat=ios)
       if ( ios /= 0 ) call die ( myname, 'cannot deallocate ptr2/ptr2f' )

    end if ! masterproc

!   All done
!   --------
    deallocate(lat, lon, stat=ios)
    if ( ios /= 0 ) call die ( myname, 'cannot deallocate lat/lon/lev')
    IF(doingMasking) THEN
     DEALLOCATE(regionNumbers, mask, STAT=ios)
     IF ( ios /= 0 ) CALL die ( myname, ': Cannot deallocate masking tape.')
    END IF
    return

CONTAINS

    subroutine Regrid_ (Gptr2file, im, jm, Gptr2bundle, im0, jm0, &
                            LONSfile, LONSbundle, verb ) 
    integer, intent(in) :: im, jm, im0, jm0
    real, pointer :: Gptr2bundle(:,:), Gptr2file(:,:)

    real, pointer :: LONSfile(:), LONSbundle(:)
    logical, intent(inout) :: verb

!   Local space to cope with MAPLE R4 conventions
!   ---------------------------------------------
    real(KIND=ESMF_KIND_R4) :: r4Gptr2bundle(im,jm), r4Gptr2file(im,jm)

 !   Code borrowed from MAPL_CFIO - Note (im,jm) here differs from above
!   Notice that since ARCTAS unless ForceBinning has been specified as
!   .FALSE. *binning* is always the regridding  method.
!                             ---
    integer :: STATUS

   
!ams    print *, 'Regrid_:  im, jm  (file)   = ', im, jm
!ams    print *, 'Regrid_: im0, jm0 (bundle) = ', im0, jm0
    if ( im /= im0 .OR. jm /= JM0 ) then
      r4Gptr2file(:,:) = Gptr2file(:,:) ! Only needed for GFS/GOCART 
      if (IM0 <  IM .or. JM0 < JM .or. ForceBinning_ ) then
          if ( verb ) print *, myname // ': Binning... '
          call MAPL_HorzTransformRun(Trans, r4Gptr2file, r4Gptr2bundle, MAPL_undef, rc=STATUS )
          VERIFY_(STATUS)
       else
          if ( verb ) print *, myname // ': Interpolating... '
          call hinterp ( r4Gptr2file,im,jm, r4Gptr2bundle,im0,jm0,1,MAPL_Undef )
       end if ! coarsening
       Gptr2bundle(:,:) = r4Gptr2bundle(:,:) ! to cope w/ MAPL r4 conventions
  else
       Gptr2bundle = Gptr2file
    end if ! change resolution
    if ( abs(LONSbundle(1)-LONSfile(1)+180.) <= 10.*tiny(1.) ) then
       if ( verb ) print *, myname // &
            ': shifting input longitudes by 180 degrees'
       call shift180Lon2D_ ( Gptr2bundle, im0, jm0 )
    end if
    end subroutine Regrid_


    subroutine shift180Lon2D_ ( c, im, jm )
    integer, intent(in) :: im, jm
    real, intent(inout) :: c(im,jm)
    real :: cj(im)
    integer :: m(4), n(4), imh, j
    imh = nint(im/2.)
    m = (/ 1,      imh, 1+imh,    im   /)
    n = (/ 1,   im-imh, 1+im-imh, im   /)
    do j = 1, jm
       cj(n(1):n(2)) = c(m(3):m(4),j)
       cj(n(3):n(4)) = c(m(1):m(2),j)
       c(:,j) = cj
    end do
    return
    end subroutine shift180Lon2D_

    subroutine queryFile_ ( fid, im, jm, km, lon, lat, lev, varn )
      integer, intent(in)  :: fid
      integer, intent(out) :: im, jm, km
      character(len=255), intent(inout) :: varn
      real, pointer :: lon(:), lat(:), lev(:)
!                    ----
      character(len=255)              :: title, source, contact, levunits
      character(len=255), allocatable :: vtitle(:), vunits(:)
      character(len=257), allocatable :: vname(:)
    
      real,    allocatable :: valid_range(:,:), packing_range(:,:)
      integer, allocatable :: kmvar(:), yyyymmdd(:), hhmmss(:)
      integer :: lm, nvars, timinc, ngatts, n
      real    :: amiss

      call GFIO_DimInquire ( fid, im, jm, km, lm, nvars, ngatts, rc )
      if ( rc /= 0 ) call die(myname,'cannot get file dimensions')
    
      allocate ( lat(jm), lon(im), lev(km), yyyymmdd(lm), hhmmss(lm),    &
              vname(nvars), vunits(nvars), vtitle(nvars), kmvar(nvars), &
              valid_range(2,nvars), packing_range(2,nvars),             &
              stat=rc )
      if ( rc /= 0 ) call die(myname,'cannot allocate file attribs' )

      call GFIO_Inquire ( fid, im, jm, km, lm, nvars,     &
         title, source, contact, amiss,  &
         lon, lat, lev, levunits,        &
         yyyymmdd, hhmmss, timinc,       &
         vname, vtitle, vunits, kmvar,   &
         valid_range , packing_range, rc )
      if ( rc /= 0 ) call die(myname,'cannot get file dimensions')
    
!     Make variable names case insensitive
!     ------------------------------------
      do n = 1, nvars
         if ( uppercase(trim(varn)) .EQ. uppercase(trim(vname(n))) ) then
              varn = trim(vname(n))
              exit
         end if
      end do 

      deallocate ( yyyymmdd, hhmmss,                  &
                  vname, vunits, vtitle, kmvar,      &
                  valid_range, packing_range,        &
                  stat=rc )
      
    end subroutine queryFile_

end subroutine Chem_UtilMPread

subroutine Chem_UtilGetTimeInfo ( fname, begDate, begTime, nTimes, incSecs )

  implicit NONE
  character(len=*), intent(in) :: fname    ! GFIO/CFIO filename
  integer, intent(out) :: begDate, begTime ! initial time/date on file
                                           ! given as YYYYMMDD and HHMMSS
  integer, intent(out) :: nTimes           ! number of time steps on file
  integer, intent(out) :: incSecs          ! time steps in seconds
 

  integer :: READ_ONLY=1
  integer :: fid, rc, im,jm,km,nvars,ngatts
  character(len=*), parameter ::  myname = 'Chem_UtilGetTimeInfo'
  
! Special case
! -----------
  if ( fname(1:9) == '/dev/null' ) then    
     begDate=0; begTime=0;  nTimes=0; incSecs=0
     return
  end if

! Open file
! ---------
  call GFIO_Open ( fname, READ_ONLY, fid, rc )
  if ( rc /= 0 ) call die(myname, 'Unable to open '// trim(fname))
  if ( rc .ne. 0 ) then
     begDate=-1; begTime=-1;  nTimes=-1; incSecs=-1
     return
  endif

! Get dimension sizes
! -------------------
  call GFIO_DimInquire (fid,im,jm,km,nTimes,nvars,ngatts,rc)
  if ( rc /= 0 ) call die(myname, 'Unable to inquire about dimension sizes for file '// trim(fname))
  if ( rc .ne. 0 ) then
     begDate=-1; begTime=-1;  nTimes=-1; incSecs=-1
     return
  endif

! Get initial time/timestep
! -------------------------
  call GetBegDateTime ( fid, begDate, begTime, incSecs, rc )
  if ( rc .ne. 0 ) then
     begDate=-1; begTime=-1;  nTimes=-1; incSecs=-1
     return
  endif

  call GFIO_close(fid,rc)

end subroutine Chem_UtilGetTimeInfo

!............................... geos4 stubs  ..........................

! Parallelized utility routine for computing/printing
! max/min of an input array
!
      subroutine pmaxmin3d ( qname, a, pmin, pmax, im, jt, fac )
      implicit none
      character*(*)  qname
      integer im, jt
      real :: a(:,:,:)
      real pmax, pmin
      real fac                     ! multiplication factor
      call pmaxmin2d ( qname, reshape(a,(/ im, jt /)), &
                             pmin, pmax, im, jt, fac )
      
      end subroutine pmaxmin3d

      subroutine pmaxmin2d ( qname, a, pmin, pmax, im, jt, fac )

      use parutilitiesmodule, only : commglobal, gid, maxop, parcollective

      implicit none

      character*(*)  qname
      integer im, jt
      real a(im,jt)
      real pmax, pmin
      real fac                     ! multiplication factor

      integer :: i, j, maxop_, two=2

      real qmin(jt), qmax(jt)
      real pm1(2)

      character(len=16) :: name

!$omp parallel do private(i, j, pmax, pmin)

      do j=1,jt
         pmax = a(1,j)
         pmin = a(1,j)
         do i=2,im
            pmax = max(pmax, a(i,j))
            pmin = min(pmin, a(i,j))
         enddo
         qmax(j) = pmax
         qmin(j) = pmin
      enddo
!
! Now find max/min of amax/amin
!
      pmax = qmax(1)
      pmin = qmin(1)
      do j=2,jt
         pmax = max(pmax, qmax(j))
         pmin = min(pmin, qmin(j))
      enddo

      pm1(1) = pmax
      pm1(2) = -pmin
      maxop_ = maxop
      call parcollective(commglobal, maxop_, two, pm1  )
      pmax=pm1(1)
      pmin=-pm1(2)
     
      if ( fac /= 0.0 ) then  ! trick to prevent printing
      if ( MAPL_am_I_root() ) then
           name = '            '
           name(1:len(qname)) = qname
           write(*,*) name, ' max = ', pmax*fac, ' min = ', pmin*fac
           return
      end if
      end if

    end subroutine pmaxmin2d


      function leap_year(ny)
!
! Determine if year ny is a leap year
!
! Author: S.-J. Lin
      implicit none
      logical leap_year
      integer ny
      integer ny00

!
! No leap years prior to 0000
!
      parameter ( ny00 = 0000 )   ! The threshold for starting leap-year 

      if( ny >= ny00 ) then
         if( mod(ny,100) == 0. .and. mod(ny,400) == 0. ) then
             leap_year = .true.
         elseif( mod(ny,4) == 0. .and. mod(ny,100) /= 0.  ) then
             leap_year = .true.
         else
             leap_year = .false.
         endif
      else
          leap_year = .false.
      endif

      return 
    end function leap_year


      integer FUNCTION INCYMD (NYMD,M)

!  PURPOSE
!     INCYMD:  NYMD CHANGED BY ONE DAY
!     MODYMD:  NYMD CONVERTED TO JULIAN DATE
!  DESCRIPTION OF PARAMETERS
!     NYMD     CURRENT DATE IN YYMMDD FORMAT
!     M        +/- 1 (DAY ADJUSTMENT)

      integer nymd, m, ny00, ny, nm, nd

      INTEGER NDPM(12)
      DATA    NDPM /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!!!      logical leap_year
      DATA    NY00     / 1900 /

      NY = NYMD / 10000
      NM = MOD(NYMD,10000) / 100
      ND = MOD(NYMD,100) + M

      IF (ND.EQ.0) THEN
      NM = NM - 1
      IF (NM.EQ.0) THEN
          NM = 12
          NY = NY - 1
      ENDIF
      ND = NDPM(NM)
      IF (NM.EQ.2 .AND. leap_year(NY))  ND = 29
      ENDIF

      IF (ND.EQ.29 .AND. NM.EQ.2 .AND. leap_year(ny))  GO TO 20

      IF (ND.GT.NDPM(NM)) THEN
      ND = 1
      NM = NM + 1
      IF (NM.GT.12) THEN
          NM = 1
          NY = NY + 1
      ENDIF
      ENDIF

   20 CONTINUE
      INCYMD = NY*10000 + NM*100 + ND
      RETURN
    END FUNCTION INCYMD

      subroutine tick (nymd, nhms, ndt)

! Input:
      integer ndt                     ! TIME-STEP
! Inpuit/Output:
      integer nymd                    ! CURRENT YYYYMMDD
      integer nhms                    ! CURRENT HHMMSS
!!!      integer incymd

! Revision:   S.-J. Lin Mar 2000
       integer nsecf, nhmsf, n, nsec

       NSECF(N)   = N/10000*3600 + MOD(N,10000)/100* 60 + MOD(N,100)
       NHMSF(N)   = N/3600*10000 + MOD(N,3600 )/ 60*100 + MOD(N, 60)

       NSEC = NSECF(NHMS) + ndt

       IF (NSEC.GT.86400)  THEN
           DO WHILE (NSEC.GT.86400)
              NSEC = NSEC - 86400
              NYMD = INCYMD (NYMD,1)
           ENDDO
       ENDIF

       IF (NSEC.EQ.86400)  THEN
           NSEC = 0
           NYMD = INCYMD (NYMD,1)
       ENDIF

       IF (NSEC .LT. 0)  THEN
           DO WHILE (NSEC .LT. 0)
               NSEC = 86400 + NSEC
               NYMD = INCYMD (NYMD,-1)
           ENDDO
        ENDIF

          NHMS = NHMSF (NSEC)
      return
    end subroutine tick

      subroutine mcalday(nymd, nhms, calday)
      implicit none

! input:
      integer nymd
      integer nhms
! Output:
      real calday                    ! Julian day (1 to 366 for non-leap year)
                                     ! Julian day (-1 to -367 for   leap year)
! Local:
      logical leapyr                 ! leap_year?
!!!      logical leap_year
      real tsec
      integer n, nsecf, m, mm
      integer dd, ds
      integer days(12)
      integer ny

      data days /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
      nsecf(n)  = n/10000*3600 + mod(n,10000)/100* 60 + mod(n,100)

      ny = nymd / 10000
      mm = mod(nymd, 10000) / 100 
      dd = mod(nymd,   100)

      ds = dd -1

      if( mm .ne. 1) then
      do m=1, mm-1
         if( m.eq.2  .and. leap_year(ny) ) then 
             ds = ds + 29
         else
             ds = ds + days(m)
         endif
      enddo
      endif

      tsec = ds * 86400 + nsecf(nhms)

      calday = tsec / 86400.  + 1.
      if( leap_year(ny) ) calday = -calday

      return
    end subroutine mcalday

      subroutine zenith(calday  ,dodiavg ,clat    ,coszrs  )

!
! Input arguments
!
      real calday              ! Calendar day, including fraction
      logical dodiavg          ! true => do diurnal averaging
      real clat                ! Current latitude (radians)
!
! Output arguments
!
      real coszrs(*)       ! Cosine solar zenith angle
!
!---------------------------Local variables-----------------------------
!

      call die ('zenith','stub only, please do not call zenith_()' )

    end subroutine zenith

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_UtilNegFiller --- Negative Filler
!
! !INTERFACE:
!
   subroutine Chem_UtilNegFiller ( q, delp, in, jn, qmin )

! !USES:

  implicit NONE

! !INPUT/OUTPUT PARAMETERS:

  integer                    :: in, jn             ! number of local lon/lat 
  real, pointer              :: delp(:,:,:)
  real, OPTIONAL, intent(in) :: qmin

! !OUTPUT PARAMETERS:

  real, pointer :: q(:,:,:)     ! 3D tracer

! !DESCRIPTION: 
!
! !REVISION HISTORY: Makes sure tracer has no negative values. This is
!                    a "flat tax" algorithm: first negative values are
!  replaced with tiny() or user specified value. Then profiles are rescaled 
!  to preserve column mass, whenever possible. No mass conservation is
!  imposed when the initial column mass is negative or zero.
!
!  18May2007  da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

  real :: mass_1(in,jn), mass_2(in,jn), ratio(in,jn)
  real :: qmin_, vmax, vmin
  integer :: k, k1, k2, km

  k1 = lbound(q,3)
  k2 = ubound(q,3)
  km = k2 - k1 + 1

! Unless specified, minimum is the smallest positive float, not zero
! ------------------------------------------------------------------
  if ( present(qmin) ) then
     qmin_ = qmin
  else  
     qmin_ = tiny(1.0)
  end if

#ifdef DEBUG
  call pmaxmin ( 'NegFill:  q_beg', q, vmax, vmin, in*jn, km, 1. )
#endif

! Column mass before fixer
! ------------------------
  mass_1 = sum ( delp * q, 3 )

! Cap q
! -----
  where ( q < qmin_ ) q = qmin_

! Enforce conservation of column mass
! -----------------------------------
  mass_2 = sum ( delp * q, 3 )
  where ( (mass_2 /= mass_1) .AND. (mass_1 > 0.0) )
          ratio = mass_1 / mass_2
  elsewhere
          ratio = 1.0
  end where

! Next correct q in each layer
! ----------------------------
  do k = k1, k2
     where ( ratio /= 1.0 )
             q(:,:,k) = ratio * q(:,:,k)
     end where
  end do

#ifdef DEBUG
  call pmaxmin ( 'NegFill: mass_1', mass_1, vmax, vmin, in*jn,1, 1. )
  call pmaxmin ( 'NegFill: mass_2', mass_2, vmax, vmin, in*jn,1, 1. )
  call pmaxmin ( 'NegFill:  ratio',  ratio, vmax, vmin, in*jn,1, 1. )
  call pmaxmin ( 'NegFill:  q_end', q, vmax, vmin, in*jn, km, 1. )
#endif

end subroutine Chem_UtilNegFiller

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_UtilTroppFixer - Repair tropopause pressure bad values
!
! !INTERFACE:
!
  SUBROUTINE Chem_UtilTroppFixer(im, jm, tropp, threshold, verbose, newTropp, rc)

! !USES:

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

  INTEGER, INTENT(IN)		:: im          !Index range, longitude
  INTEGER, INTENT(IN)		:: jm          !Index range, latitude
  REAL, POINTER , INTENT(INOUT)	:: tropp(:,:)  !Tropopause pressure (Pa)
  REAL, OPTIONAL, INTENT(IN)    :: threshold   !User-supplied threshold (Pa).  If
                                               ! not present, default is 110000 Pa.
  LOGICAL, OPTIONAL, INTENT(IN) :: verbose     !Write message when bad value is 
                                               ! encountered. DEBUG directive turns
                                               ! on the message even if verbose is not
                                               ! present or if verbose = .FALSE.
  INTEGER, INTENT(OUT)		:: rc          !Return code: 0=OK, 1=Unable to repair


! !OUTPUT PARAMETERS:

  REAL, OPTIONAL, INTENT(OUT)   :: newTropp(1:im,1:jm)
                                               !If present, fill with repaired
                                               ! tropopause pressures (Pa).  Otherwise
                                               ! tropp is overwritten.

! !DESCRIPTION: 
!
! Replace bad values in tropopause pressure array with an average of nearby good
! values.  Working outward from the afflicted cell, nearby cells are scanned for 
! the presence of good values. In the first scan, all adjacent cells, which can
! number up to eight, are considered.  If at least one of the cells has a valid
! pressure, the scanning is terminated.  Otherwise, the "radius" of the scan is
! increased by one, and up to 24 cells are considered, and so on, until, at the
! extreme, all cells on the current processor fall under consideration.
!
! After the scanning is done, the bad value is replaced with the average the 
! valid pressures found by the scan.  Thus, the accuracy of the replaced value
! drops rapidly as the scan expands outward.  At the same time, the only case
! in which a valid pressure will not be found to replace a bad one is if all 
! pressures on the current processor are invalid.  In this case the return 
! code is set to 1.
!
! If newTropp is not present, then the input array, tropp, is overwritten.  If
! it is present, it is filled, even if there are no invalid pressures in tropp.
!
! No scanning is done if all tropopause pressures are valid.
!
! Assumptions/bugs:
!
! Bad values are high, but not Indef, and the primary purpose of this routine
! is to repair the case where GEOS-5 fails to find the tropopause and assigns
! MAPL_UNDEF as the tropopause pressure.  We recommend using the blended 
! tropopause values for tropp, because the frequency of bad values is quite
! rare compared to the unblended case.
!
! !REVISION HISTORY: 
!
!  25Jan2008  Nielsen  Initial code and testing.
!
!EOP
!-------------------------------------------------------------------------
  CHARACTER(LEN=*), PARAMETER :: myName = 'Chem_UtilTroppFixer'

  LOGICAL :: tellMe

  INTEGER :: i,ier,j,l,m
  INTEGER :: ie,iw,jn,js

  INTEGER, ALLOCATABLE :: mask(:,:)
  INTEGER, ALLOCATABLE :: mx(:)
  REAL, ALLOCATABLE :: p(:,:)
  
  REAL :: badValue, r

  rc = 0

! Determine verbosity, letting the DEBUG 
! directive override local specification
! --------------------------------------
  tellMe = .FALSE.
  IF(PRESENT(verbose)) THEN
   IF(verbose) tellMe = .TRUE.
  END IF
#ifdef DEBUG
  tellMe = .TRUE.
#endif

! Set the bad value to 110000 Pa (1100 hPa)
! -----------------------------------------
  IF(PRESENT(threshold)) THEN
   badValue = threshold
  ELSE
   badValue = 1.10E+05 !Pa, 1100 hPa
  END IF

! There may be no bad values ...
! ------------------------------
  IF( ALL( tropp(1:im,1:jm) < badValue ) ) THEN
   IF(PRESENT(newTropp)) newTropp(1:im,1:jm) = tropp(1:im,1:jm)
   RETURN
  END IF

! ... or there is at least one bad value
! --------------------------------------
  ALLOCATE(mask(1:im,1:jm),STAT=ier)
  ALLOCATE(p(1:im,1:jm),STAT=ier)
  ALLOCATE(mx(4),STAT=ier)

! Loop over each cell
! -------------------
  DO j=1,jm
   DO i=1,im

! Invalid pressure found at cell(i,j)
! -----------------------------------
    IF(tropp(i,j) >= badValue) THEN

! Determine maximum "radius" of search
! ------------------------------------
     mx(1) = im-i
     mx(2) = i-1
     mx(3) = jm-j
     mx(4) = j-1

! Start search
! ------------
     DO m=1,MAXVAL(mx)

! Clear the mask
! --------------
      mask(1:im,1:jm) = 0

! Range of search
! ---------------
      iw = MAX( 1,i-m)
      ie = MIN(im,i+m)
      js = MAX( 1,j-m)
      jn = MIN(jm,j+m)

! Set mask to one for cells in range of search
! --------------------------------------------
       mask(iw:ie,js:jn) = 1

! Set mask back to zero for cells in range
! of search that have invalid pressures.
! ----------------------------------------
       WHERE(tropp(iw:ie,js:jn) >= badValue) mask(iw:ie,js:jn) = 0

! One valid pressure is enough ...
! --------------------------------
       IF(SUM(MASK) >= 1) EXIT

! ... or "radius" of search needs to be extended
! ----------------------------------------------
     END DO

! Repair bad value at cell(i,j) with average 
! of valid pressures found in range of search
! -------------------------------------------
     r = SUM(tropp,mask == 1)
     p(i,j) = r/(1.00*SUM(mask))

! For debugging
! -------------
     IF(tellMe) THEN
      WRITE(*,FMT="(A,': ',ES12.5,' becomes ',ES12.5,' Pa [',I4,2X,I4,']')") &
            TRIM(myName),tropp(i,j),p(i,j),m,SUM(mask)
     END IF

    ELSE

! Input pressure at cell(i,j) was valid
! -------------------------------------
     p(i,j) = tropp(i,j)

    END IF

! Next cell
! ---------
   END DO
  END DO

! Clean up
! --------
  DEALLOCATE(mask,STAT=ier)
  DEALLOCATE(mx,STAT=ier)

! If all cells have bad values, then 
! register a failure, but continue.
! ----------------------------------
  IF( ANY( p(1:im,1:jm) >= badValue ) ) THEN
   PRINT *, myName,": WARNING Unable to fix bad tropopause pressure(s)"
   rc = 1
  END IF
  
! Overwrite input or fill output array
! ------------------------------------
  IF(PRESENT(newTropp)) THEN
   newTropp(1:im,1:jm) = p(1:im,1:jm)
  ELSE
   tropp(1:im,1:jm) = p(1:im,1:jm)
  END IF

! Clean up some more
! ------------------
  DEALLOCATE(p,STAT=ier)

  RETURN
  END SUBROUTINE Chem_UtilTroppFixer

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !ROUTINE:  Chem_UtilExtractIntegers - Extract integers from a delimited string
!
! !INTERFACE:
!
  SUBROUTINE Chem_UtilExtractIntegers(string,iSize,iValues,delimiter,verbose,rc)

! !USES:

  IMPLICIT NONE

! !INPUT/OUTPUT PARAMETERS:

  CHARACTER(LEN=*), INTENT(IN)   :: string	  ! Character-delimited string of integers
  INTEGER, INTENT(IN)            :: iSize
  INTEGER, INTENT(INOUT)         :: iValues(iSize)! Space allocated for extracted integers
  CHARACTER(LEN=*), OPTIONAL     :: delimiter     ! 1-character delimiter
  LOGICAL, OPTIONAL, INTENT(IN)  :: verbose	  ! Let me know iValues as they are found. 
                                		  ! DEBUG directive turns on the message even 
                                		  ! if verbose is not present or if 
                                		  ! verbose = .FALSE.
  INTEGER, OPTIONAL, INTENT(OUT) :: rc            ! Return code

! !DESCRIPTION: 
!
!  Extract integers from a character-delimited string, for example, "-1,45,256,7,10".  In the context
!  of Chem_Util, this is provided for determining the numerically indexed regions over which an 
!  emission might be applied.
!
!  In multiple passes, the string is parsed for the delimiter, and the characters up to, but not
!  including the delimiter are taken as consecutive digits of an integer.  A negative sign ("-") is
!  allowed.  After the first pass, each integer and its trailing delimiter are lopped of the head of
!  the (local copy of the) string, and the process is started over.
!
!  The default delimiter is a comma (",").
!
!  "Unfilled" iValues are zero.
!  
!  Return codes:
!  1 Zero-length string.
!  2 iSize needs to be increased.
!
!  Assumptions/bugs:
!
!  A non-zero return code does not stop execution.
!  Allowed numerals are: 0,1,2,3,4,5,6,7,8,9.
!  A delimiter must be separated from another delimiter by at least one numeral.
!  The delimiter cannot be a numeral or a negative sign.
!  The character following a negative sign must be an allowed numeral.
!  The first character must be an allowed numeral or a negative sign.
!  The last character must be an allowed numeral.
!  The blank character (" ") cannot serve as a delimiter.
!
!  Examples of strings that will work:
!  "1"
!  "-1"
!  "-1,2004,-3"
!  "1+-2+3"
!  "-1A100A5"
!
!  Examples of strings that will not work:
!  "1,--2,3"
!  "1,,2,3"
!  "1,A,3"
!  "1,-,2"
!  "1,2,3,4,"
!  "+1"
!  "1 3 6"
!
! !REVISION HISTORY: 
!
!  29Feb2008  Nielsen  Initial code and testing.
!
!EOP
!-------------------------------------------------------------------------
 CHARACTER(LEN=*), PARAMETER :: myName = 'Chem_UtilExtractIntegers'

 INTEGER :: base,count,i,iDash,last,lenStr
 INTEGER :: multiplier,pos,posDelim,sign
 CHARACTER(LEN=255) :: str
 CHARACTER(LEN=1) :: char,delimChar
 LOGICAL :: Done
 LOGICAL :: tellMe

! Initializations
! ---------------
 rc = 0
 count = 1
 Done = .FALSE.
 iValues(:) = 0
 base = ICHAR("0")
 iDash = ICHAR("-")

! Determine verbosity, letting the DEBUG 
! directive override local specification
! --------------------------------------
  tellMe = .FALSE.
  IF(PRESENT(verbose)) THEN
   IF(verbose) tellMe = .TRUE.
 END IF
#ifdef DEBUG
  tellMe = .TRUE.
#endif

! Check for zero-length string
! ----------------------------
 lenStr = LEN_TRIM(string)
 IF(lenStr == 0) THEN
  rc = 1
  PRINT *,myname,": ERROR - Found zero-length string."
  RETURN
 END IF

! Default delimiter is a comma
! ----------------------------
 delimChar = ","
 IF(PRESENT(delimiter)) delimChar(1:1) = delimiter(1:1)

! Work on a local copy
! --------------------
 str = TRIM(string)

! One pass for each delimited integer
! -----------------------------------
 Parse: DO

  lenStr = LEN_TRIM(str)

! Parse the string for the delimiter
! ----------------------------------
  posDelim = INDEX(TRIM(str),TRIM(delimChar))
  IF(tellMe) PRINT *,myname,": Input string is >",TRIM(string),"<"

! If the delimiter does not exist,
! one integer remains to be extracted.
! ------------------------------------
  IF(posDelim == 0) THEN
   Done = .TRUE.
   last = lenStr
  ELSE
   last = posDelim-1
  END IF
  multiplier = 10**last

! Examine the characters of this integer
! --------------------------------------
  Extract: DO pos=1,last

   char = str(pos:pos)
   i = ICHAR(char)

! Account for a leading "-"
! -------------------------
   IF(pos == 1) THEN 
    IF(i == iDash) THEN
     sign = -1
    ELSE
     sign = 1
    END IF
   END IF

! "Power" of 10 for this character
! --------------------------------
   multiplier = multiplier/10

   IF(pos == 1 .AND. sign == -1) CYCLE Extract

! Integer comes from remaining characters
! ---------------------------------------
   i = (i-base)*multiplier
   iValues(count) = iValues(count)+i
   IF(pos == last) THEN
    iValues(count) = iValues(count)*sign
    IF(tellMe) PRINT *,myname,":Integer number ",count," is ",iValues(count)
   END IF
  
  END DO Extract

  IF(Done) EXIT

! Lop off the leading integer and try again
! -----------------------------------------
  str(1:lenStr-posDelim) = str(posDelim+1:lenStr)
  str(lenStr-posDelim+1:255) = " "
  count = count+1

! Check size
! ----------
  IF(count > iSize) THEN
   rc = 2
   PRINT *,myname,": ERROR - iValues does not have enough elements."
  END IF

 END DO Parse

 RETURN
END SUBROUTINE Chem_UtilExtractIntegers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This is a candidate for ESMFL, here for dependency reasons
!  

  subroutine GridGetLatLons_ ( grid, lons, lats )

    implicit NONE
    type(ESMF_Grid) :: grid

    real, pointer   :: lons(:), lats(:)

!                     ---

    character(len=*), parameter :: Iam = 'GridGetLatLons'

    real(KIND=8), pointer  :: R8D2(:,:)
    real, pointer          :: lons2d(:,:), lats2d(:,:)
    real, pointer          :: LONSLocal(:,:), LATSlocal(:,:)
    integer                :: IM_WORLD, JM_WORLD, dims(3), STATUS, RC

!                          ----

!      Get world dimensions
!      --------------------
       call MAPL_GridGet(GRID, globalCellCountPerDim=DIMS, RC=STATUS)

       IM_WORLD = dims(1)
       JM_WORLD = dims(2)

!      Allocate memory for output if necessary
!      ---------------------------------------
       if ( .not. associated(lons) ) then
            allocate(lons(IM_WORLD), stat=STATUS)
       else
            if(size(LONS,1) /= IM_WORLD) STATUS = 1
       end if
       VERIFY_(status)
       if ( .not. associated(lats) ) then
            allocate(lats(JM_WORLD), stat=STATUS)
       else
            if(size(LATS,1) /= JM_WORLD) STATUS = 1
       end if
       VERIFY_(status)

!      Local work space
!      ----------------
       allocate(LONS2d(IM_WORLD,JM_WORLD), LATS2d(IM_WORLD,JM_WORLD), &
                STAT=status)             
       VERIFY_(status)
       LONS2d=0
       LATS2d=0

!      Get the local longitudes and gather them into a global array
!      ------------------------------------------------------------
       call ESMF_GridGetCoord(grid, localDE=0, coordDim=1, &
             staggerloc=ESMF_STAGGERLOC_CENTER, doCopy=ESMF_DATA_REF, &
             fptr=R8D2, rc=status)

       allocate(LONSLOCAL(size(R8D2,1),size(R8D2,2)), STAT=status)             
       VERIFY_(status)

       LONSLOCAL = R8D2*(180/MAPL_PI)

       call ArrayGather(LONSLOCAL, LONS2D, GRID, RC=STATUS)

!      Get the local longitudes and gather them into a global array
!      ------------------------------------------------------------
       call ESMF_GridGetCoord(grid, localDE=0, coordDim=2, &
             staggerloc=ESMF_STAGGERLOC_CENTER, doCopy=ESMF_DATA_REF, &
             fptr=R8D2, rc=status)

       allocate(LATSLOCAL(size(R8D2,1),size(R8D2,2)), STAT=status)             
       VERIFY_(status)

       LATSlocal = R8D2*(180/MAPL_PI)

       call ArrayGather(LATSLOCAL, LATS2D, GRID, RC=STATUS)
       VERIFY_(STATUS)

!      Return 1D arrays
!      ----------------
       LONS = LONS2D(:,1)
       LATS = LATS2D(1,:)

       DEALLOCATE(LONSLOCAL, LATSLOCAL, LONS2d, LATS2d )
        
     end subroutine GridGetLatLons_

 end module Chem_UtilMod

