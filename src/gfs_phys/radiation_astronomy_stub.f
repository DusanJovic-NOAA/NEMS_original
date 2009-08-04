!!!!!  ==========================================================  !!!!!
!!!!!          'module_radiation_astronomy'  description           !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   set up astronomy quantities for solar radiation calculations.      !
!                                                                      !
!   in module 'module_radiation_astronomy', externally accessable      !
!   subroutines are listed below:                                      !
!                                                                      !
!      'solinit'    -- read in solar constant                          !
!         input:                                                       !
!           ( ISOL, iyear, me )                                        !
!         output:                                                      !
!           ( none )                                                   !
!                                                                      !
!      'astronomy'  -- get astronomy related quantities                !
!         input:                                                       !
!           ( lons_lar,glb_lats_r,sinlat,coslat,xlon,                  !
!!            fhswr,jdate,deltim,                                      !
!             fhswr,jdate,                                             !
!             LON2,LATD,LATR,IPT_LATR, lsswr, me)                      !
!         output:                                                      !
!           ( solcon,slag,sdec,cdec,coszen,coszdg)                     !
!                                                                      !
!                                                                      !
!   external modules referenced:                                       !
!       'module machine'                    in 'machine.f'             !
!       'module physcons'                   in 'physcons.f             !
!                                                                      !
!   program history log:                                               !
!     may-06-1977  ---  ray orzol,      created at gfdl                !
!     jul-07-1989  ---  kenneth campana                                !
!     may-15-1998  ---  mark iredell    y2k compliance                 !
!     dec-15-2003  ---  yu-tai hou      combined compjd and fcstim and !
!                       rewrite in fortran 90 compatable form          !
!     feb-15-2006  ---  yu-tai hou      add 11-yr solar constant cycle !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radiation_astronomy  !
!........................................!
!
      implicit   none
!
      public  astronomy


! =================
      contains
! =================
!-----------------------------------
      subroutine astronomy                                              &
!...................................

!  ---  inputs:
     &     ( lons_lar,glb_lats_r,sinlat,coslat,xlon,                    &
!    &       fhswr,jdate,deltim,                                        &
     &       fhswr,jdate,                                               &
     &       LON2,LATD,LATR,IPT_LATR, lsswr, me,                        &
!  ---  outputs:
     &       solcon,slag,sdec,cdec,coszen,coszdg                        &
     &      )

!  ===================================================================  !
!                                                                       !
!  astronomy computes solar parameters at forecast time                 !
!                                                                       !
!  inputs:                                                   dimension  !
!    lons_lar      - num of grid pts on a given lat circle        (LATR)!
!    glb_lats_r    - index for global latitudes                   (LATR)!
!    sinlat,coslat - sin and cos of latitude                      (LATR)!
!    xlon          - longitude in radians                    (LON2*LATD)!
!    fhswr         - sw radiation calling interval in hour              !
!    jdate         - current forecast date and time               (8)   !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)      !
!!   deltim        - duration of model integration time step in seconds !
!    LON2,LATD,LATR- dimensions for longitude/latitude directions       !
!    IPT_LATR      - latitude index location indecator                  !
!    lsswr         - logical control flag for sw radiation call         !
!    me            - integer control flag for diagnostic print out      !
!                                                                       !
!  outputs:                                                             !
!    solcon        - sun-earth distance adjusted solar constant (w/m2)  !
!    slag          - equation of time in radians                        !
!    sdec, cdec    - sin and cos of the solar declination angle         !
!    coszen        - avg of cosz for daytime only            (LON2,LATD)!
!    coszdg        - avg of cosz over entire sw call interval(LON2,LATD)!
!                                                                       !
!                                                                       !
!  external functions called: iw3jdn                                    !
!                                                                       !
!  ===================================================================  !
!
      implicit none
      
!  ---  input:
      integer,  intent(in) :: LON2, LATD, LATR, IPT_LATR, me
      integer,  intent(in) :: lons_lar(:), glb_lats_r(:), jdate(:)

      logical, intent(in) :: lsswr

      real (kind=8), intent(in) :: sinlat(:), coslat(:),        &
     &       xlon(:,:), fhswr
!    &       xlon(:,:), fhswr, deltim

!  ---  output:
      real (kind=8), intent(out) :: solcon, slag, sdec, cdec,   &
     &       coszen(:,:), coszdg(:,:)

!
      return
!...................................
      end subroutine astronomy
!-----------------------------------
!
!...........................................!
      end module module_radiation_astronomy !
!===========================================!
