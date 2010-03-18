
                        module module_control
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!use module_include
use module_dm_parallel,only : ids,ide,jds,jde &
                             ,ims,ime,jms,jme &
                             ,its,ite,jts,jte &
                             ,lm &
                             ,mype_share,npes,num_pts_max &
                             ,mpi_comm_comp &
                             ,dstrb
use module_exchange
use module_constants
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      public NMMB_Finalize
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---look-up tables------------------------------------------------------
!-----------------------------------------------------------------------
integer,parameter :: &
 itb=201 &                   ! convection tables, dimension 1
,jtb=601 &                   ! convection tables, dimension 2
,kexm=10001 &                ! size of exponentiation table
,kztm=10001                  ! size of surface layer stability function table

real,parameter :: &
 ph=105000. &                ! upper bound of pressure range
,thh=350. &                  ! upper bound of potential temperature range
,thl=200.                    ! upper bound of potential temperature range

integer:: &
 kexm2 &                     ! internal points of exponentiation table
,kztm2                       ! internal pts. of the stability function table

integer,private :: &
 mype

real :: &
 dex &                       ! exponentiation table step
,dzeta1 &                    ! sea table z/L step
,dzeta2 &                    ! land table z/L step 
,fh01 &                      ! prandtl number for sea stability functions
,fh02 &                      ! prandtl number for land stability functions
,pl &                        ! lower bound of pressure range
,rdp &                       ! scaling factor for pressure
,rdq &                       ! scaling factor for humidity
,rdth &                      ! scaling factor for potential temperature
,rdthe &                     ! scaling factor for equivalent pot. temperature
,rdex &                      ! exponentiation table scaling factor
,xmax &                      ! upper bound for exponent in the table
,xmin &                      ! lower bound for exponent in the table
,ztmax1 &                    ! upper bound for z/L for sea stab. functions
,ztmin1 &                    ! lower bound for z/L for sea stab. functions
,ztmax2 &                    ! upper bound for z/L for land stab. functions
,ztmin2 &                    ! lower bound for z/L for land stab. functions
,dphd &
,dlmd
real,dimension(1:itb):: &
 sthe &                      ! range for equivalent potential temperature
,the0                        ! base for equivalent potential temperature           

real,dimension(1:jtb):: &
 qs0 &                       ! base for saturation specific humidity
,sqs                         ! range for saturation specific humidity

real,dimension(1:kexm):: &
 expf                        ! exponentiation table

real,dimension(1:kztm):: &
 psih1 &                     ! sea heat stability function
,psim1 &                     ! sea momentum stability function
,psih2 &                     ! land heat stability function
,psim2                       ! land momentum stability function

real,dimension(1:itb,1:jtb):: &
 ptbl                        ! saturation pressure table

real,dimension(1:jtb,1:itb):: &
 ttbl                        ! temperature table
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---miscellaneous control parameters------------------------------------
!-----------------------------------------------------------------------
character(64):: &
 infile
character(64):: &
 input_data
 
logical:: &
 adv_standard &              ! my task is in standard advec region
,adv_upstream &              ! my task is in upstream advec region
,first &                     ! if true preparation for the first time step
,hydro &                     ! if true hydrostatic dynamics
,read_global_sums &          ! read global sums in subroutine adv2 (bit identity)
,readbc &                    ! read regional boundary conditions
,run &                       ! initial data ready, start run
,runbc &                     ! boundary data ready, start run
,write_global_sums  &         ! write global sums in subroutine adv2 (bit identity)
,global 
integer:: &
 ierr &                      ! error code
,ihr &                       ! current forecast hour
,ihrbc &                     ! boundary condition hour
,ihrend &                    ! maximum forecast length, hours
,ihrst &                     ! forecast starting time
,ihrstbc &                   ! boundary conditions starting time
,nbc &                       ! boundary data logical unit
,nboco &                     ! time steps between updating boundary conditions
,ncnvc &                     ! time steps between convection calls
,nfcst &                     ! initial data logical unit
,nhours_fcst &               ! desired forecast length, hours
,nphs &                      ! time steps between physics calls
,nprec &                     ! time steps for precip accumulation
,nradl &                     ! time steps between longwave radiation calls
,nrads &                     ! time steps between shortwave radiation calls
,nrain &                     ! time steps between grid scale precip calls
,nradpp &                    ! time steps for precipitation in radiation
,ntsti &                     ! initial time step
,ntstm &                     ! final time step
,ntstm_max &                 ! maximum final timestep
,ntsd                        ! current time step

integer,dimension(1:3):: &
 idat(3) &                   ! date of initial data, day, month, year
,idatbc(3)                   ! date of boundary data, day, month, year

integer,dimension(1:15):: & 
 nfftrh &                    ! fft working field, h points
,nfftrw                      ! fft working field, v points
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---grid constants------------------------------------------------------
!-----------------------------------------------------------------------
integer:: &
 icycle &                    ! # of independent points, x direction
,im &                        ! maximum horizontal index, x direction
,jm &                        ! maximum horizontal index, y direction
!!!,lm &                        ! maximum vertical index, p direction
,lnsad &                     ! # of boundary lines with upstream advection
,lnsbc &                     ! # of boundary lines with enhanced diffusion
,lnsh &                      ! # of boundary h lines for bc in reg. setup
,lnsv &                      ! # of boundary v lines for bc in reg. setup
,lpt2 &                      ! # of pressure layers
,n2                          ! # starting address of 3d scratch fields

logical :: &                 
 e_bdy &                     ! logical flag for tasks on eastern boundary
,n_bdy &                     ! logical flag for tasks on northern boundary
,s_bdy &                     ! logical flag for tasks on southern boundary
,w_bdy                       ! logical flag for tasks on western boundary

real:: &
 bofac &                     ! amplification of diffusion along bndrs.
,ctph0 &                     ! cos(tph0)
,ddmpv &                     ! divergence damping coefficient, y direction
,dlm &                       ! grid size, delta lambda, radians
,dph &                      ! grid size, delta phi, degrees
,dt &                        ! dynamics time step
,dtq2 &                      ! turbulence time step
,dyh &                       ! delta y, h points
,dyv &                       ! delta y, v points 
,ef4t &                      ! grid constant
,f4d &                       ! grid constant
,pdtop &                     ! depth of pressure range in hybrid v. coord.
,pt &                        ! pressure at the top of model's atmosphere
,ptsgm &                     ! approx. pressure at the top of sigma range
,rdlm &                      ! 1 / delta lambda
,rdph &                      ! 1 / delta phi
,rdyh &                      ! 1 / delta y, h points
,rdyv &                      ! 1 / delta y, v points
,sb &                        ! radians from center to southern boundary
,sbd &                       ! degrees from center to southern boundary
,stph0 &                     ! sin(tph0)
,tboco &                     ! boundary conditions interval, hours
,tlm0 &                      ! radians grid rotated in lambda direction
,tlm0d &                     ! degrees grid rotated in lambda direction
,tph0 &                      ! radians grid rotated in phi direction
,tph0d &                     ! degrees grid rotated in phi direction
,wb &                        ! radians from center to western boundary
,wbd                         ! degrees from center to western boundary

!real,allocatable,dimension(:):: &      !lm
! dsg1 &                      ! thicnesses of sigma layers in pressure range
!,dsg2 &                      ! thicnesses of sigma layers in sigma range
!,pdsg1 &                     ! thicnesses of pressure layers in press. range
!,psgml1 &                    ! pressure at midlayers in pressure range
!,sgml1 &                     ! sigma at midlayers in pressure range
!,sgml2                       ! sigma at midlayers in sigma range

!real,allocatable,dimension(:):: &      !lm+1
!,psg1 &                      ! pressure at interfaces in pressure range
!,sg1 &                       ! sigma at interfaces in pressure range
!,sg2 &                       ! sigma at interfaces in sigma range
!,sgm                         ! sigma at interfaces

integer,allocatable,dimension(:):: &      !jm
 khfilt &                    ! polar filter, truncation wave #, h points
,kvfilt &                    ! polar filter, truncation wave #, v points
,nhsmud                      ! polar smoother for unfiltered variables

!real,allocatable,dimension(:):: &      !jm
! curv &                      ! curvature term in coriolis force
!,dare &                      ! gridbox area
!,ddmpu &                     ! divergence damping coefficient, x direction
!,ddv &                       ! gridbox diagonal distance 
!,dxh &                       ! delta x, h points
!,dxv &                       ! delta x, v points
!,fad &                       ! momentum advection factor
!,fah &                       ! z, w advection factor
!,fcp &                       ! temperature advection factor
!,fdiv &                      ! divergence factor
!,rare &                      ! 1 / gridbox area
!,rddv &                      ! 1 / gridbox diagonal distance 
!,rdxh &                      ! 1 / delta x, h points
!,rdxv &                      ! 1 / delta x, v points
!,wpdar                       ! weight of grid separaton filter

real,allocatable,dimension(:):: &      !2*(im-3)
 wfftrh &                    ! fft working field, h points
,wfftrw                      ! fft working field, v points

real,allocatable,dimension(:,:):: &      !im,jm
 f &                         ! coriolis parameter
!,glat &                      ! latitudes of h points
!,glon &                      ! longitudes of h points
!,hdac &                      ! lateral diffusion coefficient, h points
!,hdacv &                     ! lateral diffusion coefficient, v points
,hfilt &                     ! polar filter, h points
,vfilt                       ! polar filter, v points
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---fixed surface fields------------------------------------------------
!-----------------------------------------------------------------------
real,allocatable,dimension(:,:):: &      !im,jm
 albedo &                    ! base albedo
!,epsr &                      ! emissivity
!,fis &                       ! surface geopotential
!,sice &                      ! sea ice mask
!,sm &                        ! sea mask
!,stdh &                      ! standard deviation of topography height
 ,stdh                        ! standard deviation of topography height
!,sst &                       ! sea surface temperature
!,vegfrc &                    ! vegetation fraction
!,z0                          ! roughness
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---surface variables---------------------------------------------------
!-----------------------------------------------------------------------
real,allocatable,dimension(:,:):: &      !im,jm
 akhs &                      ! heat exchange coeff. / sfc. layer depth
,akms &                      ! momentum exchange coeff. / sfc. layer depth
,acprec &                    ! accumulated total precipitation
!,cldefi &                    ! convective cloud efficiency
,cuppt &                     ! convective precipitation for radiation
,cuprec &                    ! accumulated convective precipitation
,deep &                      ! deep convection mask, 0. or 1.
,ghct &                      ! countergradient correction
,prcu &                      ! physics time step convective liquid precip
,prec &                      ! total physics time step precipitation
,prra &                      ! physics time step grid scale liquid precip
,prsn &                      ! physics time step grid scale/conv. snow
,q02 &                       ! 2m specific humidity
,q10 &                       ! 10m specific humidity
,qs &                        ! specific humidity at the surface
,qz0 &                       ! spec. humidity at the top of viscous sublayer
,radin &                     ! total incoming radiation at the surface
,rlwin &                     ! downward lw at the surface
,rswin &                     ! downward sw at the surface
,roff &                      ! runoff
,sno &                       ! snow water equivalent
,th02 &                      ! 2m potential temperature
,th10 &                      ! 10m potential temperature
,ths &                       ! potential temperature at the surface
,thz0 &                      ! pot. temperature at the top of viscous sublayer
,u10 &                       ! 10m u
,ustar &                     ! u star
,uz0 &                       ! u at the top of viscous sublayer
,v10 &                       ! 10m v
,vz0 &                       ! v at the top of viscous sublayer
,wliq                        ! canopy moisture

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---regional boundary conditions----------------------------------------
!-----------------------------------------------------------------------
real,allocatable,dimension(:,:,:):: &      !im,lnsh,2
 pdbn &                      ! pressure difference at northern boundary
,pdbs                        ! pressure difference at southern boundary

real,allocatable,dimension(:,:,:):: &      !lnsh,jm,2
 pdbe &                      ! pressure difference at eastern boundary
,pdbw                        ! pressure difference at western boundary

real,allocatable,dimension(:,:,:,:):: &    !im,lnsh,lm,2
 tbn &                       ! temperature at northern boundary
,tbs &                       ! temperature at southern boundary
,qbn &                       ! specific humidity at northern boundary
,qbs &                       ! specific humidity at southern boundary
,wbn &                       ! condensate at northern boundary
,wbs                         ! condensate at southern boundary

real,allocatable,dimension(:,:,:,:):: &    !lnsh,jm,lm,2
 tbe &                       ! temperature at eastern boundary
,tbw &                       ! temperature at western boundary
,qbe &                       ! specific humidity at eastern boundary
,qbw &                       ! specific humidity at western boundary
,wbe &                       ! condensate at eastern boundary
,wbw                         ! condensate at western boundary

real,allocatable,dimension(:,:,:,:):: &    !im,lnsv,lm,2
 ubn &                       ! u wind component at northern boundary
,ubs &                       ! u wind component at southern boundary
,vbn &                       ! v wind component at northern boundary
,vbs                         ! v wind component at southern boundary

real,allocatable,dimension(:,:,:,:):: &    !lnsv,im,lm,2
 ube &                       ! u wind component at eastern boundary
,ubw &                       ! u wind component at western boundary
,vbe &                       ! v wind component at eastern boundary
,vbw                         ! v wind component at western boundary
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---atmospheric variables, hydrostatic----------------------------------
!-----------------------------------------------------------------------
!real,allocatable,dimension(:,:):: &      !im,jm
! pd &                        ! pressure difference, sigma range
!,pdo &                       ! previous pressure difference, sigma range
!,psdt                        ! hydrostatic surface pressure tendency
!real,allocatable,dimension(:,:,:):: &      !im,jm,lm
! cw &                        ! condensate
!,omgalf &                    ! omega-alpha
!,q &                         ! specific humidity
!,q2 &                        ! 2tke
!,o3 &                        ! ozone
!,t &                         ! temperature
!,tp &                        ! previous temperature
!,u &                         ! u wind component
!,up &                        ! previous u wind component
!,v &                         ! v wind component
!,vp                          ! previous v wind component

!real,allocatable,dimension(:,:,:):: &      !im,jm,lm+1
! pint                        ! nonhydrostatic pressure at interfaces
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---atmospheric variables, nonhydrostatic-------------------------------
!-----------------------------------------------------------------------
!real,allocatable,dimension(:,:,:):: &      !im,jm,lm
! dwdt &                      ! vertical acceleration, correction factor
!,pdwdt &                     ! previous correction factor
!,rtop &                      ! RT/p, specific volume
!,w &                         ! vertical velocity
!,z                           ! height at midlayers
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---physics variables---------------------------------------------------
!-----------------------------------------------------------------------
integer,allocatable,dimension(:,:):: &      !im,jm
 lbot &                      ! convective cloud base level
,lpbl &                      ! pbl top level
,ltop                        ! convective cloud top level
!-----------------------------------------------------------------------
!---scratch area for passing temporary arguments------------------------
!-----------------------------------------------------------------------
real,allocatable,dimension(:,:,:):: &      !im,jm,lm
 div &                       ! horizontal mass divergence
,e2 &                        ! 2*TKE in the layers
,pcne &                      ! second term of pgf, ne direction
,pcnw &                      ! second term of pgf, nw direction
,pcx &                       ! second term of pgf, x direction
,pcy &                       ! second term of pgf, y direction
,pfne &                      ! mass flux, ne direction
,pfnw &                      ! mass flux, nw direction
,pfx &                       ! mass flux, x direction
,pfy &                       ! mass flux, y direction
,tct                         ! time change of temperature
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
       contains
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine consts &
      (global &
      ,smag2,codamp,wcor &
      ,pt &
      ,tph0d,tlm0d &
      ,dphd,dlmd &
      ,dxh,rdxh &
      ,dxv,rdxv &
      ,dyh,rdyh &
      ,dyv,rdyv &
      ,ddv,rddv &
      ,ddmpu,ddmpv &
      ,ef4t,wpdar &
      ,fcp,fdiv &
      ,curv,f &
      ,fad,fah &
      ,dare,rare &
      ,glat,glon &
      ,vlat,vlon &
      ,hdacx,hdacy &
      ,hdacvx,hdacvy &
      ,lnsh,lnsad &
      ,nboco,tboco)
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
integer,intent(in) :: &
 lnsh

integer,intent(out) :: &
lnsad &
,nboco

real,intent(in) :: &
 codamp &    ! divergence damping coefficient
,pt &        ! Pressure at top of domain (Pa)
,smag2 &     ! Smagorinsky coefficient for 2nd order diffusion 
,tlm0d &
,tph0d &
,wcor

logical,intent(in) :: &
 global
 
real,intent(out) :: &
 ddmpv &
,dlmd &
,dphd &
,dyh &
,dyv &
,ef4t &
,rdyh &
,rdyv &
,tboco
 
real,dimension(jds:jde),intent(out) :: &
 curv &
,dare &
,ddmpu &
,ddv &
,dxh &
,dxv &
,fad &
,fah &
,fcp &
,fdiv &
,rare &
,rddv &
,rdxh &
,rdxv &
,wpdar
 
real,dimension(ims:ime,jms:jme),intent(out) :: &
 f &
,glat &
,glon &
,vlat &
,vlon &
,hdacx &
,hdacy &
,hdacvx &
,hdacvy
!
!-----------------------------------------------------------------------
!
                        end subroutine consts       
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
!
                        subroutine exptbl
!     ******************************************************************
!     *                                                                *
!     *               exponential function table                       *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer:: &
 k                           ! index

real:: &
 x &                         ! argument
,xrng                        ! argument range
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
                        end subroutine exptbl
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        function zjexp(x)
!     ******************************************************************
!     *                                                                *
!     *               exponential function table                       *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
real:: &
 zjexp


!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer:: &
 k                           ! index

real:: &
 ak &                        ! position in table
,x                           ! argument
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      ak=(x-xmin)*rdex
      k=int(ak)
      k=max(k,0)
      k=min(k,kexm2)
!
      zjexp=(expf(k+2)-expf(k+1))*(ak-real(k))+expf(k+1)


!-----------------------------------------------------------------------
!
                        end function zjexp
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine tablep
!     ******************************************************************
!     *                                                                *
!     *    generates the table for finding pressure from               *
!     *    saturation specific humidity and potential temperature      *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real,parameter:: &
 eps=1.e-10
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        end subroutine tablep
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine tablet
!     ******************************************************************
!     *                                                                *
!     *    generates the table for finding temperature from            *
!     *    pressure and equivalent potential temperature               *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real,parameter:: &
 eps=1.e-10
!-----------------------------------------------------------------------
!
                        end subroutine tablet
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine spline(jtb,nold,xold,yold,y2,nnew,xnew,ynew,p,q)           
!     ******************************************************************
!     *                                                                *
!     *  this is a one-dimensional cubic spline fitting routine        *
!     *  programed for a small scalar machine.                         *
!     *                                                                *
!     *  programer: z. janjic, yugoslav fed. hydromet. inst., beograd  *
!     *                                                                *
!     *                                                                *
!     *                                                                *
!     *  nold - number of given values of the function.  must be ge 3. *
!     *  xold - locations of the points at which the values of the     *
!     *         function are given.  must be in ascending order.       *
!     *  yold - the given values of the function at the points xold.   *
!     *  y2   - the second derivatives at the points xold.  if natural *
!     *         spline is fitted y2(1)=0. and y2(nold)=0. must be      *
!     *         specified.                                             *
!     *  nnew - number of values of the function to be calculated.     *
!     *  xnew - locations of the points at which the values of the     *
!     *         function are calculated.  xnew(k) must be ge xold(1)   *
!     *         and le xold(nold).                                     *
!     *  ynew - the values of the function to be calculated.           *
!     *  p, q - auxiliary vectors of the length nold-2.                *
!     *                                                                *
!     ******************************************************************
!-----------------------------------------------------------------------
!
      integer,intent(in) :: jtb,nnew,nold
!
      real,dimension(jtb),intent(in) :: xnew,xold,yold
      real,dimension(jtb),intent(inout) :: p,q,y2
      real,dimension(jtb),intent(out) :: ynew
!
!-----------------------------------------------------------------------
!
                        endsubroutine spline    
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine psitbl
!     ******************************************************************
!     *                                                                *
!     *               surface layer integral functions                 *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real,parameter:: &
 eps=0.000001
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        endsubroutine psitbl
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      FUNCTION TIMEF()
!
#ifdef IBM
      REAL*8 TIMEF,rtc
!
      TIMEF=rtc()
#else
      REAL*8 TIMEF
      INTEGER(kind=KINT) :: IC,IR
!
      CALL SYSTEM_CLOCK(count=IC,count_rate=IR)
      TIMEF=REAL(IC)/REAL(IR)*1000.
#endif

!
!
!-----------------------------------------------------------------------
!
      END FUNCTION TIMEF
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine NMMB_Finalize
        integer irc
      end subroutine NMMB_Finalize
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
                       end module module_control
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
