
                        module module_control
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
use module_include
use module_dm_parallel,only : ids,ide,jds,jde &
                             ,ims,ime,jms,jme &
                             ,its,ite,jts,jte &
                             ,lm &
                             ,mype_share,npes,num_pts_max &
                             ,mpi_comm_comp &
                             ,dstrb &
                             ,idstrb &
                             ,dstrb_soil
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
integer(kind=kint),parameter :: &
 itb=201 &                   ! convection tables, dimension 1
,jtb=601 &                   ! convection tables, dimension 2
,kexm=10001 &                ! size of exponentiation table
,kztm=10001                  ! size of surface layer stability function table

real(kind=kfpt),parameter :: &
 ph=105000. &                ! upper bound of pressure range
,thh=350. &                  ! upper bound of potential temperature range
,thl=200.                    ! upper bound of potential temperature range

integer(kind=kint):: &
 kexm2 &                     ! internal points of exponentiation table
,kztm2                       ! internal pts. of the stability function table

integer(kind=kint),private :: &
 mype

real(kind=kfpt) :: &
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
,ztmin2                      ! lower bound for z/L for land stab. functions

real(kind=kfpt),dimension(1:itb):: &
 sthe &                      ! range for equivalent potential temperature
,the0                        ! base for equivalent potential temperature           

real(kind=kfpt),dimension(1:jtb):: &
 qs0 &                       ! base for saturation specific humidity
,sqs                         ! range for saturation specific humidity

real(kind=kfpt),dimension(1:kexm):: &
 expf                        ! exponentiation table

real(kind=kfpt),dimension(1:kztm):: &
 psih1 &                     ! sea heat stability function
,psim1 &                     ! sea momentum stability function
,psih2 &                     ! land heat stability function
,psim2                       ! land momentum stability function

real(kind=kfpt),dimension(1:itb,1:jtb):: &
 ptbl                        ! saturation pressure table

real(kind=kfpt),dimension(1:jtb,1:itb):: &
 ttbl                        ! temperature table
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---miscellaneous control parameters------------------------------------
!-----------------------------------------------------------------------
character(64):: &
 infile
character(64):: &
 input_data
 
logical(kind=klog):: &
 adv_standard &              ! my task is in standard advec region
,adv_upstream &              ! my task is in upstream advec region
,first &                     ! if true preparation for the first time step
,global &                    ! global domain
,hydro &                     ! if true hydrostatic dynamics
,read_global_sums &          ! read global sums in subroutine adv2 (bit identity)
,readbc &                    ! read regional boundary conditions
,run &                       ! initial data ready, start run
,runbc &                     ! boundary data ready, start run
,write_global_sums           ! write global sums in subroutine adv2 (bit identity)
 
integer(kind=kint):: &
 idtad &                     ! time steps between passive substance advection
,idtadt &                    ! time steps between tracer advection
,ierr &                      ! error code
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

integer(kind=kint),dimension(1:3):: &
 idat(3) &                   ! date of initial data, day, month, year
,idatbc(3)                   ! date of boundary data, day, month, year

integer(kind=kint),dimension(1:15):: & 
 nfftrh &                    ! fft working field, h points
,nfftrw                      ! fft working field, v points
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---grid constants------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
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

logical(kind=klog) :: &                 
 e_bdy &                     ! logical flag for tasks on eastern boundary
,n_bdy &                     ! logical flag for tasks on northern boundary
,s_bdy &                     ! logical flag for tasks on southern boundary
,w_bdy                       ! logical flag for tasks on western boundary

real(kind=kfpt):: &
 bofac &                     ! amplification of diffusion along bndrs.
,ctph0 &                     ! cos(tph0)
,ddmpv &                     ! divergence damping coefficient, y direction
,dlm &                       ! grid size, delta lambda, radians
,dlmd &                      ! grid size, delta lambda, degrees
,dph &                       ! grid size, delta phi, radians
,dphd &                      ! grid size, delta phi, degrees
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

!real(kind=kfpt),allocatable,dimension(:):: &      !lm
! dsg1 &                      ! thicnesses of sigma layers in pressure range
!,dsg2 &                      ! thicnesses of sigma layers in sigma range
!,pdsg1 &                     ! thicnesses of pressure layers in press. range
!,psgml1 &                    ! pressure at midlayers in pressure range
!,sgml1 &                     ! sigma at midlayers in pressure range
!,sgml2                       ! sigma at midlayers in sigma range

!real(kind=kfpt),allocatable,dimension(:):: &      !lm+1
!,psg1 &                      ! pressure at interfaces in pressure range
!,sg1 &                       ! sigma at interfaces in pressure range
!,sg2 &                       ! sigma at interfaces in sigma range
!,sgm                         ! sigma at interfaces

integer(kind=kint),allocatable,dimension(:):: &      !jm
 khfilt &                    ! polar filter, truncation wave #, h points
,kvfilt &                    ! polar filter, truncation wave #, v points
,nhsmud                      ! polar smoother for unfiltered variables

!real(kind=kfpt),allocatable,dimension(:):: &      !jm
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

real(kind=kfpt),allocatable,dimension(:):: &      !2*(im-3)
 wfftrh &                    ! fft working field, h points
,wfftrw                      ! fft working field, v points

real(kind=kfpt),allocatable,dimension(:,:):: &      !im,jm
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
integer(kind=kint),allocatable,dimension(:,:):: &      !im,jm
 insoil &                    ! soil type
,inveg                       ! vegetation type

real(kind=kfpt),allocatable,dimension(:,:):: &      !im,jm
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
integer(kind=kint),allocatable,dimension(:,:):: &      !im,jm
 kmp                         ! # of active soil layers

real(kind=kfpt),allocatable,dimension(:,:):: &      !im,jm
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

real(kind=kfpt),allocatable,dimension(:,:,:):: wg   ! soil moisture(nwetc,im,jm)
real(kind=kfpt),allocatable,dimension(:,:,:):: tg   ! soil temperature (kmsc,im,jm)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---regional boundary conditions----------------------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),allocatable,dimension(:,:,:):: &      !im,lnsh,2
 pdbn &                      ! pressure difference at northern boundary
,pdbs                        ! pressure difference at southern boundary

real(kind=kfpt),allocatable,dimension(:,:,:):: &      !lnsh,jm,2
 pdbe &                      ! pressure difference at eastern boundary
,pdbw                        ! pressure difference at western boundary

real(kind=kfpt),allocatable,dimension(:,:,:,:):: &    !im,lnsh,lm,2
 tbn &                       ! temperature at northern boundary
,tbs &                       ! temperature at southern boundary
,qbn &                       ! specific humidity at northern boundary
,qbs &                       ! specific humidity at southern boundary
,wbn &                       ! condensate at northern boundary
,wbs                         ! condensate at southern boundary

real(kind=kfpt),allocatable,dimension(:,:,:,:):: &    !lnsh,jm,lm,2
 tbe &                       ! temperature at eastern boundary
,tbw &                       ! temperature at western boundary
,qbe &                       ! specific humidity at eastern boundary
,qbw &                       ! specific humidity at western boundary
,wbe &                       ! condensate at eastern boundary
,wbw                         ! condensate at western boundary

real(kind=kfpt),allocatable,dimension(:,:,:,:):: &    !im,lnsv,lm,2
 ubn &                       ! u wind component at northern boundary
,ubs &                       ! u wind component at southern boundary
,vbn &                       ! v wind component at northern boundary
,vbs                         ! v wind component at southern boundary

real(kind=kfpt),allocatable,dimension(:,:,:,:):: &    !lnsv,im,lm,2
 ube &                       ! u wind component at eastern boundary
,ubw &                       ! u wind component at western boundary
,vbe &                       ! v wind component at eastern boundary
,vbw                         ! v wind component at western boundary
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---atmospheric variables, hydrostatic----------------------------------
!-----------------------------------------------------------------------
!real(kind=kfpt),allocatable,dimension(:,:):: &      !im,jm
! pd &                        ! pressure difference, sigma range
!,pdo &                       ! previous pressure difference, sigma range
!,psdt                        ! hydrostatic surface pressure tendency

real(kind=kfpt),allocatable,dimension(:,:,:):: &      !im,jm,lm-1
 psgdt                       ! vertical mass flux

!real(kind=kfpt),allocatable,dimension(:,:,:):: &      !im,jm,lm
! cw &                        ! condensate
!,omgalf &                    ! omega-alpha
!,q &                         ! specific humidity
!,q2 &                        ! 2tke
!,rrw &                       ! condensate species mixing ratio
!,t &                         ! temperature
!,tp &                        ! previous temperature
!,u &                         ! u wind component
!,up &                        ! previous u wind component
!,v &                         ! v wind component
!,vp                          ! previous v wind component

!real(kind=kfpt),allocatable,dimension(:,:,:):: &      !im,jm,lm+1
! pint                        ! nonhydrostatic pressure at interfaces
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---atmospheric variables, nonhydrostatic-------------------------------
!-----------------------------------------------------------------------
!real(kind=kfpt),allocatable,dimension(:,:,:):: &      !im,jm,lm
! dwdt &                      ! vertical acceleration, correction factor
!,pdwdt &                     ! previous correction factor
!,rtop &                      ! RT/p, specific volume
!,w &                         ! vertical velocity
!,z                           ! height at midlayers
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---physics variables---------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint),allocatable,dimension(:,:):: &      !im,jm
 lbot &                      ! convective cloud base level
,lpbl &                      ! pbl top level
,ltop                        ! convective cloud top level
!-----------------------------------------------------------------------
!---scratch area for passing temporary arguments------------------------
!-----------------------------------------------------------------------
real(kind=kfpt),allocatable,dimension(:,:,:):: &      !im,jm,lm
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
      (global,idtad,idtadt,secdif &
      ,smag2,smag4,codamp,wcor &
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
integer(kind=kint),intent(in) :: &
 lnsh

integer(kind=kint),intent(out) :: &
 idtad &
,idtadt &
,lnsad &
,nboco

real(kind=kfpt),intent(in) :: &
 codamp &    ! divergence damping coefficient
,pt &        ! Pressure at top of domain (Pa)
,smag2 &     ! Smagorinsky coefficient for 2nd order diffusion 
,smag4 &     ! Smagorinsky coefficient for 4th order diffusion
,tlm0d &
,tph0d &
,wcor

logical(kind=klog),intent(in) :: &
 global
 
real(kind=kfpt),intent(out) :: &
 ddmpv &
,dlmd &
,dphd &
,dyh &
,dyv &
,ef4t &
,rdyh &
,rdyv &
,tboco
 
real(kind=kfpt),dimension(jds:jde),intent(out) :: &
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
 
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out) :: &
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
logical(kind=klog),intent(in) :: &
 secdif
!
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
!
integer(kind=kint):: &
 i &                         ! index in x direction
,i_hi &                      ! max i loop limit (cannot be > ide)
,i_lo &                      ! min i loop limit (cannot be < ids)
,j &                         ! index in y direction
,j_hi &                      ! max j loop limit (cannot be > jde)
,j_lo                        ! min j loop limit (cannot be < jds)

real(kind=kfpt):: &
 acdt &                      ! diffusion coefficient parameter
,alm &                       ! lambda
,anum &                      ! numerator
,aph &                       ! phi
,ave &                       ! average
,cddamp &                    ! divergence damping factor
,coac &                      ! nonlinear diffusion coefficient
,ctlm &                      ! temporary
,denom &                     ! denominator
,relm &                      ! temporary
,fpole &                     ! factor to modify polar area
,sph &                       ! temporary
,stlm &                      ! temporary
,stph &                      ! temporary
,ctph &                      ! temporary
,tlm &                       ! rotated lambda
,tlm_base &                  ! temporary
,tph &                       ! rotated phi
,tph_base &                  ! temporary
,tpv                         ! rotated phi, v points

real(kind=kfpt),dimension(jds:jde):: &
 hdacxj &
,hdacyj &
,hdacvxj &
,hdacvyj

!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
 1000 format(100a4)
!
      if(mype==0)then
        write(0,*)' ids=',ids,' ide=',ide,' jds=',jds,' jde=',jde
      endif
!-----------------------------------------------------------------------
!***  Because subdomains that lie along the global domain boundary
!***  may have haloes that extend beyond the global limits, create
!***  limits here that keep loops from reaching beyond those
!***  global limits.
!-----------------------------------------------------------------------
!
      i_lo=max(ims,ids)
      i_hi=min(ime,ide)
      j_lo=max(jms,jds)
      j_hi=min(jme,jde)
!
!-----------------------------------------------------------------------
      if(global) then
!-----------------------------------------------------------------------
!---global branch-------------------------------------------------------
!-----------------------------------------------------------------------
!
        icycle=ide
!
        dlmd=-wbd*2./real(ide-3)
        dphd=-sbd*2./real(jde-3)
!
!!!     pt=1000.     !<-- Read from input file
!!!     ptsgm=22000. 
!
        lnsbc=lnsh
        bofac=0.
!
        lnsad=0
!-----------------------------------------------------------------------
      else
!-----------------------------------------------------------------------
!---regional branch-----------------------------------------------------
!-----------------------------------------------------------------------
! 
        icycle=ide
	write(0,*) 'wbd, sbd: ', wbd, sbd
        dlmd=-wbd*2./real(ide-1)
        dphd=-sbd*2./real(jde-1)
	write(0,*) 'dlmd, dphd: ', dlmd, dphd
!
!!!     pt=5000.     !<-- Read from input file
!!!     ptsgm=42000.
!
        lnsbc=lnsh
        bofac=4.
!
        lnsad=lnsh+2
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
      adv_upstream=.false.
      if(jts<jds+1+lnsad.or.jte>jde-1-lnsad.or. &
         its<ids+1+lnsad.or.ite>ide-1-lnsad)then
        adv_upstream=.true.
      endif
!
      adv_standard=.false.
      if(jte>=jds+1+lnsad.and.jts<=jde-1-lnsad.and. &
         ite>=ids+1+lnsad.and.its<=ide-1-lnsad)then
        adv_standard=.true.
      endif
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      do j=jds,jde
        dxh(j)=0.
        rdxh(j)=0.
        dare(j)=0.
        rare(j)=0.
        wpdar(j)=0.
        fah(j)=0.
        fcp(j)=0.
        fdiv(j)=0.
        dxv(j)=0.
        rdxv(j)=0.
        curv(j)=0.
        fad(j)=0.
        ddmpu(j)=0.
        hdacxj(j)=0.
        hdacyj(j)=0.
        hdacvxj(j)=0.
        hdacvyj(j)=0.
      enddo
!-----------------------------------------------------------------------
!!!   input_data='/gpfstmp/wx22tb/umo_data/'
!
!!!   idtad=4
!!!   idtadt=4
!!!   nphs=8
!!!   ncnvc=nphs
!!!   nrads=nint(3600./dt)
!!!   nradl=nint(3600./dt*3.)
!!!   nprec=nint(24.*3600./dt)
      dtq2=nphs*dt
!
      if(.not.global) then
        ihrbc=0
        write(infile,'(a,i3.3)')'boco.',ihrbc
        nbc=18
        open(unit=nbc,file=infile,status='old',form='unformatted')
        read (nbc) runbc,idatbc,ihrstbc,tboco
	write(0,*) 'runbc: ', runbc
	write(0,*) 'idatbc: ', idatbc
	write(0,*) 'ihrstbc: ', ihrstbc
	write(0,*) 'tboco: ', tboco
        rewind nbc
        close(unit=nbc)
        write(0,*)'*** Read tboco in consts from ',infile
        nboco=nint(tboco/dt)
      endif
!-----------------------------------------------------------------------
!---to be read from a namelist in the future----------------------------
!-----------------------------------------------------------------------
      if(secdif) then
        coac=smag2*smag2 !second order
      else
        coac=smag4*smag4 !fourth order
      endif
!
!     if(global) then
!       coac=0.01 !global
!     else
!       coac=0.16 !regional
!     endif
!
!-----------------------------------------------------------------------
!--------------derived geometrical constants----------------------------
!-----------------------------------------------------------------------
      tph0=tph0d*dtr                                                    
      wb=wbd*dtr
      sb=sbd*dtr
      dlm=dlmd*dtr                                                      
      dph=dphd*dtr                                                      
      rdlm=1./dlm                                                       
      rdph=1./dph                                                       
!                                                                       
      stph0=sin(tph0)                                                   
      ctph0=cos(tph0)                                                   
!-----------------------------------------------------------------------
!---derived horizontal grid constants-----------------------------------
!-----------------------------------------------------------------------
      dyh=a*dph
      dyv=a*dph
      rdyh=1./dyh
      rdyv=1./dyv
!-----------------------------------------------------------------------
!---derived vertical grid constants-------------------------------------
!-----------------------------------------------------------------------
      ef4t=0.5*dt/cp
      f4d=-0.5*dt
!-----------------------------------------------------------------------
!---derived physical constants------------------------------------------
!-----------------------------------------------------------------------
      acdt=coac*dt
      cddamp=codamp*dt
!-----------------------------------------------------------------------
!***  Each task identifies whether or not it is adjacent to a boundary
!***  of the full domain.
!-----------------------------------------------------------------------
!
      w_bdy=(its==ids)  ! This task is on the western boundary
      e_bdy=(ite==ide)  ! This task is on the eastern boundary
      s_bdy=(jts==jds)  ! This task is on the southern boundary
      n_bdy=(jte==jde)  ! This task is on the northern boundary
!
!-----------------------------------------------------------------------
!
      global_regional_setup: if(global) then
!
!-----------------------------------------------------------------------
!---global branch-------------------------------------------------------
!-----------------------------------------------------------------------
        tph=sb-dph
        fpole=1.0                                                                      
!-----------------------------------------------------------------------
!----south pole---------------------------------------------------------
!-----------------------------------------------------------------------
        tph=tph+dph
        tpv=tph+dph*0.5                                               
        dxh(jds+1)=0.
        dxv(jds+1)=a*dlm*cos(tpv)
        ddmpv=cddamp*dyv/(2.*dyv)                                                         
!
        rdxh(jds+1)=0.
        rdxv(jds+1)=1./dxv(jds+1)
        dare(jds+1)=dxv(jds+1)*dyv*0.5*fpole !for ghost area
        rare(jds+1)=1./dare(jds+1)
        wpdar(jds+1)=-wcor*(dyh)**2 &
                     /(dt*dxv(jds+1)*dyv*0.5*fpole)/100000.00
        curv(jds+1)=tan(tpv)/a
        fad(jds+1)=-dt/(3.*dxv(jds+1)*dyv*2.*2.)
        fah(jds+1)=-dt/(3.*dxv(jds+1)*dyv*0.5*fpole)
        fcp(jds+1)= dt/(3.*dxv(jds+1)*dyv*0.5*cp*fpole)
        fdiv(jds+1)=2./(3.*dxv(jds+1)*dyv*0.5*fpole)
!
        hdacxj(jds+1)=0.
        hdacyj(jds+1)=acdt*dyh**2 &
                  /(4.*sqrt(2.)*dxv(jds+1)*dyv*0.5*fpole)
        hdacvxj(jds+1)=acdt*dyv**2 &
                  /(4.*sqrt(2.)*dxv(jds+1)*dyv)
        hdacvyj(jds+1)=acdt*dyv**2 &
                  /(4.*sqrt(2.)*dxv(jds+1)*dyv)
!
!        ddmpu(jds+1)=cddamp*dxv(jds+1)/(2.*dxv(jds+1))
        ddmpu(jds+1)=cddamp*dyv/(2.*dxv(jds+1))
!-----------------------------------------------------------------------
!-------------between the poles----------------------------------------- 
!-----------------------------------------------------------------------
        do j=jds+2,jde-2
          tph=sb+(j-jds-1)*dph
          tpv=tph+dph*0.5
          dxh(j)=a*dlm*cos(tph)
          dxv(j)=a*dlm*cos(tpv)
          rdxh(j)=1./dxh(j)
!     if(abs(dxv(j))<1.e-5)then
!       write(0,*)' in CONSTS j=',j,' dxv=',dxv(j)
!        dxv(j)=1.0*dxv(j)
!     endif
          rdxv(j)=1./dxv(j)
          dare(j)=dxh(j)*dyh
          rare(j)=1./dare(j)
          wpdar(j)=-wcor*(dyh)**2 &
                  /(dt*dxh(j)*dyh)/100000.00
          curv(j)=tan(tpv)/a
          fad(j)=-dt/(3.*dxv(j)*dyv*2.*2.)
          fah(j)=-dt/(3.*dxh(j)*dyh)
          fcp(j)= dt/(3.*dxh(j)*dyh*cp)
          fdiv(j)=2./(3.*dxh(j)*dyh)
!
          hdacxj(j)= acdt*dyh*max(dxh(j),dyh)/(4.*sqrt(2.)*dxh(j)*dyh)
          hdacyj(j)= acdt*dyh*max(dxh(j),dyh)/(4.*sqrt(2.)*dxh(j)*dyh)
          hdacvxj(j)=acdt*dyv*max(dxv(j),dyv)/(4.*sqrt(2.)*dxv(j)*dyv)
          hdacvyj(j)=acdt*dyv*max(dxv(j),dyv)/(4.*sqrt(2.)*dxv(j)*dyv)
!
!          ddmpu(j)=cddamp*dxv(j)/(2.*dxv(j))
          ddmpu(j)=cddamp*dyv/(2.*dxv(j))
        enddo
!-----------------------------------------------------------------------
!-------------ghost line beyond the south pole---------------------------- 
!-----------------------------------------------------------------------
        dxh(jds)=dxh(jds+2)
        dxv(jds)=dxv(jds+1)
        rdxh(jds)=rdxh(jds+2)
        rdxv(jds)=rdxv(jds+1)
        dare(jds)=dare(jds+2)
        rare(jds)=rare(jds+2)
        wpdar(jds)=wpdar(jds+2)
        curv(jds)=curv(jds+1)
        fad(jds)=fad(jds+1)
        fah(jds)=fah(jds+2)
        fcp(jds)=fcp(jds+2)
        fdiv(jds)=fdiv(jds+2)
        hdacxj(jds)=hdacxj(jds+2)
        hdacyj(jds)=hdacyj(jds+2)
        hdacvxj(jds)=hdacvxj(jds+1)
        hdacvyj(jds)=hdacvyj(jds+1)
        ddmpu(jds)=ddmpu(jds+1)
!-----------------------------------------------------------------------
!-------------north pole------------------------------------------------
!-----------------------------------------------------------------------
        tph=tph+dph
        tpv=tph+dph*0.5
        dxh(jde-1)=0.
        rdxh(jde-1)=0.
        dare(jde-1)=dxv(jde-2)*dyv*0.5*fpole
        rare(jde-1)=1./dare(jde-1)
        wpdar(jde-1)=-wcor*(dyh)**2 &
                     /(dt*dxv(jde-2)*dyv*0.5*fpole)/100000.00
        fah(jde-1)=-dt/(3.*dxv(jde-2)*dyv*0.5*fpole)
        fcp(jde-1)= dt/(3.*dxv(jde-2)*dyv*0.5*fpole*cp)
        fdiv(jde-1)=2./(3.*dxv(jde-2)*dyv*0.5*fpole)
        hdacxj(jde-1)=0.
        hdacyj(jde-1)=acdt*dyh**2 &
                    /(4.*sqrt(2.)*dxv(jde-2)*dyv*0.5*fpole)
!-----------------------------------------------------------------------
!-------------ghost line beyond north pole------------------------------
!-----------------------------------------------------------------------
        dxh(jde)=dxh(jde-2)
        dxv(jde-1)=dxv(jde-2)
        dxv(jde)=dxv(jde-2)
        rdxh(jde)=rdxh(jde-2)
        rdxv(jde-1)=rdxv(jde-2)
        rdxv(jde)=rdxv(jde-2)
        dare(jde)=dare(jde-2)
        rare(jde)=rare(jde-2)
        wpdar(jde)=wpdar(jde-2)
        curv(jde-1)=curv(jde-2)
        curv(jde-1)=curv(jde-2)
        fad(jde-1)=fad(jde-2)
        fad(jde)=fad(jde-2)
        fah(jde)=fah(jde-2)
        fcp(jde)=fcp(jde-2)
        fdiv(jde)=fdiv(jde-2)
        hdacxj(jde)=hdacxj(jde-2)
        hdacyj(jde)=hdacyj(jde-2)
        hdacvxj(jde-1)=hdacvxj(jde-2)
        hdacvyj(jde-1)=hdacvyj(jde-2)
        hdacvxj(jde)=hdacvxj(jde-2)
        hdacvyj(jde)=hdacvyj(jde-2)
        ddmpu(jde-1)=ddmpu(jde-2)
        ddmpu(jde)=ddmpu(jde-2)
!-----------------------------------------------------------------------
!-------------averaging over height latitudes for accuracy-------------- 
!-----------------------------------------------------------------------
        do j=jds,jde/2
          ave=(dxh(j)+dxh(jde+1-j))*0.5
          dxh(j)=ave
          dxh(jde+1-j)=ave
          ave=(rdxh(j)+rdxh(jde+1-j))*0.5
          rdxh(j)=ave
          rdxh(jde+1-j)=ave
          ave=(dare(j)+dare(jde+1-j))*0.5
          dare(j)=ave
          dare(jde+1-j)=ave
          ave=(rare(j)+rare(jde+1-j))*0.5
          rare(j)=ave
          rare(jde+1-j)=ave
          ave=(wpdar(j)+wpdar(jde+1-j))*0.5
          wpdar(j)=ave
          wpdar(jde+1-j)=ave
          ave=(fah(j)+fah(jde+1-j))*0.5
          fah(j)=ave
          fah(jde+1-j)=ave
          ave=(fcp(j)+fcp(jde+1-j))*0.5
          fcp(j)=ave
          fcp(jde+1-j)=ave
          ave=(fdiv(j)+fdiv(jde+1-j))*0.5
          fdiv(j)=ave
          fdiv(jde+1-j)=ave
          ave=(hdacxj(j)+hdacxj(jde+1-j))*0.5
          hdacxj(j)=ave
          hdacxj(jde+1-j)=ave
          ave=(hdacyj(j)+hdacyj(jde+1-j))*0.5
          hdacyj(j)=ave
          hdacyj(jde+1-j)=ave
        enddo
!-----------------------------------------------------------------------
!-------------averaging over wind latitudes for accuracy--------------- 
!-----------------------------------------------------------------------
        do j=jds,(jde-1)/2
          ave=(dxv(j)+dxv(jde-j))*0.5
          dxv(j)=ave
          dxv(jde-j)=ave
          ave=(rdxv(j)+rdxv(jde-j))*0.5
          rdxv(j)=ave
          rdxv(jde-j)=ave
          ave=(fad(j)+fad(jde-j))*0.5
          fad(j)=ave
          fad(jde-j)=ave
          ave=(hdacvxj(j)+hdacvxj(jde-j))*0.5
          hdacvxj(j)=ave
          hdacvxj(jde-j)=ave
          ave=(hdacvyj(j)+hdacvyj(jde-j))*0.5
          hdacvyj(j)=ave
          hdacvyj(jde-j)=ave
          ave=(ddmpu(j)+ddmpu(jde-j))*0.5
          ddmpu(j)=ave
          ddmpu(jde-j)=ave
        enddo
!-----------------------------------------------------------------------
!-------------diagonal distances at v points----------------------------
!-----------------------------------------------------------------------
        if(s_bdy)then
          ddv(jds+1)=dyv
          ddv(jds)=ddv(jds+1)
          rddv(jds+1)=1./ddv(jds+1)
          rddv(jds)=rddv(jds+1)
        endif
!
        do j=max(jms,jds+2),min(jme,jde-3)
          ddv(j)=sqrt(dxv(j)**2+dyv**2)
          rddv(j)=1./ddv(j)
        enddo
!
        if(n_bdy)then
          ddv(jde-2)=dyv
          ddv(jde-1)=ddv(jde-2)
          rddv(jde-2)=1./ddv(jde-2)
          rddv(jde-1)=rddv(jde-2)
        endif
!-----------------------------------------------------------------------
!-------------defining the diffusion coefficient inside the domain------
!-----------------------------------------------------------------------
!       do j=jms,jme
!         do i=ims,ime
        do j=j_lo,j_hi
          do i=i_lo,i_hi
            hdacx(i,j)=hdacxj(j)
            hdacy(i,j)=hdacyj(j)
            hdacvx(i,j)=hdacvxj(j)
            hdacvy(i,j)=hdacvyj(j)
          enddo
        enddo
!-----------------------------------------------------------------------
!-------------coriolis parameter in tll system--------------------------
!-----------------------------------------------------------------------
        tph_base=sb-dph-dph*.5
!       do j=jms,jme
        do j=j_lo,j_hi
          tph=tph_base+(j-jds+1)*dph
          stph=sin(tph)
          ctph=cos(tph)
!
          tlm_base=wb-dlm-dlm*0.5
!         do i=ims,ime
          do i=i_lo,i_hi
            tlm=tlm_base+(i-ids+1)*dlm
            f(i,j)=twom*(ctph0*stph+stph0*ctph*cos(tlm))
          enddo
        enddo
!-----------------------------------------------------------------------
!--------------latitudes and longitudes of h points in radians----------     
!-----------------------------------------------------------------------
!       tph_base=sb-dph-dph
!       do j=j_lo,j_hi
!         tlm_base=wb-dlm-dlm
!         tph=tph_base+(j-jds+1)*dph
!         stph=sin(tph)
!         ctph=cos(tph)                                                                       
!         do i=i_lo,i_hi
!           tlm=tlm_base+(i-ids+1)*dlm
!           stlm=sin(tlm)
!           ctlm=cos(tlm)
!           sph=ctph0*stph+stph0*ctph*ctlm
!           aph=asin(sph)
!           glat(i,j)=aph
!           anum=ctph*stlm
!           denom=(ctlm*ctph-stph0*sph)/ctph0
!           relm=atan2(anum,denom)
!           alm=relm+tlm0d*dtr
!           if(alm> pi) alm=alm-pi-pi
!           if(alm< -pi) alm=alm+pi+pi
!           glon(i,j)=alm
!         enddo
!       enddo
        tph_base=sb-dph-dph
        do j=j_lo,j_hi
          tlm_base=wb-dlm-dlm
          aph=tph_base+(j-jds+1)*dph
          do i=i_lo,i_hi
            alm=tlm_base+(i-ids+1)*dlm
            if(alm> pi) alm=alm-pi-pi
            if(alm<-pi) alm=alm+pi+pi
            glat(i,j)=aph
            glon(i,j)=alm
          enddo
        enddo
!
!***  Repeat the preceding loop for the offset V points
!
        tph_base=sb-dph-0.5*dph
        do j=j_lo,j_hi
          tlm_base=wb-dlm-0.5*dlm
          aph=tph_base+(j-jds+1)*dph
          do i=i_lo,i_hi
            alm=tlm_base+(i-ids+1)*dlm
            if(alm> pi) alm=alm-pi-pi
            if(alm<-pi) alm=alm+pi+pi
            vlat(i,j)=aph
            vlon(i,j)=alm
          enddo
        enddo
!-----------------------------------------------------------------------
      else !regional
!-----------------------------------------------------------------------
!
!---regional branch-----------------------------------------------------
!
!-----------------------------------------------------------------------
!-------------between the poles----------------------------------------- 
!-----------------------------------------------------------------------
!
        ddmpv=cddamp*dyv/(2.*dyv)
!
!!!!!   do j=jms,jme
        do j=jds,jde
          tph=sb+(j-jds)*dph
          tpv=tph+dph*0.5
          dxh(j)=a*dlm*cos(tph)
          dxv(j)=a*dlm*cos(tpv)
          rdxh(j)=1./dxh(j)
          rdxv(j)=1./dxv(j)
          dare(j)=dxh(j)*dyh
          rare(j)=1./dare(j)
          wpdar(j)=-wcor*(dyh)**2 &
                  /(dt*dxh(j)*dyh)/100000.00
          curv(j)=tan(tpv)/a
          fad(j)=-dt/(3.*dxv(j)*dyv*2.*2.)
          fah(j)=-dt/(3.*dxh(j)*dyh)
          fcp(j)= dt/(3.*dxh(j)*dyh*cp)
          fdiv(j)=2./(3.*dxh(j)*dyh)
!
          hdacxj(j)= acdt*dyh*max(dxh(j),dyh)/(4.*sqrt(2.)*dxh(j)*dyh)
          hdacyj(j)= acdt*dyh*max(dxh(j),dyh)/(4.*sqrt(2.)*dxh(j)*dyh)
          hdacvxj(j)=acdt*dyv*max(dxv(j),dyv)/(4.*sqrt(2.)*dxv(j)*dyv)
          hdacvyj(j)=acdt*dyv*max(dxv(j),dyv)/(4.*sqrt(2.)*dxv(j)*dyv)
!
!          ddmpu(j)=cddamp*dxv(j)/(2.*dxv(j))
          ddmpu(j)=cddamp*dyv/(2.*dxv(j))
        enddo
!-----------------------------------------------------------------------
!-------------diagonal distances at v points----------------------------
!-----------------------------------------------------------------------
!!!!    do j=jms,jme
        do j=jds,jde
          ddv(j)=sqrt(dxv(j)**2+dyv**2)
          rddv(j)=1./ddv(j)
        enddo
!-----------------------------------------------------------------------
!---defining the diffusion coefficient inside the domain----------------
!-----------------------------------------------------------------------
!       do j=jms,jme
!         do i=ims,ime
        do j=j_lo,j_hi
          do i=i_lo,i_hi
            hdacx(i,j)=hdacxj(j)
            hdacy(i,j)=hdacyj(j)
            hdacvx(i,j)=hdacvxj(j)
            hdacvy(i,j)=hdacvyj(j)
          enddo
        enddo
!-----------------------------------------------------------------------
!---enhancing diffusion along boundaries--------------------------------
!-----------------------------------------------------------------------
        if(lnsbc>=2)then
!
          if(s_bdy)then
            do j=jts+1,jts-1+lnsbc
              do i=max(its,ids+1),min(ite,ide-1)
                hdacx(i,j)=hdacxj(j)*bofac
                hdacy(i,j)=hdacyj(j)*bofac
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do j=jte+1-lnsbc,jte-1
              do i=max(its,ids+1),min(ite,ide-1)
                hdacx(i,j)=hdacxj(j)*bofac
                hdacy(i,j)=hdacyj(j)*bofac
              enddo
            enddo
          endif
!
        endif
!
        if(w_bdy)then
          do j=max(jts,jds+lnsbc),min(jte,jde-lnsbc)
            do i=its+1,its-1+lnsbc
              hdacx(i,j)=hdacxj(j)*bofac
              hdacy(i,j)=hdacyj(j)*bofac
            enddo
          enddo
        endif
!
        if(e_bdy)then
          do j=max(jts,jds+lnsbc),min(jte,jde-lnsbc)
            do i=ite+1-lnsbc,ite-1
              hdacx(i,j)=hdacxj(j)*bofac
              hdacy(i,j)=hdacyj(j)*bofac
            enddo
          enddo
        endif
!
!-----------------------------------------------------------------------
        if(s_bdy)then
          do j=jts+1,jts-1+lnsbc
            do i=max(its,ids+1),min(ite,ide-2)
              hdacvx(i,j)=hdacvxj(j)*bofac
              hdacvy(i,j)=hdacvyj(j)*bofac
            enddo
          enddo
        endif
!
        if(n_bdy)then
          do j=jte-lnsbc,jte-2
            do i=max(its,ids+1),min(ite,ide-2)
              hdacvx(i,j)=hdacvxj(j)*bofac
              hdacvy(i,j)=hdacvyj(j)*bofac
            enddo
          enddo
        endif
!
        if(w_bdy)then
          do j=max(jts,jds+lnsbc),min(jte,jde-1-lnsbc)
            do i=its+1,its-1+lnsbc
              hdacvx(i,j)=hdacvxj(j)*bofac
              hdacvy(i,j)=hdacvyj(j)*bofac
            enddo
          enddo
        endif
!
        if(e_bdy)then
          do j=max(jts,jds+lnsbc),min(jte,jde-1-lnsbc)
            do i=ite-lnsbc,ite-2
              hdacvx(i,j)=hdacvxj(j)*bofac
              hdacvy(i,j)=hdacvyj(j)*bofac
            enddo
          enddo
        endif
!
!-----------------------------------------------------------------------
!-------------coriolis parameter in tll system--------------------------
!-----------------------------------------------------------------------
      
        tph_base=sb-dph*.5
!
        do j=jts,jte
          tph=tph_base+(j-jds+1)*dph
          stph=sin(tph)
          ctph=cos(tph)
!
          tlm_base=wb-dlm*0.5
          do i=its,ite
            tlm=tlm_base+(i-ids+1)*dlm
            f(i,j)=twom*(ctph0*stph+stph0*ctph*cos(tlm))
          enddo
        enddo
!-----------------------------------------------------------------------
!--------------latitudes and longitudes of h points in radians----------     
!-----------------------------------------------------------------------
        tph_base=sb-dph
	write(0,*) 'sb,dph: ', sb, dph
	write(0,*) 'wb,dlm: ', wb, dlm
!
        do j=jts,jte
          tph=tph_base+(j-jds+1)*dph
          stph=sin(tph)
          ctph=cos(tph)                                                                       
!
          tlm_base=wb-dlm
          do i=its,ite
            tlm=tlm_base+(i-ids+1)*dlm
            stlm=sin(tlm)
            ctlm=cos(tlm)
            sph=ctph0*stph+stph0*ctph*ctlm
            aph=asin(sph)
            glat(i,j)=aph
            anum=ctph*stlm
            denom=(ctlm*ctph-stph0*sph)/ctph0
            relm=atan2(anum,denom)
            alm=relm+tlm0d*dtr
            if(alm>  pi) alm=alm-pi-pi
            if(alm< -pi) alm=alm+pi+pi
            glon(i,j)=alm
          enddo
        enddo
!
!***  Repeat preceding loop for V lat/lon offset
!
        tph_base=sb-0.5*dph
!
        do j=jts,jte
          tph=tph_base+(j-jds+1)*dph
          stph=sin(tph)
          ctph=cos(tph)
!
          tlm_base=wb-0.5*dlm
          do i=its,ite
            tlm=tlm_base+(i-ids+1)*dlm
            stlm=sin(tlm)
            ctlm=cos(tlm)
            sph=ctph0*stph+stph0*ctph*ctlm
            aph=asin(sph)
            vlat(i,j)=aph
            anum=ctph*stlm
            denom=(ctlm*ctph-stph0*sph)/ctph0
            relm=atan2(anum,denom)
            alm=relm+tlm0d*dtr
            if(alm>  pi) alm=alm-pi-pi
            if(alm< -pi) alm=alm+pi+pi
            vlon(i,j)=alm
          enddo
        enddo
!-----------------------------------------------------------------------
!
      endif global_regional_setup
!
!-----------------------------------------------------------------------
!-------------look-up tables--------------------------------------------
!-----------------------------------------------------------------------
!
      call exptbl
!
!-----------------------------------------------------------------------
      pl=max(pt,1000.)
      call tablep
      call tablet
!-----------------------------------------------------------------------
      call psitbl
!-----------------------------------------------------------------------
!
                        end subroutine consts       
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine init &
      (global &
      ,kss,kse &
      ,pdtop,pt,lpt2 &
      ,sg1,dsg1 &
      ,psg1,pdsg1 &
      ,sg2,dsg2,sgm &
      ,sgml1,psgml1,sgml2 &
      ,fis,sm,sice &
      ,pd,pdo,pint &
      ,u,v,q2 &
      ,t,q,cw &
      ,tp,up,vp &
      ,rrw,dwdt,w &
      ,sp,indx_q,indx_cw,indx_rrw,indx_q2 &
      ,ntsti,ntstm &
      ,ihr,ihrst,idat &
      ,run &
      ,num_water,water &
      ,p_qv,p_qc,p_qr &
      ,p_qi,p_qs,p_qg &
       )
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 kse &
,kss &
,num_water &
,p_qv &
,p_qc &
,p_qr &
,p_qi &
,p_qs &
,p_qg &
,indx_q &
,indx_cw &
,indx_rrw &
,indx_q2

integer(kind=kint),intent(out) :: &
 ihr &
,ihrst &
,lpt2 &
,ntsti &
,ntstm

integer(kind=kint),dimension(3),intent(out) :: &
 idat

real(kind=kfpt),intent(out) :: &
 pdtop &
,pt

real(kind=kfpt),dimension(1:lm),intent(out) :: &
 dsg1 &
,dsg2 &
,pdsg1 &
,psgml1 &
,sgml1 &
,sgml2

real(kind=kfpt),dimension(1:lm+1),intent(out) :: &
 psg1 &
,sg1 &
,sg2 &
,sgm

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out) :: &
 fis &
,pd &
,pdo &
,sice &
,sm  

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out) :: &
 cw &
,dwdt &
,q &
,q2 &
,rrw &
,t &
,tp &
,u &
,up &
,v &
,vp &
,w

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(out) :: &
 pint

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:num_water),intent(out) :: &
 water


real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,kss:kse),intent(out) :: &
 sp 

logical(kind=klog),intent(in) :: &
 global
 
logical(kind=klog),intent(out) :: &
 run
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &                         ! index in x direction
,irtn &
,j &                         ! index in y direction
,k &                         ! index
,ks &                        ! tracer index
,l &                         ! index in p direction
,n &
,nrecs_skip_for_pt

integer(kind=kint),allocatable,dimension(:,:) :: &
 insoil &
,inveg &
,itemp

real(kind=kfpt):: &
 tend,tend_max   

real(kind=kfpt),allocatable,dimension(:,:) :: &
 temp1 

logical(kind=klog) :: opened
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
      allocate(temp1(ids:ide,jds:jde),stat=i)
      allocate(itemp(ids:ide,jds:jde),stat=i)
      allocate(stdh(ims:ime,jms:jme),stat=i)
!
!-----------------------------------------------------------------------
      infile='                                                          '
!-----------------------------------------------------------------------

      select_unit: do n=51,59
        inquire(n,opened=opened)
        if(.not.opened)then
          nfcst=n
          exit select_unit
        endif
      enddo select_unit

!
      infile='main_input_filename'
      open(unit=nfcst,file=infile,status='old',form='unformatted')
!
!-----------------------------------------------------------------------
!***  First we need to extract the value of PT (pressure at top of domain)
!-----------------------------------------------------------------------
!
      if(mype==0)then
        nrecs_skip_for_pt=6+5*lm+21 !<-- For current WPS iniput
!       nrecs_skip_for_pt=6+5*lm+23 !zj 21
!
        do n=1,nrecs_skip_for_pt
          read(nfcst)
        enddo
!
        read(nfcst)pt
        write(0,*)' PT from input file equals ',pt
        rewind nfcst
      endif
      call mpi_bcast(pt,1,mpi_real,0,mpi_comm_comp,irtn)
!
!-----------------------------------------------------------------------
!
      read(nfcst) run,idat,ihrst,ihrend,ntsd
!     read(nfcst)

!     run=.true.

      if(mype==0)then
        write(0,*) 'run, idat,ntsd: ', run, idat, ntsd
      endif
!
!-----------------------------------------------------------------------
!***  Print the time information.
!-----------------------------------------------------------------------
!
      if(mype==0)then
        write(0,*)' Start year =',idat(3)
        write(0,*)' Start month=',idat(2)
        write(0,*)' Start day  =',idat(1)
        write(0,*)' Start hour =',ihrst
        write(0,*)' Timestep   =',dt
        write(0,*)' Steps/hour =',3600./dt
        if(.not.global)write(0,*)' Max fcst hours=',ihrend
      endif
!
!-----------------------------------------------------------------------
!
      read (nfcst) pdtop,lpt2,sgm,sg1,dsg1,sgml1,sg2,dsg2,sgml2
!
      if(mype==0)then
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        fis(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,fis,1,1,1,1,1)
      call halo_exch(fis,1,2,2)
!-----------------------------------------------------------------------
!
      if(mype==0)then
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        stdh(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,stdh,1,1,1,1,1)
      call halo_exch(stdh,1,2,2)
!-----------------------------------------------------------------------
!
      if(mype==0)then
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        sm(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,sm,1,1,1,1,1)
      call halo_exch(sm,1,2,2)
!-----------------------------------------------------------------------
!
      if(mype==0)then
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        pd(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,pd,1,1,1,1,1)
      call halo_exch(pd,1,2,2)
!-----------------------------------------------------------------------
!
      call mpi_barrier(mpi_comm_comp,irtn)
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          u(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,u,1,1,1,lm,l)
      enddo
      call halo_exch(u,lm,2,2)
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          v(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,v,1,1,1,lm,l)
      enddo
      call halo_exch(v,lm,2,2)
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
          write(0,*) 'L, T extremes: ', L, minval(temp1),maxval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          t(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,t,1,1,1,lm,l)
      enddo
      call halo_exch(t,lm,2,2)
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          q(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,q,1,1,1,lm,l)
      enddo
      call halo_exch(q,lm,2,2)
!
      do l=1,lm
      do j=jms,jme
      do i=ims,ime
        water(i,j,l,p_qv)=q(i,j,l)/(1.-q(i,j,l))    ! WRF water array uses mixing ratio for vapor
      enddo
      enddo
      enddo
!
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          cw(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,cw,1,1,1,lm,l)
      enddo
      call halo_exch(cw,lm,2,2)
!
      print*,'*** Read initial conditions in init from ',infile
!
!      close(unit=nfcst)
!-----------------------------------------------------------------------
      ntsti=ntsd+1
!
      tend_max=real(ihrend)
      ntstm_max=nint(tend_max*3600./dt)+1
      tend=real(nhours_fcst)
      ntstm=nint(tend*3600./dt)+1
      if(.not.global)then
         write(0,*)' Max runtime is ',tend_max,' hours'
      endif
      if(mype==0)then
        write(0,*)' Requested runtime is ',tend,' hours'
        write(0,*)' NTSTM=',ntstm
      endif
      if(ntstm>ntstm_max.and..not.global)then
        if(mype==0)then
          write(0,*)' Requested fcst length exceeds maximum'
          write(0,*)' Resetting to maximum'
        endif
        ntstm=min(ntstm,ntstm_max)
      endif
!
      ihr=nint(ntsd*dt/3600.)
!-----------------------------------------------------------------------
      do l=1,lm
        pdsg1(l)=dsg1(l)*pdtop
        psgml1(l)=sgml1(l)*pdtop+pt
!	write(0,*) 'L, pdsg1, psgml1: ', L, pdsg1(L), psgml1(L)
      enddo
!
      do l=1,lm+1
        psg1(l)=sg1(l)*pdtop+pt
      enddo
!-----------------------------------------------------------------------
      do j=jts,jte
        do i=its,ite
          pdo(i,j)=pd(i,j)
        enddo
      enddo
      call halo_exch(pdo,1,2,2)
!
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            up(i,j,l)=u(i,j,l)
            vp(i,j,l)=v(i,j,l)
            tp(i,j,l)=t(i,j,l)
          enddo
        enddo
      enddo
      call halo_exch(tp,lm,up,lm,vp,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            q2(i,j,l)=0.02
            rrw(i,j,l)=0.
            if(i.ge.ide  /2+1- 6.and.i.le.ide  /2+1+ 6.and. &
               j.ge.jde*3/4+1- 6.and.j.le.jde*3/4+1+ 6.) then !global
!               j.ge.jde  /2+1- 6.and.j.le.jde  /2+1+ 6.) then !regional
              rrw(i,j,l)=10.
            endif
            dwdt(i,j,l)=1.
            w(i,j,l)=0.
          enddo
        enddo
      enddo
      call halo_exch(dwdt,lm,2,2)
!
      do j=jts,jte
        do i=its,ite
          pint(i,j,1)=pt
        enddo
      enddo
!
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            pint(i,j,l+1)=pint(i,j,l)+dsg2(l)*pd(i,j)+pdsg1(l)
          enddo
        enddo
      enddo
      call halo_exch(pint,lm+1,2,2)
!
      call halo_exch(q2,lm,rrw,lm,2,2)
      do l=1,lm
        do j=jms,jme
          do i=ims,ime
            sp(i,j,l,indx_q  )=sqrt(max(q  (i,j,l),0.))
            sp(i,j,l,indx_cw )=sqrt(max(cw (i,j,l),0.))
            sp(i,j,l,indx_rrw)=sqrt(max(rrw(i,j,l),0.))
            sp(i,j,l,indx_q2 )=sqrt(max(q2 (i,j,l),0.))
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!---reading surface data------------------------------------------------
!-----------------------------------------------------------------------
!
       if(mype==0)then
         read(nfcst)temp1  ! ALBEDO
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! EPSR
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! MXSNAL
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! TSKIN
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! SST
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! SNO
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! SI
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! SICE
       endif
!      write(0,*) 'maxval(sice): ', maxval(temp1)
       do j=jms,jme
       do i=ims,ime
         sice(i,j)=0.
       enddo
       enddo
       call dstrb(temp1,sice,1,1,1,1,1)
       call halo_exch(sice,1,2,2)
!
!-----------------------------------------------------------------------
!
      if(mype==0)then
        do n=1,13 !<-- For current WPS input
!       do n=1,15 !zj 13
          read(nfcst)
        enddo
      endif
!
      if(mype==0)read(nfcst)pt
      call mpi_bcast(pt,1,mpi_real,0,mpi_comm_comp,irtn)
!
!-----------------------------------------------------------------------
!
      if(mype==0)then
        close(nfcst)
      endif
!
!-----------------------------------------------------------------------
!
      deallocate(temp1)
      deallocate(itemp)
      deallocate(stdh)
      if(mype==0)then
        write(0,*)' EXIT SUBROUTINE INIT pt=',pt
      endif
!-----------------------------------------------------------------------
!
                        end subroutine init
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
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
integer(kind=kint):: &
 k                           ! index

real(kind=kfpt):: &
 x &                         ! argument
,xrng                        ! argument range
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      xmax= 30.
      xmin=-30.
!
      kexm2=kexm-2
      xrng=xmax-xmin
!
      dex=xrng/(kexm-1)
      rdex=1./dex
!--------------function definition loop---------------------------------
      x=xmin-dex
      do k=1,kexm
        x=x+dex
        expf(k)=exp(x)
      enddo
!-----------------------------------------------------------------------
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
!-----------------------------------------------------------------------
real(kind=kfpt):: &
 zjexp

!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 k                           ! index

real(kind=kfpt):: &
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
real(kind=kfpt),parameter:: &
 eps=1.e-10
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 kth &                       ! index
,kp                          ! index

real(kind=kfpt):: &
 ape &                       ! exner function
,dth &                       ! potential temperature step
,dp &                        ! pressure step
,dqs &                       ! saturation specific humidity step
,p &                         ! pressure
,qs0k &                      ! base value for saturation humidity
,sqsk &                      ! saturation spec. humidity range
,th                          ! potential temperature

real(kind=kfpt),dimension(1:itb):: &
 app &                       ! temporary
,aqp &                       ! temporary
,pnew &                      ! new pressures
,pold &                      ! old pressure
,qsnew &                     ! new saturation spec. humidity
,qsold &                     ! old saturation spec. humidity
,y2p                         ! temporary
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      dth=(thh-thl)/real(jtb-1)
      dp=(ph-pl)/real(itb-1)
      rdth=1./dth
!-----------------------------------------------------------------------
      th=thl-dth
      do kth=1,jtb
        th=th+dth
        p=pl-dp
        do kp=1,itb
          p=p+dp
          ape=(100000./p)**cappa
          qsold(kp)=pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          pold(kp)=p
        enddo
!
        qs0k=qsold(1)
        sqsk=qsold(itb)-qsold(1)
        qsold(1)=0.
        qsold(itb)=1.
!
        do kp=2,itb-1
          qsold(kp)=(qsold(kp)-qs0k)/sqsk
!wwwwwwwwwwwwww fix due to 32 bit precision limitation wwwwwwwwwwwwwwwww
          if((qsold(kp)-qsold(kp-1)).lt.eps) then
            qsold(kp)=qsold(kp-1)+eps
          endif
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        enddo
!
        qs0(kth)=qs0k
        sqs(kth)=sqsk
!
        qsnew(1)=0.
        qsnew(itb)=1.
        dqs=1./real(itb-1)
        rdq=1./dqs
!
        do kp=2,itb-1
          qsnew(kp)=qsnew(kp-1)+dqs
        enddo
!
        y2p(1)=0.
        y2p(itb)=0.
!
        call spline(jtb,itb,qsold,pold,y2p,itb,qsnew,pnew,app,aqp)
!
        do kp=1,itb
          ptbl(kp,kth)=pnew(kp)
        enddo
!-----------------------------------------------------------------------
      enddo
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
real(kind=kfpt),parameter:: &
 eps=1.e-10
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 kth &                       ! index
,kp                          ! index

real(kind=kfpt):: &
 ape &                       ! exner function
,dth &                       ! potential temperature step
,dp &                        ! pressure step
,dthe &                      ! equivalent pot. temperature step
,p &                         ! pressure
,qs &                        ! saturation specific humidity
,the0k &                     ! base value for equivalent pot. temperature
,sthek &                     ! equivalent pot. temperature range
,th                          ! potential temperature

real(kind=kfpt),dimension(1:jtb):: &
 apt &                       ! temporary
,aqt &                       ! temporary
,tnew &                      ! new temperature
,told &                      ! old temperature
,thenew &                    ! new equivalent potential temperature
,theold &                    ! old equivalent potential temperature
,y2t                         ! temporary
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      dth=(thh-thl)/real(jtb-1)
      dp=(ph-pl)/real(itb-1)
      rdp=1./dp
!-----------------------------------------------------------------------
      p=pl-dp
      do kp=1,itb
        p=p+dp
        th=thl-dth
        do kth=1,jtb
          th=th+dth
          ape=(100000./p)**cappa
          qs=pq0/p*exp(a2*(th-a3*ape)/(th-a4*ape))
          told(kth)=th/ape
          theold(kth)=th*exp(eliwv*qs/(cp*told(kth)))
        enddo
!
        the0k=theold(1)
        sthek=theold(jtb)-theold(1)
        theold(1)=0.
        theold(jtb)=1.
!
        do kth=2,jtb-1
          theold(kth)=(theold(kth)-the0k)/sthek
!wwwwwwwwwwwwww fix due to 32 bit precision limitation wwwwwwwwwwwwwwwww
          if((theold(kth)-theold(kth-1)).lt.eps) then
            theold(kth)=theold(kth-1)+eps
          endif
!mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
        enddo
!
        the0(kp)=the0k
        sthe(kp)=sthek
!
        thenew(1)=0.
        thenew(jtb)=1.
        dthe=1./real(jtb-1)
        rdthe=1./dthe
!
        do kth=2,jtb-1
          thenew(kth)=thenew(kth-1)+dthe
        enddo
!
        y2t(1)=0.
        y2t(jtb)=0.
!
        call spline(jtb,jtb,theold,told,y2t,jtb,thenew,tnew,apt,aqt)
!
        do kth=1,jtb
          ttbl(kth,kp)=tnew(kth)
        enddo
!-----------------------------------------------------------------------
      enddo
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
      integer :: k,k1,k2,kold,noldm1
!
      real :: ak,bk,ck,den,dx,dxc,dxl,dxr,dydxl,dydxr                   &
             ,rdx,rtdxc,x,xk,xsq,y2k,y2kp1
!
!-----------------------------------------------------------------------
      noldm1=nold-1                                                     
!                                                                       
      dxl=xold(2)-xold(1)                                               
      dxr=xold(3)-xold(2)                                               
      dydxl=(yold(2)-yold(1))/dxl                                       
      dydxr=(yold(3)-yold(2))/dxr                                       
      rtdxc=.5/(dxl+dxr)                                                
!                                                                       
      p(1)= rtdxc*(6.*(dydxr-dydxl)-dxl*y2(1))                          
      q(1)=-rtdxc*dxr                                                   
!                                                                       
      if(nold.eq.3) go to 700                                           
!-----------------------------------------------------------------------
      k=3                                                               
!                                                                       
 100  dxl=dxr                                                           
      dydxl=dydxr                                                       
      dxr=xold(k+1)-xold(k)                                             
      dydxr=(yold(k+1)-yold(k))/dxr                                     
      dxc=dxl+dxr                                                       
      den=1./(dxl*q(k-2)+dxc+dxc)                                       
!                                                                       
      p(k-1)= den*(6.*(dydxr-dydxl)-dxl*p(k-2))                         
      q(k-1)=-den*dxr                                                   
!                                                                       
      k=k+1                                                             
      if(k.lt.nold) go to 100                                           
!-----------------------------------------------------------------------
 700  k=noldm1                                                          
!                                                                       
 200  y2(k)=p(k-1)+q(k-1)*y2(k+1)                                       
!                                                                       
      k=k-1                                                             
      if(k.gt.1) go to 200                                              
!-----------------------------------------------------------------------
      k1=1                                                              
!                                                                       
 300  xk=xnew(k1)                                                       
!                                                                       
      do 400 k2=2,nold                                                  
      if(xold(k2).le.xk) go to 400                                      
      kold=k2-1                                                         
      go to 450                                                         
 400  continue                                                          
      ynew(k1)=yold(nold)                                               
      go to 600                                                         
!                                                                       
 450  if(k1.eq.1)   go to 500                                           
      if(k.eq.kold) go to 550                                           
!                                                                       
 500  k=kold                                                            
!                                                                       
      y2k=y2(k)                                                         
      y2kp1=y2(k+1)                                                     
      dx=xold(k+1)-xold(k)                                              
      rdx=1./dx                                                         
!                                                                       
      ak=.1666667*rdx*(y2kp1-y2k)                                       
      bk=.5*y2k                                                         
      ck=rdx*(yold(k+1)-yold(k))-.1666667*dx*(y2kp1+y2k+y2k)            
!                                                                       
 550  x=xk-xold(k)                                                      
      xsq=x*x                                                           
!                                                                       
      ynew(k1)=ak*xsq*x+bk*xsq+ck*x+yold(k)                             
!                                                                       
 600  k1=k1+1                                                           
      if(k1.le.nnew) go to 300                                          
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
real(kind=kfpt),parameter:: &
 eps=0.000001
!--local variables------------------------------------------------------
integer(kind=kint):: &
 k                           ! index
 
real(kind=kfpt):: &
 x &                         ! temporary
,zeta1 &                     ! z/L, sea
,zeta2 &                     ! z/L, land
,zrng1 &                     ! z/L range, sea
,zrng2                       ! z/L range, land
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      kztm2=kztm-2
!
      fh01=1.
      fh02=1.
!
      ztmin1=-5.0
      ztmax1= 1.0
!
      ztmin2=-5.0
      ztmax2= 1.0
!
      zrng1=ztmax1-ztmin1
      zrng2=ztmax2-ztmin2
!
      dzeta1=zrng1/(kztm-1)
      dzeta2=zrng2/(kztm-1)
!--------------function definition loop---------------------------------
      zeta1=ztmin1
      zeta2=ztmin2
      do k=1,kztm
!--------------unstable range-------------------------------------------
        if(zeta1.lt.0.)then
!--------------paulson 1970 functions-----------------------------------
          x=sqrt(sqrt(1.-16.*zeta1))
          psim1(k)=-2.*log((x+1.)/2.)-log((x*x+1.)/2.)+2.*atan(x)-pihf
          psih1(k)=-2.*log((x*x+1.)/2.)
!--------------stable range---------------------------------------------
        else
!--------------holtslag and de bruin 1988-------------------------------
          psim1(k)=0.7*zeta1+0.75*zeta1*(6.-0.35*zeta1)*exp(-0.35*zeta1)
          psih1(k)=0.7*zeta1+0.75*zeta1*(6.-0.35*zeta1)*exp(-0.35*zeta1)
!-----------------------------------------------------------------------
        endif
!--------------unstable range-------------------------------------------
        if(zeta2.lt.0.)then
!--------------paulson 1970 functions-----------------------------------
          x=sqrt(sqrt(1.-16.*zeta2))
          psim2(k)=-2.*log((x+1.)/2.)-log((x*x+1.)/2.)+2.*atan(x)-pihf
          psih2(k)=-2.*log((x*x+1.)/2.)
!--------------stable range---------------------------------------------
        else
!--------------holtslag and de bruin 1988-------------------------------
          psim2(k)=0.7*zeta2+0.75*zeta2*(6.-0.35*zeta2)*exp(-0.35*zeta2)
          psih2(k)=0.7*zeta2+0.75*zeta2*(6.-0.35*zeta2)*exp(-0.35*zeta2)
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
        if(k.eq.kztm)then
          ztmax1=zeta1
          ztmax2=zeta2
        endif
!
        zeta1=zeta1+dzeta1
        zeta2=zeta2+dzeta2
!-----------------------------------------------------------------------
      enddo
!-----------------------------------------------------------------------
      ztmax1=ztmax1-eps
      ztmax2=ztmax2-eps
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
!-----------------------------------------------------------------------
!
      REAL*8 TIMEF
      INTEGER :: IC,IR
!
!-----------------------------------------------------------------------
!
      CALL SYSTEM_CLOCK(COUNT=IC,COUNT_RATE=IR)
      TIMEF=REAL(IC)/REAL(IR)*1000.
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
        CALL MPI_Finalize(irc)
        stop 911
!
!!!     call ESMF_Finalize(terminationflag=ESMF_ABORT)
!
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
