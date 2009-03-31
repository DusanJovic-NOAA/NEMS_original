!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        module module_surface
!-----------------------------------------------------------------------
use module_control,only:klog,kint,kfpt
!-----------------------------------------------------------------------
! land, ice, sea surface (liss) model
!
! z. janjic, august 2002, may 2008, bethesda
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------
integer(kind=kint),parameter:: &
 nwet=4 &           ! number of moisture layers
,ksnomx=1 &         ! maximum number of snow temperature layers
,nosno=nwet+1 &     ! number of ground/ice temperature layers
,kmx=nosno+ksnomx   ! maximum number of temperature layers

integer(kind=kint),parameter:: &
 nstm=30 &          ! max. allowed number of soil types
,nvtm=30 &          ! max. allowed number of vegetation types
,nsoil=21 &         ! actual number of soil types
,nveg=27            ! actual number of vegetation types

real(kind=kfpt),parameter:: &
 t0=273.15 &        ! freezing temperaure
,tice=273.15-1.7 &  ! salt water freezing temperature
,tf1=t0-3. &        ! freezing of soil water finished
,tf2=t0+1. &        ! freezing of soil water begins
,Tfm=(tf1+tf2)*0.5 &!
,rsgr=50 &          ! soil evaporation resistence
,rowliw=.334e9 &    ! volumetric water fusion heat capacity
,roci=2.05e6 &      ! ice volumetric heat capacity
,rocw=4.19e6 &      ! water volumetric heat capacity
,rroci=1./roci &    !
,vaki=2.260 &       ! ice volumetric heat conductivity
,aki=vaki*rroci &   ! ice heat conductivity
!vakm=3.44 &        ! soil matrix heat conductivity - replaced with peters-lidard
,vakq=7.7 &         ! quartz heat conductivity
,vako=2.0 &         ! heat conductivity of material other than quartz
,vakw=0.570 &       ! water volumetric heat conductivity
,roi=920.0 &        ! ice density
,rosno=300. &       ! snow density
,snofc=roi/rosno &  !
,rsnofc=1./snofc &  !
,rgmin=50. &        ! minimum ground resistance
,wlaymx=0.0002 &    ! maximum water amount on single leaf
,epssno=0.00001 &   ! minimum depth of snow layer (m)
,epsw=1.e-8 &       ! minimum soil moisture content
,pq0=379.90516 &    !
,a2=17.2693882 &    !
,a3=273.15 &        !
,a4=35.86 &         !
,a23m4=a2*(a3-a4) & !
,pi=3.14159 &       !
,hpi=pi*0.5 &       !
,rtfpi=pi/(tf2-tf1)&!
,r=287.04 &         ! dry air gas constant
,cp=1004.6 &        ! air heat capacity at constant pressure
,capa=r/cp &        ! r/cp
,row=1000. &        ! water density
,tresh=1.0 &        !
,pq0lnd=pq0*tresh & !
,pq0sea=pq0*0.98 &  ! correction for salt water
,elwv=2.50e6 &      !
,eliv=2.834e6 &     !
,reliw=eliv/elwv &  !
,eliw=.334e6 &      !
,stbol=5.67e-8 &    ! Stefan-Bolzman constant
,eintcp=0.25 &      ! interception efficiency
,tsno=273.15        ! snow temperature

!--soil properties------------------------------------------------------
real(kind=kfpt),save:: &
 roctbl(nstm) &     ! volumetric heat capacity
,vaktbl(nstm) &     ! volumetric heat diffusivity
,qtztbl(nstm) &     ! volumetric heat diffusivity
,akgtbl(nstm) &     ! quartz content
,betbl(nstm) &      ! Clapp and Hornberger nondimensional exponent
,psstbl(nstm) &     ! matric potential at saturation
,gamtbl(nstm) &     ! hydraulic conductivity at saturation
,akwtbl(nstm) &     ! hydraulic diffusivity at saturation
,wgstbl(nstm) &     ! soil moisture at saturation
,wgctbl(nstm) &     ! soil moisture at field capacity
,wgwtbl(nstm)       ! soil moisture at permanent wilting point

!--vegetation properties------------------------------------------------
real(kind=kfpt),save:: &
 cvetbl(nvtm) &     ! maximum vegetation cover
,elatbl(nvtm) &     ! leaf area index
,rsmtbl(nvtm) &     ! minimum stomatal resistence
,hvetbl(nvtm) &     ! air moisture defficit factor
,arftbl(nvtm) &     ! root fraction parameter a
,brftbl(nvtm)       ! root fraction parameter b

real(kind=kfpt),save:: &
 rfrtbl(nwet,nvtm)  ! root fraction

!--discretization-------------------------------------------------------
real(kind=kfpt),save:: &
 dgztbl(nosno) &    ! basic soil layer depths
,dgttbl(kmx,nstm) & ! actual soil/snow layer depths for temperature
,dgwtbl(nwet,nstm)  ! soil layer depths for wetness
!-----------------------------------------------------------------------
! soil  statsgo
! type  class
! ----  -------
!   1   sand
!   2   loamy sand
!   3   sandy loam
!   4   silt loam
!   5   silt
!   6   loam
!   7   sandy clay loam
!   8   silty clay loam
!   9   clay loam
!  10   sandy clay
!  11   silty clay
!  12   clay
!  13   organic material
!  14   water
!  15   bedrock
!  16   other(land-ice)
!  17   playa
!  18   lava
!  19   white sand
!--------------------
!  20   generic ecmwf
!  21   sea ice
!-----------------------------------------------------------------------
!data roctbl &
!/    1.27e6,    1.33e6,    1.35e6,    1.48e6,    1.48e6,   1.42e6 &
!,    1.44e6,    1.48e6,    1.50e6,    1.51e6,    1.63e6,   1.63e6 &
!,    0.19e6,    4.18e6,    2.10e6,    2.05e6,    2.10e6,   2.10e6 &
!,    1.27e6,    2.19e6,    2.05e6,    0.00  ,    0.00  ,   0.00   &
!,    0.00  ,    0.00  ,    0.00  ,    0.00  ,    0.00  ,   0.00  /

data roctbl & ! noah costant volumetric heat capacity
/    2.19e6,    2.19e6,    2.19e6,    2.19e6,    2.19e6,   2.19e6 &
,    2.19e6,    2.19e6,    2.19e6,    2.19e6,    2.19e6,   2.19e6 &
,    2.19e6,    4.18e6,    2.19e6,    2.05e6,    2.19e6,   2.19e6 &
,    2.19e6,    2.19e6,    2.05e6,    2.19e6,    2.19e6,   2.19e6 &
,    2.19e6,    2.19e6,    2.19e6,    2.19e6,    2.19e6,   2.19e6/

data vaktbl & ! independent on wgs, oterwise peters-lidard
/    0.000,     0.000,     0.000,     0.000,     0.000,    0.000  &
,    0.000,     0.000,     0.000,     0.000,     0.000,    0.000  &
,    0.000,     0.570,     0.000,     2.260,     0.000,    0.000  &
,    0.000,     0.000,     2.260,     0.000,     0.000,    0.000  &
,    0.000,     0.000,     0.000,     0.000,     0.000,    0.000 /

data qtztbl & ! independent on wgs
/    0.920,     0.820,     0.600,     0.250,     0.100,    0.400  &
,    0.600,     0.100,     0.350,     0.520,     0.100,    0.250  &
,    0.050,     0.000,     0.070,     0.000,     0.600,    0.520  &
,    0.920,     0.600,     0.000,     0.000,     0.000,    0.000  &
,    0.000,     0.000,     0.000,     0.000,     0.000,    0.000 /

data betbl & !(noah bb)
/    4.05,      4.26,      4.74,      5.33,      5.33,     5.25   &
,    6.77,      8.72,      8.17,     10.73,     10.39,    11.55   &
,    5.25,      0.00,      4.05,      4.26,     11.55,     4.05   &
,    4.05,      6.04,      0.00,      0.00,      0.00,     0.00   &
,    0.00,      0.00,      0.00,      0.00,      0.00,     0.00  /

data psstbl & !(noah -satpsi)
/   -0.0350,   -0.0363,   -0.1413,   -0.7586,   -0.7586,  -0.3548 &
,   -0.1349,   -0.6166,   -0.2630,   -0.0977,   -0.3236,  -0.4677 &
,   -0.3548,    0.0000,   -0.0350,   -0.0363,   -0.4677,  -0.0350 &
,   -0.0350,   -0.3380,   -0.0000,   -0.0000,   -0.0000,  -0.0000 &
,   -0.0000,   -0.0000,   -0.0000,   -0.0000,   -0.0000,  -0.0000/

data gamtbl & !(noah satdk)
/1.7600e-4, 1.4078e-5, 5.2304e-6, 2.8089e-6, 2.8089e-6, 3.3770e-6 &
,4.4518e-6, 2.0348e-6, 2.4464e-6, 7.2199e-6, 1.3444e-6, 9.7394e-7 &
,3.3770e-6, 0.0000   , 1.4078e-5, 1.4078e-5, 9.7394e-7, 1.4078e-5 &
,1.7600e-4, 4.5700e-4, 0.0000   , 0.0000   , 0.0000   , 0.0000    &
,0.0000   , 0.0000   , 0.0000   , 0.0000   , 0.0000   , 0.0000   /

data wgstbl & !(noah maxsmc)
/    0.395,     0.421,     0.434,     0.476,     0.476,    0.439  &
,    0.404,     0.464,     0.465,     0.406,     0.468,    0.457  &
,    0.464,     0.000,     0.200,     0.000,     0.457,    0.200  &
,    0.395,     0.472,     0.000,     0.000,     0.000,    0.000  &
,    0.000,     0.000,     0.000,     0.000,     0.000,    0.000 /
!-----------------------------------------------------------------------
! class usgs-wrf vegetation/surface type
!   1   Urban and Built-Up Land 11
!   2   Dryland Cropland and Pasture 1
!   3   Irrigated Cropland and Pasture 10
!   4   Mixed Dryland/Irrigated Cropland and Pasture (1+10)/2
!   5   Cropland/Grassland Mosaic [1+(2+7)/2]/2
!   6   Cropland/Woodland Mosaic (1+18)/2
!   7   Grassland (2+7)/2
!   8   Shrubland (16+17)/2
!   9   Mixed Shrubland/Grassland (2+7+16+17)/4
!  10   Savanna (2+7+18)/3
!  11   Deciduous Broadleaf Forest   5
!  12   Deciduous Needleleaf Forest  4
!  13   Evergreen Broadleaf Forest   6
!  14   Evergreen Needleleaf Forest  3
!  15   Mixed Forest 18
!  16   Water Bodies 14,15
!  17   Herbaceous Wetland 13
!  18   Wooded Wetland (18+13)/2
!  19   Barren or Sparsely Vegetated 11
!  20   Herbaceous Tundra (2+9)/2
!  21   Wooded Tundra (18+9)/2
!  22   Mixed Tundra 9
!  23   Bare Ground Tundra 9nocveg
!  24   Snow or Ice
!  25   Playa
!  26   Lava
!  27   White Sand
!-----------------------------------------------------------------------
data rsmtbl & !(noah rsmtbl)
/  200.,      180.,      180.,      180.,      140.,      200.     &
,  100.,      225.,      170.,      120.,      175.,      500.     &
,  240.,      500.,      250.,     9999.,      040.,      100.     &
,  300.,      080.,      080.,      080.,      080.,     9999.     &
,  040.,      100.,      300.,     9999.,     9999.,     9999.    /

!data rsmtbl & !(noah rsmtbl, noah values)
!/  200.,      070.,      070.,      070.,      070.,      070.     &
!,  070.,      300.,      170.,      070.,      100.,      150.     &
!,  150.,      250.,      150.,     9999.,      040.,      100.     &
!,  300.,      150.,      150.,      150.,      200.,     9999.     &
!,  040.,      100.,      300.,     9999.,     9999.,     9999.    /

data elatbl & !(noah lai_data)
/    2.5,       3.0,       3.0,       3.0,       3.0,       4.0    &
,    2.0,       2.3,       2.1,       3.0,       5.0,       5.0    &
,    6.0,       5.0,       5.0,       0.0,       2.0,       4.5    &
,    4.0,       1.0,       2.0,       1.0,       1.0,       0.0    &
,    3.0,       3.0,       3.0,       0.0,       0.0,       0.0   /

!data elatbl & !(noah lai_data, noah values)
!/    4.0,       4.0,       4.0,       4.0,       4.0,       4.0    &
!,    4.0,       4.0,       4.0,       4.0,       4.0,       4.0    &
!,    4.0,       4.0,       4.0,       4.0,       4.0,       4.0    &
!,    4.0,       4.0,       4.0,       4.0,       4.0,       0.0    &
!,    4.0,       4.0,       4.0,       0.0,       0.0,       0.0   /

data cvetbl & !(noah shdfac)
/    0.10,      0.90,      0.90,      0.90,      0.80,      0.80   &
,    0.80,      0.50,      0.70,      0.50,      0.90,      0.90   &
,    0.99,      0.90,      0.90,      0.00,      0.60,      0.60   &
,    0.01,      0.60,      0.60,      0.60,      0.01,      0.50   &
,    0.00,      0.00,      0.00,      0.00,      0.00,      0.00  /

data hvetbl &
/     .0000,     .0000,     .0000,     .0000,     .0000,     .0001 &
,     .0000,     .0000,     .0000,     .0000,     .0003,     .0003 &
,     .0003,     .0003,     .0003,     .0000,     .0000,     .0003 &
,     .0000,     .0000,     .0003,     .0000,     .0000,     .0000 &
,     .0000,     .0000,     .0000,     .0000,     .0000,     .0000/

data arftbl &
/    10.74,     5.558,     5.558,     5.558,     6.326,     6.326  &
,    6.326,     6.326,     6.326,     6.326,     5.990,     7.066  &
,    7.344,     6.706,     4.453,     0.000,     8.992,     8.992  &
,   10.740,     8.992,     8.992,     8.992,     8.992,    10.740  &
,   10.740,    10.740,    10.740,     0.000,     0.000,     0.000 /

data brftbl &
/    2.608,     2.614,     2.614,     2.614,     1.567,     1.567  &
,    1.567,     1.567,     1.567,     1.567,     1.955,     1.953  &
,    1.303,     2.175,     1.631,     0.000,     8.992,     8.992  &
,    2.608,     8.992,     8.992,     8.992,     8.992,     2.608  &
,    2.608,     2.608,     2.608,     0.000,     0.000,     0.000 /
!-----------------------------------------------------------------------
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        subroutine setsfc
!-----------------------------------------------------------------------
! constants for multiple layer svat land surface model
! z. janjic, august 2002, may 2008, bethesda
!-----------------------------------------------------------------------
implicit none
!--local variables------------------------------------------------------
integer(kind=kint):: &
 k &                ! index
,n                  ! index

real(kind=kfpt):: &
 wgc &              !
,wgg &              !
,wgw &              !
,zgu &              !
,zgd &              !
,srfr               !
!---set resolution parameters-------------------------------------------
!dgztbl(1)=0.02 ! ecmwf
!dgztbl(2)=0.07-dgztbl(1) ! ecmwf
!dgztbl(3)=0.21 ! ecmwf
!dgztbl(4)=0.72 ! ecmwf
!dgztbl(5)=1.89 ! ecmwf

dgztbl(1)=0.00
dgztbl(2)=0.10-dgztbl(1) ! noah
dgztbl(3)=0.30
dgztbl(4)=0.60
dgztbl(5)=1.00

!dgztbl(1)=.02
!do k=2,nosno
!  dgztbl(k)=dgztbl(k-1)*1.25
!enddo
!---compute dependent soil parameters-----------------------------------
do n=1,nsoil
  if(wgstbl(n).gt.epsw) then ! those dependent on wgs
    wgg=(1.-wgstbl(n))*2700.
    vaktbl(n)=(wgg*0.135+64.7)/(2700.-wgg*0.947) ! peters-lidard
  endif
enddo
!---compute dependent soil moisture parameters--------------------------
do n=1,nsoil
  if(wgstbl(n).gt.epsw) then
    akwtbl(n)=betbl(n)*gamtbl(n)*(-psstbl(n)/wgstbl(n)) ! noah satdw
  else
    akwtbl(n)=0.
  endif
  if(gamtbl(n).gt.epsw) then
    wgc=wgstbl(n)*(5.79e-9/gamtbl(n))**(1.0/(2.0*betbl(n)+3.0)) ! noah refsmc1
  else
    wgc=0.
  endif
  wgctbl(n)=wgc+(wgstbl(n)-wgc)/3.0 ! noah refsmc
  if(betbl(n).gt.epsw) then
    wgw=wgstbl(n)*(-200.0/psstbl(n))**(-1.0/betbl(n)) ! noah wltsmc1
  else
    wgw=0.
  endif
  wgwtbl(n)=wgw-0.5*wgw ! noah wltsmc
enddo
!--layer thicknesses----------------------------------------------------
do n=1,nsoil
  akgtbl(n)=vaktbl(n)/roctbl(n)
  do k=1,kmx
    dgttbl(k,n)=999.
  enddo
  do k=1,nosno
    dgttbl(k,n)=dgztbl(k)
  enddo
enddo

do n=1,nveg
  dgwtbl(1,n)=dgztbl(1)+dgztbl(2)
  do k=2,nwet
    dgwtbl(k,n)=dgztbl(k+1)
  enddo
enddo
!--compute root fraction------------------------------------------------
do n=1,nveg
  zgu=0.
  do k=1,nwet
    zgd=zgu+dgwtbl(k,n)
    rfrtbl(k,n)=0.5*(exp(-arftbl(n)*zgu)+exp(-brftbl(n)*zgu) &
                    -exp(-arftbl(n)*zgd)-exp(-brftbl(n)*zgd))
    zgu=zgd
  enddo
  srfr=0.
  do k=1,nwet
    srfr=rfrtbl(k,n)+srfr
  enddo
  if(srfr.gt.0.) then
    srfr=1./srfr
  else
    srfr=0.
  endif
  do k=1,nwet
    rfrtbl(k,n)=rfrtbl(k,n)*srfr
  enddo
enddo
!-----------------------------------------------------------------------
                        endsubroutine setsfc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        subroutine runsfc &
(insoil,inveg,sm,sst,sice,vegfrc &
,kmp,dtphys &
,epsr,radin,rswin,akhs &
,plm,rxnrlm,tlm,qlm,clm &
,ps,rxnrs,ths,qs,tlbc &
,prra,prsn,sno,smelt,roff &
,tflux,gflux,qflux,pflux &
,tg,wliq,wg)
!-----------------------------------------------------------------------
implicit none
!--global variables-----------------------------------------------------
integer(kind=kint),intent(in):: &
 insoil &           ! soil type index
,inveg              ! vegetation type index

real(kind=kfpt),intent(in):: &
 sm &               ! sea mask
,sst &              ! SST
,sice &             ! sea-ice mask, 0. no, 1. yes
,vegfrc &           ! greenness
,dtphys &           ! physics time step
,epsr &             ! emissivity
,radin &            ! total radiation at surface
,rswin &            ! short wave radiation W/m2
,akhs &             ! surface layer exchange coeff. / sfc. layer depth
,plm &              ! pressure at the lowest atmospheric level
,rxnrlm &           ! 1./exner function at the lowest atmospheric level
,tlm &              ! temperature at the lowest atmospheric level
,qlm &              ! specific humidity at the lowest atmospheric level
,clm &              ! condensate at the lowest atmospheric level
,ps &               ! surface pressure
,rxnrs &            ! 1./exner function at the surface
,tlbc               ! fixed lower BC for ground temperature

integer(kind=kint),intent(inout):: &
 kmp                ! # of layers, variable if snow pack on top

real(kind=kfpt),intent(inout):: &
 wliq &             ! water in interception reservoir (m)
,prra &             ! time step rain precipitation (m)
,prsn &             ! time step snow precipitation (m)
,sno &              ! sno accumulation (m)
,smelt &            ! amount of snow melted
,roff &             ! surface runoff (m)
,ths &              ! surface potential temperature
,qs &               ! surface saturation specific humidity
,tflux &            ! surface sensible heat flux, met. grad. sign
,gflux &            ! ground heat flux          , met. grad. sign
,qflux &            ! latent heat flux          , met. grad. sign
,pflux              ! potential evaporation flux, met. grad. sign

real(kind=kfpt),intent(inout):: &
 tg(kmx) &          ! soil/ice/snow temperature, midle of layers
,wg(nwet)           ! soil moisture, midle of layers

!--local variables------------------------------------------------------
integer(kind=kint):: &
 k &                ! working index, increases downward
,kmnew &            !
,ksno               ! # of snow layers

real(kind=kfpt):: &
 wgs &              ! saturation soil moisture
,wgcap &            ! soil moisture at field capacity
,wgpwp &            ! soil moisture at permanent wilting point
,be      &          ! Clapp and Hornberger nondimensional exponent
,gammas &           ! hydraulic conductivity at saturation
,akws &             ! hydraulic diffusivity at saturation
,cveg &             ! vegetation fraction
,elai &             ! leaf area index
,rsmin &            ! minimum stomatal resistance
,hveg &             ! high vegetation parameter
,dsno &             ! depth of snow layers
,qas &              ! saturation specific humidity at the lowest model layer
,avglm &            ! weight in balance eq., 0.<= avllm <=1.
,avgsf &            ! weight in balance eq., 0.<= avlsf <=1.
,agrlm &            ! weight in balance eq., 0.<= avllm <=1.
,agrsf &            ! weight in balance eq., 0.<= avlsf <=1.
,avllm &            ! weight in balance eq., 0.<= avllm <=1.
,avlsf &            ! weight in balance eq., 0.<= avlsf <=1.
,trf &              ! reference temperature for linearization
,pevap &            ! potential evaporation (m)
,evap &             ! total evaporation (m)
,eveg &             ! evapotranspiration (m)
,egrnd &            ! bare ground evaporation (m)
,qtz &              ! quartz content
,snoht              ! snow height

real(kind=kfpt):: &
 dgw(nwet) &        ! depth of layers for soil moisture
,rfr(nwet) &        ! root fraction
,wgi(nwet) &        ! frozen soil water
,wgl(nwet) &        ! liquid soil water
,dgt(kmx) &         ! depth of soil/ice/snow layers, if dgt(1)=0., tg(1) is skin t.
,roc(kmx) &         ! volumetric heat capacity
,akg(kmx) &         ! volumetric heat diffusivity
,akt(kmx)           ! actual heat diffusivity

!-----------------------------------------------------------------------
if(sm.lt.0.5) then ! solid surface, land or ice
!-----------------------------------------------------------------------
! soil properties
  qtz=qtztbl(insoil)
  wgs=wgstbl(insoil)
  wgcap=wgctbl(insoil)
  wgpwp=wgwtbl(insoil)

  be=betbl(insoil)
  gammas=gamtbl(insoil)
  akws=akwtbl(insoil)

  do k=1,nwet
    dgw(k)=dgwtbl(k,insoil)
  enddo

! vegetation properties
  cveg=cvetbl(inveg)*vegfrc
  elai=elatbl(inveg)
  rsmin=rsmtbl(inveg)
  hveg=hvetbl(inveg)

  do k=1,nwet
    rfr(k)=rfrtbl(k,inveg)
  enddo
!-----------------------------------------------------------------------
  if(sno.le.epssno) then ! land/sea-ice without snow cover
!-----------------------------------------------------------------------
    do k=1,nosno
      dgt(k)=dgttbl(k,insoil)
      roc(k)=roctbl(insoil)
      akg(k)=akgtbl(insoil)
    enddo

    if(kmp.gt.nosno) then ! snow just dissapeared
      ksno=kmp-nosno
      do k=2,nosno ! tg(1) remains the same as in the melted snow slab
        tg(k)=tg(k+ksno)
      enddo
    endif ! end of snow just dissapeared

    kmp=nosno
    ksno=0
!-----------------------------------------------------------------------
  else ! land/sea-ice with snow cover
!-----------------------------------------------------------------------
    snoht=sno*snofc
    ksno=min(int(snoht/dgztbl(2))+1,ksnomx)
    kmnew=nosno+ksno

    dgt(1)=min(sno*snofc*0.1,dgztbl(1)) ! snow surface slab
    dsno=(snoht-dgt(1))/float(ksno)

    do k=1,ksno
      dgt(k+1)=dsno
    enddo

    do k=1,ksno+1
      akg(k)=aki
      roc(k)=roci*rsnofc
    enddo

    do k=2,nosno
      dgt(k+ksno)=dgttbl(k,insoil)
      roc(k+ksno)=roctbl(insoil)
      akg(k+ksno)=akgtbl(insoil)
    enddo

    if(kmp.eq.nosno) then ! new snow
      do k=nosno,2,-1
        tg(k+ksno)=tg(k)
      enddo
      do k=2,ksno+1
        tg(k)=tg(1) ! tg(1) remains the same as in top soil slab
      enddo
    else ! old snow
      if    (kmnew.gt.kmp) then ! more layers
        do k=kmp,2,-1
          tg(k+kmnew-kmp)=tg(k) ! move old layer values
        enddo
        do k=3,kmnew-kmp+1
          tg(k)=tg(2) ! add values
        enddo
      elseif(kmnew.lt.kmp) then ! less layers
        do k=3,kmnew
          tg(k)=tg(k+kmp-kmnew)
        enddo
      endif
    endif

    kmp=kmnew
!-----------------------------------------------------------------------
  endif ! end of snow/nosnow branch
!-----------------------------------------------------------------------
  call preheat &
  (kmp,ksno &
  ,sno,rswin,akhs,wliq &
  ,cveg,elai,rsmin,hveg &
  ,wgs,wgcap,wgpwp,qtz &
  ,dgw,rfr,wg,wgi,wgl &
  ,dgt,roc,akg,tg &
  ,plm,tlm,qlm,qas &
  ,avglm,avgsf,agrlm,agrsf,avllm,avlsf,akt)
!
  trf=tg(1)
!
  call heat &
  (kmp,dtphys &
  ,sice,wgs,sno &
  ,plm,rxnrlm,tlm,qlm,qas,clm,tlbc,akhs &
  ,epsr,radin &
  ,trf,ps,rxnrs &
  ,avglm,avgsf,agrlm,agrsf,avllm,avlsf &
  ,ths,qs,dgt,roc,akt,tg &
  ,pevap,evap,eveg,egrnd,smelt &
  ,tflux,gflux,qflux,pflux)
!
  call water &
  (dtphys,cveg,elai &
  ,wgs,wgpwp,be,gammas,akws &
  ,pevap,evap,eveg,smelt &
  ,dgw,rfr &
  ,prra,prsn,sno,roff &
  ,wliq,wg,wgi,wgl)
!-----------------------------------------------------------------------
else ! liquid surface
!-----------------------------------------------------------------------
  ths=sst*rxnrs
!
  qs=pq0sea/ps*exp(a2*(sst-a3)/(sst-a4))
  sno=0.
  prsn=0.
! marine surface fluxes computed elsewhere because of thz0 and qz0
  tflux=9999.
  gflux=9999.
  qflux=9999.
  pflux=9999.
!-----------------------------------------------------------------------
endif ! end of solid/liquid surface branch
!-----------------------------------------------------------------------
                        endsubroutine runsfc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        subroutine preheat &
(kmp,ksno &
,sno,rswin,akhs,wliq &
,cveg,elai,rsmin,hveg &
,wgs,wgcap,wgpwp,qtz &
,dgw,rfr,wg,wgi,wgl &
,dgt,roc,akg,tg &
,plm,tlm,qlm,qas &
,avglm,avgsf,agrlm,agrsf,avllm,avlsf,akt)
!-----------------------------------------------------------------------
!  Multiple layer thermal and moisture surface processes
!  Janjic, 2002
!-----------------------------------------------------------------------
implicit none
!--global variables-----------------------------------------------------
integer(kind=kint),intent(in):: &
 kmp &              ! # of layers, variable if snow pack on top
,ksno               ! # of sno layers

real(kind=kfpt),intent(in):: &
 sno &              ! snow amount in equivalent water depth (m)
,rswin &            ! short wave radiation W/m2
,akhs &             ! surface layer exchange coeff. divided by sfc. layer depth
,wliq &             ! water in interception reservoir
,cveg &             ! vegetation fraction
,elai &             ! leaf area index
,rsmin &            ! minimum stomatal resistance
,hveg &             ! high vegetation parameter
,wgs &              ! saturation soil moisture
,wgcap &            ! soil moisture at field capacity
,wgpwp &            ! soil moisture at permanent wilting point
,plm &              ! pressure at the lowest atmospheric level
,tlm &              ! temperature at the lowest atmospheric level
,qlm &              ! specific humidity at the lowest atmospheric level
,qtz                ! quartz content

real(kind=kfpt),intent(in):: &
 rfr(nwet) &        ! root fraction
,dgw(nwet) &        ! depth of layers for soil moisture
,dgt(kmx) &         ! depth of soil/ice/snow layers, if dgt(1)=0., tg(1) is skin t.
,akg(kmx) &         ! volumetric heat diffusivity
,tg(kmx)            ! soil/ice/snow temperature, midle of layers

real(kind=kfpt),intent(inout):: &
 wg(nwet) &         ! soil moisture, midle of layers
,roc(kmx)           ! volumetric heat capacity

real(kind=kfpt),intent(out):: &
 qas &              ! saturation specific humidity at the lowest model layer
,avglm &            ! weight in balance eq., 0.<= avllm <=1.
,avgsf &            ! weight in balance eq., 0.<= avlsf <=1.
,agrlm &            ! weight in balance eq., 0.<= avllm <=1.
,agrsf &            ! weight in balance eq., 0.<= avlsf <=1.
,avllm &            ! weight in balance eq., 0.<= avllm <=1.
,avlsf              ! weight in balance eq., 0.<= avlsf <=1.

real(kind=kfpt),intent(out):: &
 wgi(nwet) &        ! frozen part of soil moisture
,wgl(nwet) &        ! liquid part of soil moisture
,akt(kmx)           ! soil/ice/snow thermal diffusivity, middle of the layers
!--local variables------------------------------------------------------
integer(kind=kint):: &
 k &                ! working index, increases downwards
,koff               ! offset index

real(kind=kfpt):: &
 tgk &              ! tg(k)
,arg &              !
,aks &              ! dry soil heat conductivity
,fcs &              !
,wlmx &             !
,rwlmx &            !
,cliq &             !
,rroc &             ! 1./roc
,rwg &              !
,ce &               !
,akm &              !
,cavl &             !
,rf1 &              ! resistence function
,rf2 &              ! resistence function
,rf3 &              ! resistence function
,wgbar &            ! mean wg
,rc &               ! resistence
,alpha &            !
,vaks               !

real(kind=kfpt):: &
 tfr(kmx) &         ! frozen soil fraction
,fr1(kmx)           ! heat capacity increase due to fusion
!-----------------------------------------------------------------------
qas=pq0/plm*exp(a2*(tlm-a3)/(tlm-a4))
!-----------------------------------------------------------------------
if(wgs.gt.epsw) then ! soil points
!---limit wg by saturation soil moisture, if needed in the first step---
  do k=1,nwet
    wg(k)=min(wg(k),wgs)
  enddo
!---modify heat capacity for water content and freezing soil water------
  if(kmp.eq.nosno) then ! no sno
    do k=1,nosno
      tgk=tg(k)
      if    (tgk.gt.tf2) then
        tfr(k)=0.
        fr1(k)=0.
      elseif(tgk.lt.tf1) then
        tfr(k)=1.
        fr1(k)=0.
      else
        arg=(tgk-tfm)*rtfpi
        tfr(k)=(1.-sin(arg))*0.5
        fr1(k)=cos(arg)*rtfpi*0.5
      endif
    enddo
    wgi(1)=(tfr(1)*dgt(1)+tfr(2)*dgt(2))/dgw(1)*wg(1)
    wgl(1)=wg(1)-wgi(1)
    do k=2,nwet
      wgi(k)=tfr(k+1)*wg(k)
      wgl(k)=wg(k)-wgi(k)
    enddo
    roc(1)=(1.-wgs)*roc(1) &
          +wgl(1)*rocw+wgi(1)*roci+wg(1)*fr1(1)*rowliw
    do k=2,nosno
      roc(k)=(1.-wgs)*roc(k) &
            +wgl(k-1)*rocw+wgi(k-1)*roci+wg(k-1)*fr1(k)*rowliw
    enddo
  else ! sno
    do k=ksno+2,kmp
      tgk=tg(k)
      if    (tgk.gt.tf2) then
        tfr(k)=0.
        fr1(k)=0.
      elseif(tgk.lt.tf1) then
        tfr(k)=1.
        fr1(k)=0.
      else
        arg=(tgk-tfm)*rtfpi
        tfr(k)=(1.-sin(arg))*0.5
        fr1(k)=cos(arg)*rtfpi*0.5
      endif
    enddo
    do k=1,nwet
      wgi(k)=tfr(k+ksno+1)*wg(k)
      wgl(k)=wg(k)-wgi(k)
    enddo
    do k=ksno+2,kmp
      roc(k)=(1.-wgs)*roc(k) &
            +wgl(k-ksno-1)*rocw+wgi(k-ksno-1)*roci &
            +wg(k-ksno-1)*fr1(k)*rowliw
    enddo
  endif ! end of sno/no sno branch
!-----------------------------------------------------------------------
else ! ice, sea or glacier
!-----------------------------------------------------------------------
  do k=1,nwet
    wgi(k)=0.
    wgl(k)=0.
  enddo
!-----------------------------------------------------------------------
endif ! end of soil/ice branch
!--thermal diffusivity profile------------------------------------------
if(wgs.gt.epsw) then ! soil points
!-----------------------------------------------------------------------
  if(sno.le.epssno) then ! soil diffusivity, peters-lidard
    koff=1
  else
    koff=kmp-nwet
  endif
  do k=1,nwet
    if(wgs.gt.epsw) then
      rwg=wg(k)/wgs
    else
      rwg=0.
    endif
    if(wgi(k).le.epsw) then
      if(rwg.gt.0.1) then ! Kersten #
        ce=log10(rwg)+1.
      else
        ce=0.
      endif
    else
      ce=rwg
    endif
    rroc=1./roc(k+koff)
    aks=akg(k+koff)
!    akm=(vakm**(1.-wgs))*(vakw**wgl(k))*(vaki**wgi(k))*rroc ! retired
    vaks=(vakq**qtz)*(vako**(1.-qtz))
    akm=(vaks**(1.-wgs))*(vakw**wgl(k))*(vaki**wgi(k))*rroc
    akt(k+koff)=(akm-aks)*ce+aks
  enddo
  if(sno.le.epssno) then
    akt(koff)=akt(koff+1)
  else ! snow conductivity
    fcs=rsnofc !**1.88
    do k=1,ksno+1
      akt(k)=aki*fcs
    enddo
  endif
!-----------------------------------------------------------------------
else ! ice points
!-----------------------------------------------------------------------
  if(sno.le.epssno) then ! no sno
    do k=1,nosno
      akt(k)=akg(k)
    enddo
  else ! snow
    fcs=rsnofc !**1.88
    do k=1,ksno+1
      akt(k)=aki*fcs
    enddo
    do k=ksno+2,kmp
      akt(k)=akg(k)
    enddo
  endif
!-----------------------------------------------------------------------
endif ! end of soil/sea-ice branching
!--surface and lowest model level moisture availabilities---------------
if(wgs.gt.epsw.and.sno.le.epssno) then ! snow free land points
!-----------------------------------------------------------------------
  wlmx=(cveg*elai+(1.-cveg))*wlaymx ! interception area
  rwlmx=1./wlmx
  cliq=min(wliq*rwlmx,1.) ! interception reservoir availability

  avllm=cliq
    if(tg(1).le.t0) then
      avllm=avllm*reliw
    endif
  avlsf=avllm

  if(wgs.gt.epsw) then ! soil that takes water
    rf1=min((rswin*0.004+0.05)/(0.81*(rswin*0.004+1.)),1.)
    arg=(qlm-qas)*hveg*plm
    if(abs(arg).lt.2.) then
      rf3=exp(arg) ! high veg. moisture deficit
    else
      rf3=1.
    endif
    wgbar=0.
    do k=1,nwet
      wgbar=max(wgl(k),wgpwp)*rfr(k)+wgbar
    enddo
    if(wgbar.gt.wgpwp.and.elai.gt.0.) then
      rf2=min((wgbar-wgpwp)/(wgcap-wgpwp),1.)
      rc=rsmin/(elai*rf1*rf2*rf3)
      cavl=1./(1.+rc*akhs)
    else
      cavl=0.
    endif

    avglm=(1.-cliq)*cveg*cavl ! evapotranspiration
    avgsf=avglm

    avllm=avllm+avglm
    avlsf=avlsf+avgsf

    if(wgl(1).gt.wgpwp) then
      rf2=min((wgl(1)-wgpwp)/(wgcap-wgpwp),1.)
      rc=rsgr/rf2
      cavl=1./(1.+rc*akhs)
    else
      cavl=0.
    endif

    agrlm=(1.-cliq)*(1.-cveg)*cavl
    agrsf=agrlm

    avllm=avllm+agrlm
    avlsf=avlsf+agrsf
  else ! surface that doesn't take water
    avglm=0.
    avgsf=0.
    agrlm=0.
    agrsf=0.
  endif ! end of takes/doesn't take water branch
!-----------------------------------------------------------------------
else ! snow and ice points
!-----------------------------------------------------------------------
  avglm=0.
  avgsf=0.
  agrlm=0.
  agrsf=0.
  avllm=reliw
  avlsf=reliw
!-----------------------------------------------------------------------
endif
!-----------------------------------------------------------------------
                        endsubroutine preheat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        subroutine heat &
(kmp,dtphys &
,sice,wgs,sno &
,plm,rxnrlm,tlm,qlm,qas,clm,tlbc,akhs &
,epsr,radin &
,trf,ps,rxnrs &
,avglm,avgsf,agrlm,agrsf,avllm,avlsf &
,ths,qs,dgt,roc,akt,tg &
,pevap,evap,eveg,egrnd,smelt &
,tflux,gflux,qflux,pflux)
!-----------------------------------------------------------------------
!  Multiple layer thermal land/ice/snow surface processes
!  Janjic, 2001
!-----------------------------------------------------------------------
implicit none
!--global variables-----------------------------------------------------
integer(kind=kint),intent(in):: &
 kmp                ! # of layers, variable if snow pack on top

real(kind=kfpt),intent(in):: &
 dtphys &           ! physics time step
,sice &             ! sea-ice mask, 0. no, 1. yes
,wgs &              ! saturation soil moisture
,sno &              ! snow amount in equivalent water depth (m)
,plm &              ! pressure at the lowest atmospheric level
,rxnrlm &           ! pressure **R/cp at the lowest atmospheric level
,tlm &              ! temperature at the lowest atmospheric level
,qlm &              ! specific humidity at the lowest atmospheric level
,qas &              ! saturation spec. humidity in the air
,clm &              ! condensate at the lowest atmospheric level
,akhs &             ! surface layer exchange coeff. divided by sfc. layer depth
,epsr &             ! emissivity of the surface
,radin &            ! total incoming radiation at the surface
,trf &              ! temperature to linearize surface balance eq. around
,ps &               ! surface pressure
,rxnrs &            ! surface pressure **R/cp
,avglm &            ! weight in balance eq., 0.<= avllm <=1.
,avgsf &            ! weight in balance eq., 0.<= avlsf <=1.
,agrlm &            ! weight in balance eq., 0.<= avllm <=1.
,agrsf &            ! weight in balance eq., 0.<= avlsf <=1.
,tlbc               ! fixed lower BC for ground temperature

real(kind=kfpt),intent(in):: &
 dgt(kmx) &         ! depth of soil/ice/snow layers, if dgt(1)=0., tg(1) is skin t.
,roc(kmx)           ! volumetric heat capacity

real(kind=kfpt),intent(inout):: &
 ths &              ! skin[for dgt(1)=0.]/surface slab potential temperature
,qs &               ! equivqlent specific humidity at the surface
,avllm &            ! weight in balance eq., 0.<= avllm <=1.
,avlsf              ! weight in balance eq., 0.<= avlsf <=1.

real(kind=kfpt),intent(inout):: &
 akt(kmx) &         ! soil/ice/snow thermal diffusivity, middle of the layers
,tg(kmx)            ! soil/ice/snow temperature, midle of the layers

real(kind=kfpt),intent(out):: &
 pevap &            ! potential evaporation (m)
,evap &             ! total evaporation (m)
,eveg &             ! evapotranspiration (m)
,egrnd &            ! bare ground evaporation (m)
,smelt              ! melted snow (m)

real(kind=kfpt),intent(inout):: &
 tflux &            ! surface sensible heat flux, met. grad. sign
,gflux &            ! ground heat flux          , met. grad. sign
,qflux &            ! latent heat flux          , met. grad. sign
,pflux              ! potential evaporation flux, met. grad. sign
!--local variables------------------------------------------------------
integer(kind=kint):: &
 k                  ! working index, increases downwards

real(kind=kfpt):: &
 roa &              ! air density
,thlm &             ! potential temperature at the lowest atmospheric level
,xnrs &             ! temporary
,ts &               ! temporary
,qsrf &             ! basic spec. humidity for linearization around
,qss &              ! saturation specific humidity at the surface
,rdt &              ! 1./dt
,shfa &             !
,shfacp &           !
,shfael &           !
,elhlm &            !
,elhsf &            !
,trf3 &             !
,stefc &            !
,radot &            !
,radfc &            !
,cf &               !
,evfc &             !
,ts2 &              !
,ts4 &              !
,elhf2 &            !
,trf4

real(kind=kfpt):: &
 ti(kmx) &          ! buffer for initial tg
,rocdot(kmx) &      ! roc*dz/dt
,cm(kmx) &          ! central coeff, in the tridiagonal system
,akp(kmx) &         ! temporary
,cr(kmx) &          ! right coeff. in the tridiagonal system
,rst(kmx)           ! right hand side in the tridiagonal system
!--store initial tg in case of ice/sno melting--------------------------
if(sno.gt.epssno.or.wgs.le.epsw) then
  do k=2,kmp
    ti(k)=tg(k)
  enddo
endif
!--air variables--------------------------------------------------------
roa=plm/((qlm*0.608-clm+1.)*tlm*r)
thlm=tlm*rxnrlm
!--surface variables----------------------------------------------------
xnrs=1./rxnrs
ts=tg(1)
!--humidity to linearize around-----------------------------------------
if(trf.ne.tlm)  then
  qsrf=pq0lnd/ps*exp(a2*(trf-a3)/(trf-a4))
else
  qsrf=qas ! Penman for trf=tlm
endif
!-----------------------------------------------------------------------
rdt=1./dtphys
do k=1,kmp
  akt(k)=akt(k)*roc(k)
  rocdot(k)=roc(k)*dgt(k)*rdt
enddo

do k=2,kmp
  akp(k-1)=2.*akt(k-1)*akt(k)/(akt(k-1)*dgt(k)+akt(k)*dgt(k-1))
enddo
!--no evaporation if the air is already saturated-----------------------
if(qlm.ge.qas*tresh.and.qlm*avllm.lt.qs*avlsf) then
  avllm=0.
  avlsf=0.
endif
!--derived air variables------------------------------------------------
shfa=akhs*roa
shfacp=shfa*cp
shfael=shfa*elwv ! availabilities take care of fusion heat

elhlm=avllm*shfael
elhsf=avlsf*shfael
!-----------------------------------------------------------------------
elhf2=a23m4/(trf-a4)**2
trf3=trf*trf
trf3=trf3*trf
trf4=trf3*trf
stefc=epsr*stbol
radot=stefc*trf4*3.
radfc=stefc*trf3*4.
!--------------surface without ice/snow melting-------------------------
cm(1)=rocdot(1)+radfc+shfacp+elhsf*qsrf*elhf2+akp(1)
cr(1)=-akp(1)
rst(1)=rocdot(1)*ts+(radin+radot)+shfacp*thlm*xnrs &
      +elhlm*qlm-(1.-elhf2*trf)*elhsf*qsrf
!-----------------------------------------------------------------------
if(sice.lt.0.5) then ! land
!-----------------------------------------------------------------------
!  akp(kmp)=0. ! zero bottom heat flux
  akp(kmp)=akt(kmp)/dgt(kmp) ! fixed lower BC for temerature
  do k=2,kmp
    cf=-akp(k-1)/cm(k-1)
    cm(k)=-cr(k-1)*cf+akp(k-1)+akp(k)+rocdot(k)
    cr(k)=-akp(k)
    rst(k)=-rst(k-1)*cf+rocdot(k)*tg(k)
  enddo
!  tg(kmp)=rst(kmp)/cm(kmp) ! zero bottom heat flux
  tg(kmp)=(-cr(kmp)*tlbc+rst(kmp))/cm(kmp) ! fixed lower BC temerature
  do k=kmp-1,1,-1
    tg(k)=(-cr(k)*tg(k+1)+rst(k))/cm(k)
  enddo
!-----------------------------------------------------------------------
else ! sea ice
!-----------------------------------------------------------------------
  akp(kmp)=akt(kmp)*2./dgt(kmp) ! constant water temperature below ice
  do k=2,kmp
    cf=-akp(k-1)/cm(k-1)
    cm(k)=-cr(k-1)*cf+akp(k-1)+akp(k)+rocdot(k)
    cr(k)=-akp(k)
    rst(k)=-rst(k-1)*cf+rocdot(k)*tg(k)
  enddo
  tg(kmp)=(-cr(kmp)*tice+rst(kmp))/cm(kmp) ! sea water freezing point
  do k=kmp-1,1,-1
    tg(k)=(-cr(k)*tg(k+1)+rst(k))/cm(k)
  enddo
!-----------------------------------------------------------------------
endif ! end of land/sea ice branching
!-----------------------------------------------------------------------
ts=tg(1)
qss=pq0lnd/ps*exp(a2*(ts-a3)/(ts-a4))
evfc=shfa*dtphys/row
evap=(avllm*qlm-avlsf*qss)*evfc
eveg=(avglm*qlm-avgsf*qss)*evfc
egrnd=(agrlm*qlm-agrsf*qss)*evfc
pevap=(qlm-qss)*evfc
if(abs(pevap).gt.epsw) then
  qs=-evap/evfc+qlm ! equivalent specific humidity at the surface
  if(qs.lt.epsw) then
    qs=epsw
    evap=(qlm-qs)*evfc
  endif
else
  qs=qlm
  evap=0.
  eveg=0.
  egrnd=0.
  pevap=0.
endif
smelt=0.
!--------------surface with ice/snow melting----------------------------
if(tg(1).gt.t0.and.(sno.gt.epssno.or.wgs.le.epsw)) then
!-----------------------------------------------------------------------
  ts=t0
  tg(1)=t0
  qs=pq0lnd/ps
  qss=qs
  evap=(qlm-qs)*evfc
  eveg=0.
  egrnd=0.
  pevap=evap
!--------------recompute thermal diffusion with t0 at the surface-------
  do k=2,kmp
    tg(k)=ti(k)
  enddo
!-----------------------------------------------------------------------
  cm(1)=1.
  cr(1)=0.
  rst(1)=ts
!-----------------------------------------------------------------------
  if(sice.lt.0.5) then ! land
!-----------------------------------------------------------------------
!    akp(kmp)=0. ! zero bottom heat flux
    akp(kmp)=akt(kmp)/dgt(kmp) ! fixed lower BC for temerature
    do k=2,kmp
      cf=-akp(k-1)/cm(k-1)
      cm(k)=-cr(k-1)*cf+akp(k-1)+akp(k)+rocdot(k)
      cr(k)=-akp(k)
      rst(k)=-rst(k-1)*cf+rocdot(k)*tg(k)
    enddo
!    tg(kmp)=rst(kmp)/cm(kmp) ! zero bottom heat flux
    tg(kmp)=(-cr(kmp)*tlbc+rst(kmp))/cm(kmp) ! fixed lower BC temerature
    do k=kmp-1,2,-1
      tg(k)=(-cr(k)*tg(k+1)+rst(k))/cm(k)
    enddo
!-----------------------------------------------------------------------
  else ! sea ice
!-----------------------------------------------------------------------
    akp(kmp)=akt(kmp)*2./dgt(kmp) ! constant water temperature below ice
    do k=2,kmp
      cf=-akp(k-1)/cm(k-1)
      cm(k)=-cr(k-1)*cf+akp(k-1)+akp(k)+rocdot(k)
      cr(k)=-akp(k)
      rst(k)=-rst(k-1)*cf+rocdot(k)*tg(k)
    enddo
    tg(kmp)=(-cr(kmp)*tice+rst(kmp))/cm(kmp) ! sea water freezing point
    do k=kmp-1,2,-1
      tg(k)=(-cr(k)*tg(k+1)+rst(k))/cm(k)
    enddo
!-----------------------------------------------------------------------
  endif ! end of land/sea ice branching
!--------------melted snow amount---------------------------------------
  if(sno.gt.epssno) then
!-----------------------------------------------------------------------
    ts2=ts*ts
    ts4=ts2*ts2
    radot=stefc*ts4
!
    smelt=(radin-radot &
          +shfacp*(thlm*xnrs-ts)+elhlm*qlm-elhsf*qs+akp(1)*(ts-tg(2))) &
         *dtphys/rowliw
!-----------------------------------------------------------------------
  endif ! end of melted snow amount
!-----------------------------------------------------------------------
endif ! end of ice/snow melting case
!-----------------------------------------------------------------------
ths=ts*rxnrs
!--fluxes---------------------------------------------------------------
tflux=shfacp*(thlm*xnrs-tg(1))
gflux=akp(1)*(tg(1)-tg(2))
qflux=(elhlm*qlm-elhsf*qss)
if(ts.gt.t0) then
  pflux=shfa*elwv*(qlm-qss)
else
  pflux=shfa*eliv*(qlm-qss)
endif
!-----------------------------------------------------------------------
                        endsubroutine heat
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        subroutine water &
(dtphys,cveg,elai &
,wgs,wgpwp,be,gammas,akws &
,pevap,evap,eveg,smelt &
,dgw,rfr &
,prra,prsn,sno,roff &
,wliq,wg,wgi,wgl)
!-----------------------------------------------------------------------
! multiple layer soil moisture conduction
! z. janjic, june 2002, dwd
!-----------------------------------------------------------------------
implicit none
!--global variables-----------------------------------------------------
real(kind=kfpt),intent(in):: &
 dtphys &           ! physics time step
,cveg &             ! vegetation fraction
,elai &             ! leaf area index
,wgs &              ! soil moisture at saturation
,wgpwp &            ! soil moisture at permanent wilting point
,be &               ! Clapp and Hornberger nondimensional exponent
,gammas &           ! hydraulic conductivity at saturation
,akws &             ! hydraulic diffusivity at saturation
,pevap &            ! potential evaporation (m)
,evap &             ! total evaporation (m)
,eveg &             ! evapotranspiration (m)
,smelt              ! melted snow (m)

real(kind=kfpt),intent(in):: &
 dgw(nwet) &        ! depths of soil layers
,rfr(nwet) &        ! root fraction
,wgi(nwet)          ! frozen water fraction

real(kind=kfpt),intent(inout):: &
 wliq &             ! water in interception reservoir (m)
,prra &             ! time step rain precipitation (m)
,prsn &             ! time step snow precipitation (m)
,sno &              ! sno accumulation (m)
,roff               ! surface runoff (m)

real(kind=kfpt),intent(inout):: &
 wg(nwet) &         ! soil moisture
,wgl(nwet)          ! liquid part of soil moisture
!--local variables------------------------------------------------------
integer(kind=kint):: &
 k                  ! working index

real(kind=kfpt):: &
 cliq &             ! liquid fraction
,eliq &             ! evaporation (m)
,etop &             ! flux at the top (m/s)
,wintcp &           ! precipitation interception (m)
,thru &             ! precipitation throughfall (m)
,gamma1 &           !
,droff &            !
,fill &             ! infiltration
,rdt &              ! 1/dtphys

,wlmx &             ! capacity of interception reservoir (m)
,wliqn &            ! temporary
,rwlmx &            ! temporary
,cf &               ! temporary
,wgbar &            ! temporary
,rwgbar &           ! temporary
,coff &             ! temporary
,rwg &              ! temporary
,rwgs &             ! temporary
,tfrk               ! frozen soil water fraction in layer k

real(kind=kfpt):: &
 akw(nwet) &        ! diffusion coefficients
,dgwodt(nwet) &     ! dgw/dt
,cm(nwet) &         ! central term
,akp(nwet) &        ! diffusivity contribution
,cr(nwet) &         ! right term
,rst(nwet) &        ! right hand side
,gamma(nwet) &      ! conductivity
,rex(nwet)          ! root extraction
!--modify snow pack for melting, evaporation and new snow---------------
if(sno.gt.epssno) then
  sno=sno-smelt+evap
  if(sno.lt.smelt) sno=0.
endif ! end of snow pack modification for melting and evap.
sno=sno+prsn
!--interception reservoir, interception and throughfall-----------------
if(sno.le.epssno.and.wgs.gt.epsw) then
  wlmx=(cveg*elai+(1.-cveg))*wlaymx
  rwlmx=1./wlmx
  cliq=min(wliq*rwlmx,1.)

  if(pevap.lt.0.) then ! evaporation
    wliqn=cliq*pevap/(1.-pevap*rwlmx)+wliq
  else ! dew deposition
    wliqn=cliq*pevap+wliq
  endif
  wliqn=max(wliqn,0.)

  eliq=wliqn-wliq
  wintcp=min(eintcp*cveg*prra,wlmx-wliqn)
  wliq=wliqn+wintcp
  thru=max(prra-wintcp,0.)
else
  eliq=0.
  wliq=0.
  thru=0.
endif

prra=0.
prsn=0.
!-----------------------------------------------------------------------
if(wgs.gt.epsw) then ! normal soil that takes water
!--diffusivity, conductivity, root extraction---------------------------
  rdt=1./dtphys
  wgbar=0.
  do k=1,nwet
    wgbar=wgl(k)*rfr(k)+wgbar
  enddo
  if(wgbar.gt.epsw) then
    rwgbar=1./wgbar
  else
    rwgbar=0.
  endif
  coff=eveg*rdt*rwgbar

  do k=1,nwet
    if(wgs.gt.epsw) then
      rwgs=1./wgs
      if(wgl(k).ge.epsw) then ! liquid water part
        rwg=wgl(k)*rwgs
        gamma(k)=gammas*rwg**(2.*be+3.)
        akw(k)=akws*rwg**(be+2.)
      else
        gamma(k)=0.
        akw(k)=0.
      endif
      if(wgi(k).ge.epsw) then ! frozen water part, ecmwf style fix
        rwg=wgpwp*rwgs
        tfrk=wgi(k)/wg(k)
        gamma(k)=(1.-tfrk)*gamma(k)+tfrk*gammas*rwg**(2.*be+3.)
        akw(k)=(1.-tfrk)*akw(k)+tfrk*akws*rwg**(be+2.)
      endif
    endif
    rex(k)=min(wgl(k)*rfr(k)*coff,0.)
  enddo
  gamma1=gamma(1)

  dgwodt(1)=dgw(1)*rdt
  do k=2,nwet
    dgwodt(k)=dgw(k)*rdt
    if(akw(k-1).gt.epsw.and.akw(k).gt.epsw) then
      akp(k-1)=2.*akw(k-1)*akw(k)/(akw(k-1)*dgw(k)+akw(k)*dgw(k-1))
    else
      akp(k-1)=0.
    endif
    gamma(k-1)=(gamma(k-1)*dgw(k-1)+gamma(k)*dgw(k))/(dgw(k-1)+dgw(k))
  enddo
  akp(nwet)=0.
!--infiltration and surface runoff--------------------------------------
  fill=(2.*akw(1)/dgw(1)*(wgs-wg(1))+gamma1)*dtphys
  if(fill.gt.0.) then
    droff=max(thru+smelt-fill*dgw(1),0.)
    roff=roff+droff
  else
    droff=0.
  endif
!--upper boundary condition---------------------------------------------
  if(sno.le.epssno) then ! no sno
    if(eveg.le.0.) then ! evaporation into atmosphere
      etop=(evap-eliq-eveg+thru+smelt-droff)*rdt
    else ! dew deposition
      etop=(evap-eliq+thru+smelt-droff)*rdt
    endif
  else ! sno
    etop=(smelt-droff)*rdt
  endif
!--tridiagonal system solver, only liquid water-------------------------
  cm(1)=akp(1)+dgwodt(1)
  cr(1)=-akp(1)
  rst(1)=dgwodt(1)*wgl(1)-gamma(1)+dgw(1)*rex(1)+etop

  do k=2,nwet
    cf=-akp(k-1)/cm(k-1)
    cm(k)=-cr(k-1)*cf+akp(k-1)+akp(k)+dgwodt(k)
    cr(k)=-akp(k)
    rst(k)=-rst(k-1)*cf+dgwodt(k)*wgl(k)  &
           +gamma(k-1)-gamma(k)+dgw(k)*rex(k)
  enddo
  wgl(nwet)=min(max(rst(nwet)/cm(nwet),epsw),wgs-wgi(nwet))
  do k=nwet-1,1,-1
    wgl(k)=min(max((-cr(k)*wgl(k+1)+rst(k))/cm(k),epsw),wgs-wgi(k))
  enddo
  do k=1,nwet
    wg(k)=wgl(k)+wgi(k)
  enddo
!--tridiagonal system solver, all moisture------------------------------
!  cm(1)=akp(1)+dgwodt(1)
!  cr(1)=-akp(1)
!  rst(1)=dgwodt(1)*wg(1)-gamma(1)+dgw(1)*rex(1)+etop
!
!  do k=2,nwet
!    cf=-akp(k-1)/cm(k-1)
!    cm(k)=-cr(k-1)*cf+akp(k-1)+akp(k)+dgwodt(k)
!    cr(k)=-akp(k)
!    rst(k)=-rst(k-1)*cf+dgwodt(k)*wg(k)  &
!           +gamma(k-1)-gamma(k)+dgw(k)*rex(k)
!  enddo
!  wg(nwet)=min(max(rst(nwet)/cm(nwet),epsw),wgs)
!  do k=nwet-1,1,-1
!    wg(k)=min(max((-cr(k)*wg(k+1)+rst(k))/cm(k),epsw),wgs)
!  enddo
!-----------------------------------------------------------------------
else ! impermeable surface, such as ice, pavement
!-----------------------------------------------------------------------
  droff=max(thru+smelt+evap,0.)
  roff=roff+droff
!-----------------------------------------------------------------------
endif ! end of soil/rock branch
!-----------------------------------------------------------------------
                        endsubroutine water
!-----------------------------------------------------------------------
                        end module module_surface

