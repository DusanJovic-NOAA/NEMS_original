!-----------------------------------------------------------------------
!
      MODULE MODULE_SURFACE_LAYER
!
!-----------------------------------------------------------------------
!
!***  THE SURFACE LAYER PACKAGES AND THE SURFACE DRIVER (FOR THE
!***  SURFACE LAYER AND LANDSURFACE SCHEMES).
!
!-----------------------------------------------------------------------
!
      USE MODULE_INCLUDE
!
!     USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
!                                  ,IMS,IME,JMS,JME                     &
!                                  ,ITS,ITE,JTS,JTE                     &
      USE MODULE_DM_PARALLEL,ONLY : MPI_COMM_COMP,MYPE_SHARE,NUM_TILES
!
      USE MODULE_LANDSURFACE
!
      USE MODULE_CONSTANTS,ONLY: A2,A3,A4,CAPPA,CP,ELWV,EPSQ2,G         &
                                ,P608,PQ0,R_D
      USE MODULE_PHYSICS_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: MYJSFC,MYJSFC_INIT,SURFACE_DRIVER
!
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE
!
      REAL,PARAMETER :: RCP=CAPPA,XLV=ELWV
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE SURFACE LAYER OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),PARAMETER :: MYJSFCSCHEME=2                    &
                                     ,GFSSFCSCHEME=3                    &
                                     ,SFCLAYSCHEME=1
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  THE LANDSURFACE OPTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER(KIND=KINT),PARAMETER :: NMMLSMSCHEME=99                   &
                                     ,SLABSCHEME  =1                    &
                                     ,LSMSCHEME   =2                    &
                                     ,RUCLSMSCHEME=3                    &  
                                     ,lissscheme  =101
!
      INTEGER(KIND=KINT),PARAMETER :: JF = 1000000
      INTEGER(KIND=KINT),PARAMETER :: ILIM=97                           &
                                     ,JLIM=101                          &
                                     ,MAXLEV=39
      INTEGER(KIND=KINT),PARAMETER :: ITOT=ILIM*JLIM
      INTEGER(KIND=KINT),PARAMETER :: MBUF=2000000
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  FOR MYJSFC SCHEME
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      INTEGER :: ITRMX=5 ! Iteration count for sfc layer computations
!
      REAL,PARAMETER :: VKARMAN=0.4
      REAL,PARAMETER :: CAPA=R_D/CP,ELOCP=2.72E6/CP,RCAP=1./CAPA
      REAL,PARAMETER :: GOCP02=G/CP*2.,GOCP10=G/CP*10.
!        REAL,PARAMETER :: EPSU2=1.E-4,EPSUST=0.07,EPSZT=1.E-5
!       ECMWF sets lower limit on ln(Z0h)=-20 --> EPSZT=2.e-9
      REAL,PARAMETER :: EPSU2=1.E-6,EPSUST=1.e-9,EPSZT=1.E-28
      REAL,PARAMETER :: A2S=17.2693882,A3S=273.16,A4S=35.86
      REAL,PARAMETER :: SEAFC=0.98,PQ0SEA=PQ0*SEAFC
      REAL,PARAMETER :: BETA=1./273.,CZIL=0.1,EXCML=0.0001,EXCMS=0.0001  &
     &                 ,GLKBR=10.,GLKBS=30.,PI=3.1415926               &
     &                 ,QVISC=2.1E-5,RIC=0.505,SMALL=0.35              &
     &                 ,SQPR=0.84,SQSC=0.84,SQVISC=258.2,TVISC=2.1E-5  &
     &                 ,USTC=0.7,USTR=0.225,VISC=1.5E-5,FH=1.10        &
     &                 ,WWST=1.2,ZTFC=1.,TOPOFAC=9.0e-6
!
      REAL,PARAMETER :: BTG=BETA*G,CZIV=SMALL*GLKBS                    &
     &                 ,GRRS=GLKBR/GLKBS                               &
     &                 ,RTVISC=1./TVISC,RVISC=1./VISC                  &
     &                 ,ZQRZT=SQSC/SQPR                                &
     &                 ,FZQ1=RTVISC*QVISC*ZQRZT                        &
     &                 ,FZQ2=RTVISC*QVISC*ZQRZT                        &
     &                 ,FZT1=RVISC *TVISC*SQPR                         &
     &                 ,FZT2=CZIV*GRRS*TVISC*SQPR                      &
     &                 ,FZU1=CZIV*VISC                                 &
     &                 ,PIHF=0.5*PI                                    &
     &                 ,RQVISC=1./QVISC                                &
     &                 ,USTFC=0.018/G                                  &
     &                 ,WWST2=WWST*WWST                                &
     &                 ,ZILFC=-CZIL*VKARMAN*SQVISC
!
!----------------------------------------------------------------------
      INTEGER, PARAMETER :: KZTM=10001,KZTM2=KZTM-2
!
      REAL :: DZETA1,DZETA2,FH01,FH02,ZTMAX1,ZTMAX2,ZTMIN1,ZTMIN2
!
      REAL,DIMENSION(KZTM) :: PSIH1,PSIH2,PSIM1,PSIM2
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
   SUBROUTINE surface_driver(                                         &
     &           acsnom,acsnow,akhs,akms,albedo,br,canwat             &
     &          ,chklowq,dt,dx,dz8w,dzs,glw                           &
     &          ,grdflx,gsw,swdown,gz1oz0,hfx,ht,ifsnow,isfflx        &
     &          ,isltyp,itimestep,ivgtyp,lowlyr,mavail,rmol           &
     &          ,num_soil_layers,p8w,pblh,pi_phy,pshltr,psih          &
     &          ,psim,p_phy,q10,q2,qfx,qsfc,qshltr,qz0                &
     &          ,raincv,rho,sfcevp,sfcexc,sfcrunoff                   &
     &          ,smois,smstav,smstot,snoalb,snow,snowc,snowh,stepbl   &
     &          ,th10,th2,thz0,th_phy,tmn,tshltr,tsk,tslb             &
     &          ,t_phy,u10,udrunoff,ust,uz0,u_frame,u_phy,v10,vegfra  &
     &          ,vz0,v_frame,v_phy,warm_rain,wspd,xice,xland,z,znt,zs &
     &          ,ct,tke_myj                                           &
     &          ,albbck,lh,sh2o,shdmax,shdmin,z0                      &
     &          ,flqc,flhc,psfc,sst,sst_update,t2,emiss                        &
     &          ,sf_sfclay_physics,sf_surface_physics,ra_lw_physics   &
            !  Optional urban 
     &          ,declin_urb,cosz_urb2d,omg_urb2d,xlat_urb2d           & !I urban
     &          ,num_roof_layers, num_wall_layers                     & !I urban
     &          ,num_road_layers, dzr, dzb, dzg                       & !I urban
     &          ,tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d         & !H urban
     &          ,uc_urb2d                                             & !H urban
     &          ,xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d          & !H urban
     &          ,trl_urb3d,tbl_urb3d,tgl_urb3d                        & !H urban
     &          ,sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d,ts_urb2d          & !H urban
     &          ,frc_urb2d, utype_urb2d                               & !H urban
     &          ,ucmcall                                              & ! urban
     &          , ids,ide,jds,jde,kds,kde                             &
     &          , ims,ime,jms,jme,kms,kme                             &
     &          , i_start,i_end,j_start,j_end,kts,kte,num_tiles       &
             !  Optional moisture tracers
     &           ,qv_curr, qc_curr, qr_curr                           &
     &           ,qi_curr, qs_curr, qg_curr                           &
             !  Optional moisture tracer flags
     &           ,f_qv,f_qc,f_qr                                      &
     &           ,f_qi,f_qs,f_qg                                      &
             !  Other optionals (more or less em specific)
     &          ,capg,hol,mol                                   &
     &          ,rainncv,rainbl,regime,thc                         &
     &          ,qsg,qvg,qcg,soilt1,tsnav                             &
             !  Other optionals (more or less nmm specific)
     &          ,potevp,snopcx,soiltb,sr                              &
             !  Optional observation nudging
     &          ,uratx,vratx,tratx                                    &
                                                                      )

!nmmb#if ( ! NMM_CORE == 1 )
!nmmb  USE module_state_description, ONLY : SFCLAYSCHEME              &
!nmmb                                      ,MYJSFCSCHEME              &
!nmmb                                      ,GFSSFCSCHEME              &
!nmmb                                      ,SLABSCHEME                &
!nmmb                                      ,LSMSCHEME                 &
!nmmb                                      ,RUCLSMSCHEME
!nmmb#else
!nmmb  USE module_state_description, ONLY : SFCLAYSCHEME              &
!nmmb                                      ,MYJSFCSCHEME              &
!nmmb                                      ,GFSSFCSCHEME              &
!nmmb                                      ,SLABSCHEME                &
!nmmb                                      ,NMMLSMSCHEME              &
!nmmb                                      ,LSMSCHEME                 &
!nmmb                                      ,RUCLSMSCHEME
!nmmb#endif
!nmmb   USE module_model_constants
! *** add new modules of schemes here

!nmmb   USE module_sf_sfclay
!nmmb   USE module_sf_myjsfc
!nmmb   USE module_sf_gfs
!nmmb   USE module_sf_noahlsm
!nmmb   USE module_sf_ruclsm
!nmmb#if ( NMM_CORE == 1 )
!nmmb   USE module_sf_lsm_nmm
!nmmb#endif

!nmmb   USE module_sf_slab
!
!nmmb   USE module_sf_sfcdiags
!

   !  This driver calls subroutines for the surface parameterizations.
   !
   !  surface layer: (between surface and pbl)
   !      1. sfclay
   !      2. myjsfc
   !  surface: ground temp/lsm scheme:
   !      1. slab
   !      2. Noah LSM
   !      99. NMM LSM (NMM core only)
!------------------------------------------------------------------
   IMPLICIT NONE
!======================================================================
! Grid structure in physics part of WRF
!----------------------------------------------------------------------
! The horizontal velocities used in the physics are unstaggered
! relative to temperature/moisture variables. All predicted
! variables are carried at half levels except w, which is at full
! levels. Some arrays with names (*8w) are at w (full) levels.
!
!----------------------------------------------------------------------
! In WRF, kms (smallest number) is the bottom level and kme (largest
! number) is the top level.  In your scheme, if 1 is at the top level,
! then you have to reverse the order in the k direction.
!
!         kme      -   half level (no data at this level)
!         kme    ----- full level
!         kme-1    -   half level
!         kme-1  ----- full level
!         .
!         kms+2    -   half level
!         kms+2  ----- full level
!         kms+1    -   half level
!         kms+1  ----- full level
!         kms      -   half level
!         kms    ----- full level
!
!======================================================================
! Definitions
!-----------
! Theta      potential temperature (K)
! Qv         water vapor mixing ratio (kg/kg)
! Qc         cloud water mixing ratio (kg/kg)
! Qr         rain water mixing ratio (kg/kg)
! Qi         cloud ice mixing ratio (kg/kg)
! Qs         snow mixing ratio (kg/kg)
!-----------------------------------------------------------------
!-- itimestep     number of time steps
!-- GLW           downward long wave flux at ground surface (W/m^2)
!-- GSW           net short wave flux at ground surface (W/m^2)
!-- SWDOWN        downward short wave flux at ground surface (W/m^2)
!-- EMISS         surface emissivity (between 0 and 1)
!-- TSK           surface temperature (K)
!-- TMN           soil temperature at lower boundary (K)
!-- XLAND         land mask (1 for land, 2 for water)
!-- ZNT           time-varying roughness length (m)
!-- Z0            background roughness length (m)
!-- MAVAIL        surface moisture availability (between 0 and 1)
!-- UST           u* in similarity theory (m/s)
!-- MOL           T* (similarity theory) (K)
!-- HOL           PBL height over Monin-Obukhov length
!-- PBLH          PBL height (m)
!-- CAPG          heat capacity for soil (J/K/m^3)
!-- THC           thermal inertia (Cal/cm/K/s^0.5)
!-- SNOWC         flag indicating snow coverage (1 for snow cover)
!-- HFX           net upward heat flux at the surface (W/m^2)
!-- QFX           net upward moisture flux at the surface (kg/m^2/s)
!-- LH            net upward latent heat flux at surface (W/m^2)
!-- REGIME        flag indicating PBL regime (stable, unstable, etc.)
!-- tke_myj       turbulence kinetic energy from Mellor-Yamada-Janjic (MYJ) (m^2/s^2)
!-- akhs          sfc exchange coefficient of heat/moisture from MYJ
!-- akms          sfc exchange coefficient of momentum from MYJ
!-- thz0          potential temperature at roughness length (K)
!-- uz0           u wind component at roughness length (m/s)
!-- vz0           v wind component at roughness length (m/s)
!-- qsfc          specific humidity at lower boundary (kg/kg)
!-- uratx         ratio of u over u10 (Added for obs-nudging)
!-- vratx         ratio of v over v10 (Added for obs-nudging)
!-- tratx         ratio of t over th2 (Added for obs-nudging)
!-- u10           diagnostic 10-m u component from surface layer
!-- v10           diagnostic 10-m v component from surface layer
!-- th2           diagnostic 2-m theta from surface layer and lsm
!-- t2            diagnostic 2-m temperature from surface layer and lsm
!-- q2            diagnostic 2-m mixing ratio from surface layer and lsm
!-- tshltr        diagnostic 2-m theta from MYJ
!-- th10          diagnostic 10-m theta from MYJ
!-- qshltr        diagnostic 2-m specific humidity from MYJ
!-- q10           diagnostic 10-m specific humidity from MYJ
!-- lowlyr        index of lowest model layer above ground
!-- rr            dry air density (kg/m^3)
!-- u_phy         u-velocity interpolated to theta points (m/s)
!-- v_phy         v-velocity interpolated to theta points (m/s)
!-- th_phy        potential temperature (K)
!-- moist         moisture array (4D - last index is species) (kg/kg)
!-- p_phy         pressure (Pa)
!-- pi_phy        exner function (dimensionless)
!-- pshltr        diagnostic shelter (2m) pressure from MYJ (Pa)
!-- p8w           pressure at full levels (Pa)
!-- t_phy         temperature (K)
!-- dz8w          dz between full levels (m)
!-- z             height above sea level (m)
!-- DX            horizontal space interval (m)
!-- DT            time step (second)
!-- PSFC          pressure at the surface (Pa)
!-- SST           sea-surface temperature (K)
!-- TSLB          
!-- ZS
!-- DZS
!-- num_soil_layers number of soil layer
!-- IFSNOW      ifsnow=1 for snow-cover effects
!
!-- ids           start index for i in domain
!-- ide           end index for i in domain
!-- jds           start index for j in domain
!-- jde           end index for j in domain
!-- kds           start index for k in domain
!-- kde           end index for k in domain
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-- its           start index for i in tile
!-- ite           end index for i in tile
!-- jts           start index for j in tile
!-- jte           end index for j in tile
!-- kts           start index for k in tile
!-- kte           end index for k in tile
!
!******************************************************************
!------------------------------------------------------------------ 

   INTEGER, INTENT(IN) ::                                             &
     &           ids,ide,jds,jde,kds,kde                              &
     &          ,ims,ime,jms,jme,kms,kme                              &
     &          ,kts,kte,num_tiles

   INTEGER, INTENT(IN) :: sf_sfclay_physics,sf_surface_physics,ra_lw_physics,sst_update

   INTEGER, INTENT(IN) :: ucmcall                                     !urban

   INTEGER, DIMENSION(num_tiles), INTENT(IN) ::                       &
     &           i_start,i_end,j_start,j_end

   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   ISLTYP
   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   IVGTYP
   INTEGER, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   LOWLYR
   INTEGER, INTENT(IN )::   IFSNOW
   INTEGER, INTENT(IN )::   ISFFLX
   INTEGER, INTENT(IN )::   ITIMESTEP
   INTEGER, INTENT(IN )::   NUM_SOIL_LAYERS
   INTEGER, INTENT(IN )::   STEPBL
   LOGICAL, INTENT(IN )::   WARM_RAIN
   REAL , INTENT(IN )::   U_FRAME
   REAL , INTENT(IN )::   V_FRAME
   REAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), INTENT(INOUT)::   SMOIS
   REAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), INTENT(INOUT)::   TSLB
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   GLW
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   GSW,SWDOWN
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   HT
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   RAINCV
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   SST
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   TMN
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   VEGFRA
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   XICE
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   XLAND
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT)::   MAVAIL
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT)::   SNOALB
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   ACSNOW
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   AKHS
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   AKMS
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   ALBEDO
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   CANWAT


   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   GRDFLX
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   HFX
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   RMOL
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   PBLH
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   Q2
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   QFX
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   QSFC
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   QZ0
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SFCRUNOFF
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SMSTAV
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SMSTOT
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SNOW
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SNOWC
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   SNOWH
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   TH2
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   THZ0
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   TSK
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   UDRUNOFF
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   UST
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   UZ0
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   VZ0
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   WSPD
   REAL, DIMENSION( ims:ime, jms:jme ) , INTENT(INOUT)::   ZNT
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   BR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   CHKLOWQ
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   GZ1OZ0
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   PSHLTR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   PSIH
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   PSIM
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   Q10
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   QSHLTR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   TH10
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   TSHLTR
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   U10
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   V10
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(OUT)::   PSFC
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   ACSNOM
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   SFCEVP
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   SFCEXC
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   FLHC
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)::   FLQC
   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::   CT
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   DZ8W
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   P8W
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   PI_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   P_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   RHO
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   TH_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   T_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   U_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   V_PHY
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN )::   Z
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(INOUT) ::   TKE_MYJ
   REAL, DIMENSION(1:num_soil_layers), INTENT(IN)::   DZS
   REAL, DIMENSION(1:num_soil_layers), INTENT(IN)::   ZS
   REAL, INTENT(IN )::   DT
   REAL, INTENT(IN )::   DX

!  arguments for NCAR surface physics

   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT )::   ALBBCK  ! INOUT needed for NMM
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(INOUT )::   LH
   REAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), INTENT(INOUT)::   SH2O
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   SHDMAX
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   SHDMIN
   REAL, DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   Z0

!
! Optional
!

!
!  Observation nudging
!
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(OUT)::   uratx  !Added for obs-nudging
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(OUT)::   vratx  !Added for obs-nudging
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(OUT)::   tratx  !Added for obs-nudging
!
! Flags relating to the optional tendency arrays declared above
! Models that carry the optional tendencies will provdide the
! optional arguments at compile time; these flags all the model
! to determine at run-time whether a particular tracer is in
! use or not.
!
   LOGICAL, INTENT(IN), OPTIONAL ::                             &
                                                      f_qv      &
                                                     ,f_qc      &
                                                     ,f_qr      &
                                                     ,f_qi      &
                                                     ,f_qs      &
                                                     ,f_qg

   REAL, DIMENSION( ims:ime, kms:kme, jms:jme ),                 &
         OPTIONAL, INTENT(INOUT) ::                              &
                      ! optional moisture tracers
                      ! 2 time levels; if only one then use CURR
                      qv_curr, qc_curr, qr_curr                  &
                     ,qi_curr, qs_curr, qg_curr

   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   capg
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   emiss
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   hol
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   mol
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   regime
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(IN )::     rainncv
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   RAINBL
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   t2
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(IN )::     thc
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   qsg
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   qvg
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   qcg
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   soilt1
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   tsnav
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   potevp ! NMM LSM
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   snopcx ! NMM LSM
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   soiltb ! NMM LSM
   REAL, DIMENSION( ims:ime, jms:jme ), OPTIONAL, INTENT(INOUT)::   sr ! NMM and RUC LSM

!  LOCAL  VAR

   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ) ::v_phytmp
   REAL,       DIMENSION( ims:ime, kms:kme, jms:jme ) ::u_phytmp

   REAL,       DIMENSION( ims:ime, jms:jme )          ::  ZOL

   REAL,       DIMENSION( ims:ime, jms:jme )          ::          &
                                                             QGH, &
                                                             CHS, &
                                                             CPM, &
                                                            CHS2, &
                                                            CQS2

   REAL    :: DTMIN,DTBL
!
   INTEGER :: i,J,K,NK,jj,ij
   LOGICAL :: radiation, myj, frpcpn
!-------------------------------------------------
! urban related variables are added to declaration
!-------------------------------------------------
     REAL, OPTIONAL, INTENT(IN) :: DECLIN_URB                                 !urban
     REAL, OPTIONAL , DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: COSZ_URB2D  !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: OMG_URB2D   !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(IN) :: XLAT_URB2D  !urban
     INTEGER,  OPTIONAL, INTENT(IN) :: num_roof_layers                         !urban
     INTEGER,  OPTIONAL, INTENT(IN) :: num_wall_layers                         !urban
     INTEGER,  OPTIONAL, INTENT(IN) :: num_road_layers                         !urban
     REAL, OPTIONAL, DIMENSION(1:num_soil_layers), INTENT(IN) :: DZR          !urban
     REAL, OPTIONAL, DIMENSION(1:num_soil_layers), INTENT(IN) :: DZB          !urban
     REAL, OPTIONAL, DIMENSION(1:num_soil_layers), INTENT(IN) :: DZG          !urban

     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: TR_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: TB_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: TG_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: TC_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: QC_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: UC_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: XXXR_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: XXXB_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: XXXG_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: XXXC_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), &       !urban
           INTENT(INOUT)  :: TRL_URB3D                                 !urban
     REAL, OPTIONAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), &       !urban
           INTENT(INOUT)  :: TBL_URB3D                                 !urban
     REAL, OPTIONAL, DIMENSION( ims:ime , 1:num_soil_layers, jms:jme ), &       !urban
           INTENT(INOUT)  :: TGL_URB3D                                 !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: SH_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: LH_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: G_URB2D  !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: RN_URB2D !urban
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT):: TS_URB2D !urban
!
     REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: FRC_URB2D  !urban 
     INTEGER, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT)  :: UTYPE_URB2D  !urban 

     INTEGER :: NRL,NWL,NRDL,ihr, LUGB, LUGI,KF,IRET, NUMVAL,IGDNUM,LGRIB,LSKIP,IMAX
     INTEGER :: NUMLEV,JMAX,NNUM,NLEN,KR,IRETGI,IRETGB,II,KK,KSKIP 
     INTEGER :: IRETGB2, JR, IRGS, KMAX, FHR1, FHR2, IRGI, IRETGI2, ISTAT
     INTEGER :: JJINC, LUGI2, LUGB2, LUGB3, JJ1, LUGB5, LUGB4
     REAL :: F(JF)
     INTEGER, DIMENSION(200) :: JPDS, JGDS, KPDS, KGDS,JENS,KENS,IV,L 
     LOGICAL :: LB 
     REAL,  DIMENSION(ITOT) :: GRID,APCP,APCP2,CAPCP,CAPCP2,APCP3HR,CAPCP3HR 
     LOGICAL, DIMENSION(ITOT) :: MASK, MASK2
     CHARACTER, DIMENSION(MBUF) :: CBUF, CBUF2
     CHARACTER(11) :: ENVVAR
     CHARACTER(80) :: FNAME
     
     REAL,  DIMENSION( ims:ime, jms:jme )  :: PSIM_URB2D  !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: PSIH_URB2D  !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: GZ1OZ0_URB2D  !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: AKMS_URB2D  !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: U10_URB2D   !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: V10_URB2D   !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: TH2_URB2D   !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: Q2_URB2D    !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: UST_URB2D  !urban local var
     REAL,  DIMENSION( ims:ime, jms:jme )  :: RIB   ! Bulk Richardson Number

!------------------------------------------------------------------
   CHARACTER*256 :: message
   CHARACTER(64) :: infile
!------------------------------------------------------------------
!

  if (sf_sfclay_physics .eq. 0) return

  v_phytmp = 0.
  u_phytmp = 0.
  ZOL = 0.
  QGH = 0.
  CHS = 0.
  CPM = 0.
  CHS2 = 0.
  DTMIN = 0.
  DTBL = 0.

!  RAINBL in mm (Accumulation between PBL calls)
  IF ( PRESENT( rainncv ) .AND. PRESENT( rainbl ) ) THEN
!jaa    !$omp parallel do   &
!jaa    !$omp private ( ij, i, j, k )
    DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         RAINBL(i,j) = RAINBL(i,j) + RAINCV(i,j) + RAINNCV(i,j) 
         RAINBL(i,j) = MAX (RAINBL(i,j), 0.0) 
      ENDDO
      ENDDO
    ENDDO
!jaa     !$omp end parallel do
  ELSE IF ( PRESENT( rainbl ) ) THEN
!jaa     !$omp parallel do   &
!jaa     !$omp private ( ij, i, j, k )
    DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
         RAINBL(i,j) = RAINBL(i,j) + RAINCV(i,j)
         RAINBL(i,j) = MAX (RAINBL(i,j), 0.0)
      ENDDO
      ENDDO
    ENDDO
!jaa    !$omp end parallel do
  ENDIF
  
! Update SST
  IF (sst_update .EQ. 1) THEN
!jaa    !$omp parallel do   &
!jaa    !$omp private ( ij, i, j, k )
    DO ij = 1 , num_tiles
      DO j=j_start(ij),j_end(ij)
      DO i=i_start(ij),i_end(ij)
        IF(XLAND(i,j) .GT. 1.5)TSK(i,j)=SST(i,j)
      ENDDO
      ENDDO
    ENDDO
!jaa     !$omp end parallel do
  ENDIF

  IF (itimestep .eq. 1 .or. mod(itimestep,STEPBL) .eq. 0) THEN

  radiation = .false.
  myj = .false.
  frpcpn = .false.

  IF (ra_lw_physics .gt. 0) radiation = .true.

!---- 
! CALCULATE CONSTANT
 
     DTMIN=DT/60.
! Surface schemes need PBL time step for updates and accumulations
! Assume these schemes provide no tendencies
     DTBL=DT*STEPBL

! SAVE OLD VALUES


!jaa      !$omp parallel do   &
!jaa      !$omp private ( ij, i, j, k )
     DO ij = 1 , num_tiles
       DO j=j_start(ij),j_end(ij)
       DO i=i_start(ij),i_end(ij)
! PSFC : in Pa
          PSFC(I,J)=p8w(I,kts,J)
! REVERSE ORDER IN THE VERTICAL DIRECTION
          DO k=kts,kte
            v_phytmp(i,k,j)=v_phy(i,k,j)+v_frame
            u_phytmp(i,k,j)=u_phy(i,k,j)+u_frame
          ENDDO
       ENDDO
       ENDDO
     ENDDO
!jaa      !$omp end parallel do

!jaa      !$omp parallel do   &
 !jaa    !$omp private ( ij, i, j, k )
      DO ij = 1 , num_tiles
     sfclay_select: SELECT CASE(sf_sfclay_physics)

     CASE (SFCLAYSCHEME)
!nmmb#if (NMM_CORE != 1)
! DX varies spatially in NMM, therefore, SFCLAY cannot be called
! because it takes a scalar DX. NMM passes in a dummy value for this
! scalar.  NEEDS FURTHER ATTENTION. JM 20050215
!nmmb   IF (PRESENT(qv_curr)                            .AND.    &
!nmmb       PRESENT(mol)        .AND.  PRESENT(regime)  .AND.    &
!nmmb                                                  .TRUE. ) THEN
!nmmb     CALL wrf_debug( 100, 'in SFCLAY' )
!nmmb     CALL SFCLAY(u_phytmp,v_phytmp,t_phy,qv_curr,&
!nmmb           p_phy,dz8w,cp,g,rcp,r_d,xlv,psfc,chs,chs2,cqs2,cpm, &
!nmmb           znt,ust,pblh,mavail,zol,mol,regime,psim,psih,       &
!nmmb           xland,hfx,qfx,lh,tsk,flhc,flqc,qgh,qsfc,rmol,       &
!nmmb           uratx,vratx,tratx,                                  &
!nmmb           u10,v10,th2,t2,q2,                                  &
!nmmb           gz1oz0,wspd,br,isfflx,dx,                           &
!nmmb           svp1,svp2,svp3,svpt0,ep_1,ep_2,karman,eomeg,stbolt, &
!nmmb           ids,ide, jds,jde, kds,kde,                          &
!nmmb           ims,ime, jms,jme, kms,kme,                          &
!nmmb           i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )
!nmmb   ELSE
!nmmb     CALL wrf_error_fatal('Lacking arguments for SFCLAY in surface driver')
!nmmb   ENDIF

!nmmb#else
!nmmb   CALL wrf_error_fatal('SFCLAY cannot be used with NMM')
!nmmb#endif
      CASE (MYJSFCSCHEME)
       IF (PRESENT(qv_curr)    .AND.  PRESENT(qc_curr) .AND.    &
                                                      .TRUE. ) THEN

        myj =.true.

!nmmb        CALL wrf_debug(100,'in MYJSFC')

            CALL MYJSFC(itimestep,ht,dz8w,                         &
              p_phy,p8w,th_phy,t_phy,                              &
              qv_curr,qc_curr,                                      &
              u_phy,v_phy,tke_myj,                                 &
              tsk,qsfc,thz0,qz0,uz0,vz0,                           &
              lowlyr,                                              &
              xland,                                               &
              ust,znt,z0,pblh,mavail,rmol,                         &
              akhs,akms,                                           &
              chs,chs2,cqs2,hfx,qfx,lh,flhc,flqc,qgh,cpm,ct,       &
              u10,v10,t2,th2,tshltr,th10,q2,qshltr,q10,pshltr,rib, &  ! Added Bulk Richardson No.
              ids,ide, jds,jde, kds,kde,                           &
              ims,ime, jms,jme, kms,kme,                           &
              i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )
!nmmb   ELSE
!nmmb     CALL wrf_error_fatal('Lacking arguments for MYJSFC in surface driver')
       ENDIF

     CASE (GFSSFCSCHEME)
!nmmb   IF (PRESENT(qv_curr) .AND. .TRUE. ) THEN
!nmmb   CALL wrf_debug( 100, 'in GFSSFC' )
!nmmb     CALL SF_GFS(u_phytmp,v_phytmp,t_phy,qv_curr,              &
!nmmb           p_phy,CP,RCP,R_d,XLV,PSFC,CHS,CHS2,CQS2,CPM,        &
!nmmb           ZNT,UST,PSIM,PSIH,                                  &
!nmmb           XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,                     &
!nmmb           QGH,QSFC,U10,V10,                                   &
!nmmb           GZ1OZ0,WSPD,BR,ISFFLX,                              &
!nmmb           EP_1,EP_2,KARMAN,itimestep,                         &
!nmmb           ids,ide, jds,jde, kds,kde,                          &
!nmmb           ims,ime, jms,jme, kms,kme,                          &
!nmmb           i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )
!nmmb    CALL wrf_debug(100,'in SFCDIAGS')
!nmmb   ELSE
!nmmb     CALL wrf_error_fatal('Lacking arguments for SF_GFS in surface driver')
!nmmb   ENDIF

     CASE DEFAULT
        
!nmmb   WRITE( message , * )                                &
       WRITE(0,*) &
   'The sfclay option does not exist: sf_sfclay_physics = ', sf_sfclay_physics
!nmmb   CALL wrf_error_fatal ( message )

     END SELECT sfclay_select
     ENDDO
!jaa      !$omp end parallel do

     IF (ISFFLX.EQ.0 ) GOTO 430
!jaa      !$omp parallel do   &
!jaa      !$omp private ( ij, i, j, k )
     DO ij = 1 , num_tiles

     sfc_select: SELECT CASE(sf_surface_physics)

     CASE (SLABSCHEME)

!nmmb   IF (PRESENT(qv_curr)                            .AND.    &
!nmmb       PRESENT(capg)        .AND.    &
!nmmb                                                  .TRUE. ) THEN
!nmmb       DO j=j_start(ij),j_end(ij)
!nmmb       DO i=i_start(ij),i_end(ij)
!!nmmb      CQS2 ACCOUNTS FOR MAVAIL FOR SFCDIAGS 2M Q
!nmmb          CQS2(I,J)= CQS2(I,J)*MAVAIL(I,J)
!nmmb       ENDDO
!nmmb       ENDDO

!nmmb    CALL wrf_debug(100,'in SLAB')
!nmmb      CALL SLAB(t_phy,qv_curr,p_phy,flhc,flqc,  &
!nmmb         psfc,xland,tmn,hfx,qfx,lh,tsk,qsfc,chklowq,          &
!nmmb         gsw,glw,capg,thc,snowc,emiss,mavail,                 &
!nmmb         dtbl,rcp,xlv,dtmin,ifsnow,                           &
!nmmb         svp1,svp2,svp3,svpt0,ep_2,karman,eomeg,stbolt,       &
!nmmb         tslb,zs,dzs,num_soil_layers,radiation,               &
!nmmb         ids,ide, jds,jde, kds,kde,                           &
!nmmb         ims,ime, jms,jme, kms,kme,                           &
!nmmb         i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )

!nmmb       DO j=j_start(ij),j_end(ij)
!nmmb       DO i=i_start(ij),i_end(ij)
!nmmb          SFCEVP(I,J)= SFCEVP(I,J) + QFX(I,J)*DTBL
!nmmb       ENDDO
!nmmb       ENDDO

!nmmb    CALL wrf_debug(100,'in SFCDIAGS')
!nmmb      CALL SFCDIAGS(hfx,qfx,tsk,qsfc,chs2,cqs2,t2,th2,q2,      &
!nmmb                 psfc,cp,r_d,rcp,                              &
!nmmb                 ids,ide, jds,jde, kds,kde,                    &
!nmmb                 ims,ime, jms,jme, kms,kme,                    &
!nmmb         i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )

!nmmb   ELSE
!nmmb     CALL wrf_error_fatal('Lacking arguments for SLAB in surface driver')
!nmmb   ENDIF

!nmmb#if ( NMM_CORE == 1 )
     CASE (lissscheme)
       IF (PRESENT(qv_curr)    .AND.  PRESENT(rainbl)  .AND.    &
           PRESENT(potevp)     .AND.  PRESENT(snopcx)  .AND.    &
           PRESENT(soiltb)     .AND.  PRESENT(sr)      .AND.    &
                                                      .TRUE. ) THEN
!nmmb      CALL wrf_debug(100,'in NMM LSM')
!       if (ims .lt. 5 .and. jms .lt. 5) then
!       write(0,*) 't_phy(1,KDE-1,1): ', t_phy(1,KDE-1,1)
!       write(0,*) 't_phy(1,1,1): ', t_phy(1,1,1)
!       endif


       CALL liss(dz8w,qv_curr,p8w,rho,                      &
            t_phy,th_phy,tsk,chs,                           &
            hfx,qfx,qgh,swdown,glw,lh,rmol,                 &
            smstav,smstot,sfcrunoff,                        &
            udrunoff,ivgtyp,isltyp,vegfra,sfcevp,potevp,    &
            grdflx,sfcexc,acsnow,acsnom,snopcx,             &
            albbck,tmn,xland,xice,qz0,                      &
            th2,q2,snowc,cqs2,qsfc,soiltb,chklowq,rainbl,   &
            num_soil_layers,dtbl,dzs,itimestep,             &
            smois,tslb,snow,canwat,cpm,rcp,sr,              &    !tslb
            albedo,snoalb,sh2o,snowh,                       &
            ids,ide, jds,jde, kds,kde,                      &
            ims,ime, jms,jme, kms,kme,                      &
            i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte )


!nmmb      CALL wrf_debug(100,'back from NMM LSM')
       ELSE
         WRITE(0,*)'Lacking arguments for liss in surface driver'
!nmmb     CALL wrf_error_fatal('Lacking arguments for NMMLSM in surface driver')
       ENDIF
!nmmb#endif

     CASE (LSMSCHEME)

       IF ( PRESENT(qv_curr) .AND. PRESENT(rainbl) ) THEN
!------------------------------------------------------------------
         IF( PRESENT(sr) ) THEN
           frpcpn=.true.
         ENDIF

         nrl=1
         IF ( PRESENT(num_roof_layers) ) THEN
             nrl=num_roof_layers
         ENDIF
         nwl=1
         IF ( PRESENT(num_wall_layers) ) THEN
             nwl=num_wall_layers
         ENDIF
         nrdl=1
         IF ( PRESENT(num_road_layers) ) THEN
             nrdl=num_road_layers
         ENDIF

!nmmb     CALL wrf_debug(100,'in NOAH LSM')
           CALL noahlsm(dz8w,qv_curr,p8w,t_phy,tsk,             &
                hfx,qfx,lh,grdflx,qgh,gsw,swdown,glw,smstav,smstot,    &
                sfcrunoff,udrunoff,ivgtyp,isltyp,vegfra,        &
                albedo,albbck,znt,z0, tmn,xland,xice, emiss,    &
                snowc,qsfc,rainbl,                              &
                num_soil_layers,dtbl,dzs,itimestep,             &
                smois,tslb,snow,canwat,                         &
                chs, chs2, cqs2, cpm,rcp,SR,chklowq,qz0,        &
                myj,frpcpn,                                     &
                sh2o,snowh,                                     & !h
                u_phy,v_phy,                                    & !I
                snoalb,shdmin,shdmax,                           & !i
                acsnom,acsnow,                                  & !o
                snopcx,                                         & !o
                potevp, rib,                                    & !o Added Bulk Richardson No.
                ids,ide, jds,jde, kds,kde,                      &
                ims,ime, jms,jme, kms,kme,                      &
                i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte,    &
                ucmcall,                                         &
!Optional urban
                tr_urb2d,tb_urb2d,tg_urb2d,tc_urb2d,qc_urb2d,  & !H urban
                uc_urb2d,                                       & !H urban
                xxxr_urb2d,xxxb_urb2d,xxxg_urb2d,xxxc_urb2d,    & !H urban
                trl_urb3d,tbl_urb3d,tgl_urb3d,                  & !H urban
                sh_urb2d,lh_urb2d,g_urb2d,rn_urb2d,ts_urb2d,    & !H urban
                psim_urb2d,psih_urb2d,u10_urb2d,v10_urb2d,      & !O urban
                GZ1OZ0_urb2d, AKMS_URB2D,                       & !O urban
                th2_urb2d,q2_urb2d,ust_urb2d,                   & !O urban
                declin_urb,cosz_urb2d,omg_urb2d,                & !I urban
                xlat_urb2d,                                     & !I urban
                nrl, nwl,                                       & !I urban
                nrdl, DZR, DZB, DZG,                            & !I urban
                FRC_URB2D, UTYPE_URB2D                          & ! urban
                )

           DO j=j_start(ij),j_end(ij)
           DO i=i_start(ij),i_end(ij)
               SFCEVP(I,J)= SFCEVP(I,J) + QFX(I,J)*DTBL
               SFCEXC(I,J)= CHS(I,J)
               SOILTB(I,J)= TSLB(I,num_soil_layers,J) !  nmmlsm vrbl., here only for output
           ENDDO
           ENDDO

          CALL SFCDIAGS(HFX,QFX,TSK,QSFC,CHS2,CQS2,T2,TH2,Q2,      &
                     PSFC,CP,R_d,RCP,                              &
                     ids,ide, jds,jde, kds,kde,                    &
                     ims,ime, jms,jme, kms,kme,                    &
             i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )


!urban
     IF(UCMCALL.eq.1) THEN
       DO j=j_start(ij),j_end(ij)                             !urban
         DO i=i_start(ij),i_end(ij)                           !urban
          IF( IVGTYP(I,J) == 1 .or. IVGTYP(I,J) == 31 .or. &  !urban
              IVGTYP(I,J) == 32 .or. IVGTYP(I,J) == 33 ) THEN !urban
             T2(I,J)   = FRC_URB2D(i,j)*TH2_URB2D(I,J) + (1-FRC_URB2D(i,j))*T2(I,J) !urban
             TH2(I,J) = T2(I,J)*(1.E5/PSFC(I,J))**RCP                               !urban
             Q2(I,J)   = FRC_URB2D(i,j)*Q2_URB2D(I,J) +(1-FRC_URB2D(i,j))* Q2(I,J)  !urban
             U10(I,J)  = U10_URB2D(I,J)                       !urban
             V10(I,J)  = V10_URB2D(I,J)                       !urban
             PSIM(I,J) = PSIM_URB2D(I,J)                      !urban
             PSIH(I,J) = PSIH_URB2D(I,J)                      !urban
             GZ1OZ0(I,J) = GZ1OZ0_URB2D(I,J)                  !urban
             AKHS(I,J) = CHS(I,J)                             !urban
             AKMS(I,J) = AKMS_URB2D(I,J)                      !urban
           END IF                                             !urban
         ENDDO                                                !urban
       ENDDO                                                  !urban
     ENDIF
!------------------------------------------------------------------

!nmmb   ELSE
!nmmb     CALL wrf_error_fatal('Lacking arguments for LSM in surface driver')
        ENDIF

     CASE (RUCLSMSCHEME)
!nmmb   IF (PRESENT(qv_curr)    .AND.  PRESENT(qc_curr) .AND.    &
!           PRESENT(emiss)      .AND.  PRESENT(t2)      .AND.    &
!nmmb       PRESENT(qsg)        .AND.  PRESENT(qvg)     .AND.    &
!nmmb       PRESENT(qcg)        .AND.  PRESENT(soilt1)  .AND.    &
!nmmb       PRESENT(tsnav)      .AND.  PRESENT(smfr3d)  .AND.    &
!nmmb       PRESENT(keepfr3dflag) .AND. PRESENT(rainbl) .AND.    &
!nmmb                                                  .TRUE. ) THEN

!nmmb       IF( PRESENT(sr) ) THEN
!nmmb           frpcpn=.true.
!nmmb       ELSE
!nmmb           SR = 1.
!nmmb       ENDIF

!nmmb       CALL wrf_debug(100,'in RUC LSM')
!nmmb       CALL LSMRUC(dtbl,itimestep,num_soil_layers,          &
!nmmb            zs,rainbl,snow,snowh,snowc,sr,frpcpn,           &
!nmmb            dz8w,p8w,t_phy,qv_curr,qc_curr,rho,             & !p8w in [pa]
!nmmb            glw,gsw,emiss,chklowq,                          &
!nmmb            flqc,flhc,mavail,canwat,vegfra,albedo,znt,      &
!nmmb            qsfc,qsg,qvg,qcg,soilt1,tsnav,                  &
!nmmb            tmn,ivgtyp,isltyp,xland,xice,                   &
!nmmb            cp,g,xlv,stbolt,                                &
!nmmb            smois,smstav,smstot,tslb,tsk,hfx,qfx,lh,        &
!nmmb            sfcrunoff,udrunoff,sfcexc,                      &
!nmmb            sfcevp,grdflx,acsnow,                           &
!nmmb            smfr3d,keepfr3dflag,                            &
!nmmb            myj,                                            &
!nmmb            ids,ide, jds,jde, kds,kde,                      &
!nmmb            ims,ime, jms,jme, kms,kme,                      &
!nmmb            i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )

!nmmb     IF(.not. MYJ) then

!nmmb      CALL SFCDIAGS(HFX,QFX,TSK,QVG,CHS2,CQS2,T2,TH2,Q2,      &
!nmmb                 PSFC,CP,R_d,RCP,                              &
!nmmb                 ids,ide, jds,jde, kds,kde,                    &
!nmmb                 ims,ime, jms,jme, kms,kme,                    &
!nmmb         i_start(ij),i_end(ij), j_start(ij),j_end(ij), kts,kte    )
!nmmb     ENDIF
 

!nmmb   ELSE
!nmmb     CALL wrf_error_fatal('Lacking arguments for RUCLSM in surface driver')
!nmmb   ENDIF

     CASE DEFAULT

!nmmb   WRITE( message , * ) &
       WRITE( 0 , * ) &
        'The surface option does not exist: sf_surface_physics = ', sf_surface_physics
!nmmb   CALL wrf_error_fatal ( message )

     END SELECT sfc_select
     ENDDO
!jaa      !$omp end parallel do

 430 CONTINUE


! Reset RAINBL in mm (Accumulation between PBL calls)

     IF ( PRESENT( rainbl ) ) THEN
!jaa        !$omp parallel do   &
!jaa        !$omp private ( ij, i, j, k )
       DO ij = 1 , num_tiles
         DO j=j_start(ij),j_end(ij)
         DO i=i_start(ij),i_end(ij)
            RAINBL(i,j) = 0.
         ENDDO
         ENDDO
       ENDDO
!jaa        !$omp end parallel do
     ENDIF

   ENDIF

   END SUBROUTINE surface_driver
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE MYJSFC(ITIMESTEP,HT,DZ                                & 
     &                 ,PMID,PINT,TH,T,QV,QC,U,V,Q2                    &
     &                 ,TSK,QSFC,THZ0,QZ0,UZ0,VZ0                      &
     &                 ,LOWLYR,XLAND                                   &
     &                 ,USTAR,ZNT,Z0BASE,PBLH,MAVAIL,RMOL              &
     &                 ,AKHS,AKMS                                      &
     &                 ,CHS,CHS2,CQS2,HFX,QFX,FLX_LH,FLHC,FLQC         &
     &                 ,QGH,CPM,CT                                     &
     &                 ,U10,V10,T02,TH02,TSHLTR,TH10,Q02,QSHLTR,Q10    &
     &                 ,PSHLTR,RIB                                     & ! Added Bulk Richardson No.
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                 ,IMS,IME,JMS,JME,KMS,KME                        &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,INTENT(IN) :: ITIMESTEP
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: LOWLYR
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: HT,MAVAIL,TSK      &
     &                                             ,XLAND,Z0BASE
!
      REAL,DIMENSION(IMS:IME,KMS:KME,JMS:JME),INTENT(IN) :: DZ         &
     &                                                     ,PMID,PINT  &
     &                                                     ,Q2,QC,QV   &
     &                                                     ,T,TH       &
     &                                                     ,U,V   
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: FLX_LH,HFX,PSHLTR &
     &                                              ,QFX,Q10,QSHLTR    &
     &                                              ,TH10,TSHLTR,T02       &
     &                                              ,U10,V10,TH02,Q02
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: AKHS,AKMS       &
     &                                                ,PBLH,QSFC,RIB
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: QZ0,RMOL,THZ0   &
     &                                                ,USTAR,UZ0,VZ0   &
     &                                                ,ZNT
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(OUT) :: CHS,CHS2,CQS2     &
     &                                              ,CPM,CT,FLHC,FLQC  &
     &                                              ,QGH 
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: I,J,K,KFLIP,LMH,LPBL,NTSD
!
      REAL :: A,APESFC,B,BTGX,CWMLOW                                   &
     &       ,DQDT,DTDIF,DTDT,DUDT,DVDT                                &
     &       ,FIS                                                      &
     &       ,P02P,P10P,PLOW,PSFC,PTOP,QLOW,QS02,QS10                  &
     &       ,RAPA,RAPA02,RAPA10,RATIOMX,RDZ,SEAMASK,SM                &
     &       ,T02P,T10P,TEM,TH02P,TH10P,THLOW,THELOW,THM               &
     &       ,TLOW,TZ0,ULOW,VLOW,ZSL
!
      REAL,DIMENSION(KTS:KTE) :: CWMK,PK,Q2K,QK,THEK,THK,TK,UK,VK
!
      REAL,DIMENSION(KTS:KTE-1) :: EL,ELM
!
      REAL,DIMENSION(KTS:KTE+1) :: ZHK
!
      REAL,DIMENSION(ITS:ITE,JTS:JTE) :: THSK
!
      REAL,DIMENSION(ITS:ITE,KTS:KTE+1,JTS:JTE) :: ZINT
!
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
!
!***  MAKE PREPARATIONS
!
!----------------------------------------------------------------------
      DO J=JTS,JTE
      DO K=KTS,KTE+1
      DO I=ITS,ITE
        ZINT(I,K,J)=0.
      ENDDO
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO I=ITS,ITE
        ZINT(I,KTE+1,J)=HT(I,J)     ! Z at bottom of lowest sigma layer
        PBLH(I,J)=-1.
!
!!!!!!!!!
!!!!!! UNCOMMENT THESE LINES IF USING ETA COORDINATES
!!!!!!!!!
!!!!!!  ZINT(I,KTE+1,J)=1.E-4       ! Z of bottom of lowest eta layer
!!!!!!  ZHK(KTE+1)=1.E-4            ! Z of bottom of lowest eta layer
!
      ENDDO
      ENDDO
!
      DO J=JTS,JTE
      DO K=KTE,KTS,-1
        KFLIP=KTE+1-K
        DO I=ITS,ITE
          ZINT(I,K,J)=ZINT(I,K+1,J)+DZ(I,KFLIP,J)
        ENDDO
      ENDDO
      ENDDO
!
      NTSD=ITIMESTEP
!
!#if ( NMM_CORE == 1 )
      IF(NTSD==0) then
!#else
!     IF(NTSD==1)THEN
!#endif
!tgs       IF(NTSD==1)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          USTAR(I,J)=0.1
          FIS=HT(I,J)*G
          SM=XLAND(I,J)-1.
!!!       Z0(I,J)=SM*Z0SEA+(1.-SM)*(Z0(I,J)*Z0MAX+FIS*FCM+Z0LAND)
!!!       ZNT(I,J)=SM*Z0SEA+(1.-SM)*(ZNT(I,J)*Z0MAX+FIS*FCM+Z0LAND)
        ENDDO
        ENDDO
      ENDIF
!
!!!!  IF(NTSD==1)THEN
        DO J=JTS,JTE
        DO I=ITS,ITE
          CT(I,J)=0.
        ENDDO
        ENDDO
!!!!  ENDIF
!
!......................................................................
!$omp parallel do &
!$omp private (j,i,lmh,ptop,psfc,seamask,k,kflip,thk,tk,ratiomx,qk,pk, &
!$omp          cwmk,thek,q2k,zhk,uk,vk,lpbl,plow,tlow,thlow,thelow,    &
!$omp          qlow,cwmlow,ulow,vlow,zsl,apesfc,tz0,rapa,th02p,th10p,  &
!$omp          rapa02,rapa10,t02p,t10p,p02p,p10p,qs02,qs10)
!......................................................................
!----------------------------------------------------------------------
        setup_integration:  DO J=JTS,JTE
!----------------------------------------------------------------------
!
        DO I=ITS,ITE
!
!***  LOWEST LAYER ABOVE GROUND MUST BE FLIPPED
!
          LMH=KTE-LOWLYR(I,J)+1
!
          PTOP=PINT(I,KTE+1,J)      ! KTE+1=KME
          PSFC=PINT(I,LOWLYR(I,J),J)
! Define THSK here (for first timestep mostly)
          THSK(I,J)=TSK(I,J)/(PSFC*1.E-5)**CAPA
!
!***  CONVERT LAND MASK (1 FOR SEA; 0 FOR LAND)
!
          SEAMASK=XLAND(I,J)-1.
!
!***  FILL 1-D VERTICAL ARRAYS
!***  AND FLIP DIRECTION SINCE MYJ SCHEME
!***  COUNTS DOWNWARD FROM THE DOMAIN'S TOP
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            THK(K)=TH(I,KFLIP,J)
            TK(K)=T(I,KFLIP,J)
            RATIOMX=QV(I,KFLIP,J)
            QK(K)=RATIOMX/(1.+RATIOMX)
            PK(K)=PMID(I,KFLIP,J)
            CWMK(K)=QC(I,KFLIP,J)
            THEK(K)=(CWMK(K)*(-ELOCP/TK(K))+1.)*THK(K)
            Q2K(K)=2.*Q2(I,KFLIP,J)
!
!
!***  COMPUTE THE HEIGHTS OF THE LAYER INTERFACES
!
            ZHK(K)=ZINT(I,K,J)
!
          ENDDO
          ZHK(KTE+1)=HT(I,J)          ! Z at bottom of lowest sigma layer
!
          DO K=KTE,KTS,-1
            KFLIP=KTE+1-K
            UK(K)=U(I,KFLIP,J)
            VK(K)=V(I,KFLIP,J)
          ENDDO
!
!***  FIND THE HEIGHT OF THE PBL
!
          LPBL=LMH
          DO K=LMH-1,1,-1
            IF(Q2K(K)<=EPSQ2*FH) THEN
              LPBL=K
              GO TO 110
            ENDIF
          ENDDO
!
          LPBL=1
!
!-----------------------------------------------------------------------
!--------------THE HEIGHT OF THE PBL------------------------------------
!-----------------------------------------------------------------------
!
 110      PBLH(I,J)=ZHK(LPBL)-ZHK(LMH+1)
!
!----------------------------------------------------------------------
!***
!***  FIND THE SURFACE EXCHANGE COEFFICIENTS
!***
!----------------------------------------------------------------------
          PLOW=PK(LMH)
          TLOW=TK(LMH)
          THLOW=THK(LMH)
          THELOW=THEK(LMH)
          QLOW=QK(LMH)
          CWMLOW=CWMK(LMH)
          ULOW=UK(LMH)
          VLOW=VK(LMH)
          ZSL=(ZHK(LMH)-ZHK(LMH+1))*0.5
          APESFC=(PSFC*1.E-5)**CAPA
!tgs - in ARW THZ0 is not initialized when MYJSFC is called first time
!#if ( NMM_CORE == 1 )
      if(NTSD==0) then
!#else
!     IF(NTSD==1)THEN
!#endif
!       if(itimestep.le.1) then
          TZ0=TSK(I,J)
       else
          TZ0=THZ0(I,J)*APESFC
       endif
!
          CALL SFCDIF(NTSD,SEAMASK,THSK(I,J),QSFC(I,J),PSFC            &
     &               ,UZ0(I,J),VZ0(I,J),TZ0,THZ0(I,J),QZ0(I,J)         &
     &               ,USTAR(I,J),ZNT(I,J),Z0BASE(I,J),CT(I,J),RMOL(I,J) &
     &               ,AKMS(I,J),AKHS(I,J),PBLH(I,J),MAVAIL(I,J)        &
     &               ,CHS(I,J),CHS2(I,J),CQS2(I,J)                     &
     &               ,HFX(I,J),QFX(I,J),FLX_LH(I,J)                    &
     &               ,FLHC(I,J),FLQC(I,J),QGH(I,J),CPM(I,J)            &
     &               ,ULOW,VLOW,TLOW,THLOW,THELOW,QLOW,CWMLOW          &
     &               ,ZSL,PLOW                                         &
     &               ,U10(I,J),V10(I,J),TSHLTR(I,J),TH10(I,J)          &
     &               ,QSHLTR(I,J),Q10(I,J),PSHLTR(I,J)                 &
     &               ,IDS,IDE,JDS,JDE,KDS,KDE                          &
     &               ,IMS,IME,JMS,JME,KMS,KME                          &
     &               ,ITS,ITE,JTS,JTE,KTS,KTE,i,j,ZHK(LMH+1),RIB(I,J))   ! Added Bulk Richardson No.
!
!***  REMOVE SUPERATURATION AT 2M AND 10M
!
          RAPA=APESFC
          TH02P=TSHLTR(I,J)
          TH10P=TH10(I,J)
          TH02(I,J)=TSHLTR(I,J)
!
          RAPA02=RAPA-GOCP02/TH02P
          RAPA10=RAPA-GOCP10/TH10P
!
          T02P=TH02P*RAPA02
          T10P=TH10P*RAPA10
! 1 may 06 tgs         T02(I,J) = T02P
          T02(I,J) = TH02(I,J)*APESFC
!
          P02P=(RAPA02**RCAP)*1.E5
          P10P=(RAPA10**RCAP)*1.E5
!
          QS02=PQ0/P02P*EXP(A2*(T02P-A3)/(T02P-A4))
          QS10=PQ0/P10P*EXP(A2*(T10P-A3)/(T10P-A4))
!
          IF(QSHLTR(I,J)>QS02)QSHLTR(I,J)=QS02
          IF(Q10   (I,J)>QS10)Q10   (I,J)=QS10
          Q02(I,J)=QSHLTR(I,J)/(1.-QSHLTR(I,J))
!----------------------------------------------------------------------
!
        ENDDO
!
!----------------------------------------------------------------------
!
      ENDDO setup_integration
!
!......................................................................
!$omp end parallel do
!......................................................................
!----------------------------------------------------------------------
 
      END SUBROUTINE MYJSFC
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!----------------------------------------------------------------------
      SUBROUTINE SFCDIF(NTSD,SEAMASK,THS,QS,PSFC                       &
     &                 ,UZ0,VZ0,TZ0,THZ0,QZ0                           &
     &                 ,USTAR,Z0,Z0BASE,CT,RLMO,AKMS,AKHS,PBLH,WETM    &
     &                 ,CHS,CHS2,CQS2,HFX,QFX,FLX_LH,FLHC,FLQC,QGH,CPM &
     &                 ,ULOW,VLOW,TLOW,THLOW,THELOW,QLOW,CWMLOW        &
     &                 ,ZSL,PLOW                                       &
     &                 ,U10,V10,TH02,TH10,Q02,Q10,PSHLTR               &
     &                 ,IDS,IDE,JDS,JDE,KDS,KDE                        &
     &                 ,IMS,IME,JMS,JME,KMS,KME                        &
     &                 ,ITS,ITE,JTS,JTE,KTS,KTE,i,j,ZSFC,RIB)            ! Added Bulk Richardson No.
!     ****************************************************************
!     *                                                              *
!     *                       SURFACE LAYER                          *
!     *                                                              *
!     ****************************************************************
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE,i,j
!
      INTEGER,INTENT(IN) :: NTSD
!
      REAL,INTENT(IN) :: CWMLOW,PBLH,PLOW,QLOW,PSFC,SEAMASK,ZSFC       &
     &                  ,THELOW,THLOW,THS,TLOW,TZ0,ULOW,VLOW,WETM,ZSL  &
     &                  ,Z0BASE
!
      REAL,INTENT(OUT) :: CHS,CHS2,CPM,CQS2,CT,FLHC,FLQC,FLX_LH,HFX    &
     &              ,RIB,PSHLTR,Q02,Q10,QFX,QGH,RLMO,TH02,TH10,U10,V10
!
      REAL,INTENT(INOUT) :: AKHS,AKMS,QZ0,THZ0,USTAR,UZ0,VZ0,Z0,QS
!----------------------------------------------------------------------
!***
!***  LOCAL VARIABLES
!***
      INTEGER :: ITR,K
!
      REAL :: A,B,BTGH,BTGX,CXCHL,CXCHS,DTHV,DU2,ELFC,FCT              &
     &       ,HLFLX,HSFLX,HV,PSH02,PSH10,PSHZ,PSHZL,PSM10,PSMZ,PSMZL   &
     &       ,RDZ,RDZT,RLMA,RLMN,RLMP                                  &
     &       ,RLOGT,RLOGU,RWGH,RZ,RZST,RZSU,SIMH,SIMM,TEM,THM          &
     &       ,UMFLX,USTARK,VMFLX,WGHT,WGHTT,WGHTQ,WSTAR2               &
     &       ,X,XLT,XLT4,XLU,XLU4,XT,XT4,XU,XU4,ZETALT,ZETALU          &
     &       ,ZETAT,ZETAU,ZQ,ZSLT,ZSLU,ZT,ZU,TOPOTERM,ZZIL
!
!***  DIAGNOSTICS
!
      REAL :: AKHS02,AKHS10,AKMS02,AKMS10,EKMS10,QSAT10,QSAT2          &
     &       ,RLNT02,RLNT10,RLNU10,SIMH02,SIMH10,SIMM10,T02,T10        &
     &       ,TERM1,RLOW,U10E,V10E,WSTAR,XLT02,XLT024,XLT10            &
     &       ,XLT104,XLU10,XLU104,XU10,XU104,ZT02,ZT10,ZTAT02,ZTAT10   &
     &       ,ZTAU,ZTAU10,ZU10,ZUUZ
!----------------------------------------------------------------------
!**********************************************************************
!----------------------------------------------------------------------
      RDZ=1./ZSL
      CXCHL=EXCML*RDZ
      CXCHS=EXCMS*RDZ
!
      BTGX=G/THLOW
      ELFC=VKARMAN*BTGX
!
      IF(PBLH>1000.)THEN
        BTGH=BTGX*PBLH
      ELSE
        BTGH=BTGX*1000.
      ENDIF
!
!----------------------------------------------------------------------
!
!***  SEA POINTS
!
!----------------------------------------------------------------------
!
      IF(SEAMASK>0.5)THEN 
!
!----------------------------------------------------------------------
        DO ITR=1,ITRMX
!----------------------------------------------------------------------
          Z0=MAX(USTFC*USTAR*USTAR,1.59E-5)
!
!***  VISCOUS SUBLAYER, JANJIC MWR 1994
!
!----------------------------------------------------------------------
          IF(USTAR<USTC)THEN
!----------------------------------------------------------------------
!
            IF(USTAR<USTR)THEN
!
!#if ( NMM_CORE == 1 )
      IF(NTSD==0) then
!#else
!     IF(NTSD==1)THEN
!#endif
!tgs              IF(NTSD==1)THEN
                AKMS=CXCHS
                AKHS=CXCHS
                QS=QLOW
              ENDIF
!
              ZU=FZU1*SQRT(SQRT(Z0*USTAR*RVISC))/USTAR
              WGHT=AKMS*ZU*RVISC
              RWGH=WGHT/(WGHT+1.)
              UZ0=(ULOW*RWGH+UZ0)*0.5
              VZ0=(VLOW*RWGH+VZ0)*0.5
!
              ZT=FZT1*ZU
              ZQ=FZQ1*ZT
              WGHTT=AKHS*ZT*RTVISC
              WGHTQ=AKHS*ZQ*RQVISC
!
!if ( NMM_CORE == 1 )
      IF(NTSD>0)THEN
!#else
!     IF(NTSD>1)THEN
!#endif
!tgs              IF(NTSD>1)THEN
                THZ0=((WGHTT*THLOW+THS)/(WGHTT+1.)+THZ0)*0.5
                QZ0=((WGHTQ*QLOW+QS)/(WGHTQ+1.)+QZ0)*0.5
              ELSE
                THZ0=(WGHTT*THLOW+THS)/(WGHTT+1.)
                QZ0=(WGHTQ*QLOW+QS)/(WGHTQ+1.)
              ENDIF
!
            ENDIF
!
            IF(USTAR>=USTR.AND.USTAR<USTC)THEN
              ZU=Z0
              UZ0=0.
              VZ0=0.
!
              ZT=FZT2*SQRT(SQRT(Z0*USTAR*RVISC))/USTAR
              ZQ=FZQ2*ZT
              WGHTT=AKHS*ZT*RTVISC
              WGHTQ=AKHS*ZQ*RQVISC
!
!if ( NMM_CORE == 1 )
      IF(NTSD>0)THEN
!#else
!     IF(NTSD>1)THEN
!#endif

!tgs              IF(NTSD>1)THEN
                THZ0=((WGHTT*THLOW+THS)/(WGHTT+1.)+THZ0)*0.5
                QZ0=((WGHTQ*QLOW+QS)/(WGHTQ+1.)+QZ0)*0.5
              ELSE
                THZ0=(WGHTT*THLOW+THS)/(WGHTT+1.)
                QZ0=(WGHTQ*QLOW+QS)/(WGHTQ+1.)
              ENDIF
!
            ENDIF
!----------------------------------------------------------------------
          ELSE
!----------------------------------------------------------------------
            ZU=Z0
            UZ0=0.
            VZ0=0.
!
            ZT=Z0
            THZ0=THS
!
            ZQ=Z0
            QZ0=QS
!----------------------------------------------------------------------
          ENDIF
!----------------------------------------------------------------------
          TEM=(TLOW+TZ0)*0.5
          THM=(THELOW+THZ0)*0.5
!
          A=THM*P608
          B=(ELOCP/TEM-1.-P608)*THM
!
          DTHV=((THELOW-THZ0)*((QLOW+QZ0+CWMLOW)*(0.5*P608)+1.)        &
     &        +(QLOW-QZ0+CWMLOW)*A+CWMLOW*B) 
!
          DU2=MAX((ULOW-UZ0)**2+(VLOW-VZ0)**2,EPSU2)
          RIB=BTGX*DTHV*ZSL/DU2
!----------------------------------------------------------------------
!         IF(RIB>=RIC)THEN
!----------------------------------------------------------------------
!           AKMS=MAX( VISC*RDZ,CXCHS)
!           AKHS=MAX(TVISC*RDZ,CXCHS)
!----------------------------------------------------------------------
!         ELSE  !  turbulent branch
!----------------------------------------------------------------------
            ZSLU=ZSL+ZU
            ZSLT=ZSL+ZT
!
            RZSU=ZSLU/ZU
            RZST=ZSLT/ZT
!
            RLOGU=LOG(RZSU)
            RLOGT=LOG(RZST)
!
!----------------------------------------------------------------------
!***  1./MONIN-OBUKHOV LENGTH
!----------------------------------------------------------------------
!
            RLMO=ELFC*AKHS*DTHV/USTAR**3
!
            ZETALU=ZSLU*RLMO
            ZETALT=ZSLT*RLMO
            ZETAU=ZU*RLMO
            ZETAT=ZT*RLMO
!
            ZETALU=MIN(MAX(ZETALU,ZTMIN1),ZTMAX1)
            ZETALT=MIN(MAX(ZETALT,ZTMIN1),ZTMAX1)
            ZETAU=MIN(MAX(ZETAU,ZTMIN1/RZSU),ZTMAX1/RZSU)
            ZETAT=MIN(MAX(ZETAT,ZTMIN1/RZST),ZTMAX1/RZST)
!
!----------------------------------------------------------------------
!***   WATER FUNCTIONS
!----------------------------------------------------------------------
!
            RZ=(ZETAU-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZ=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
!
            RZ=(ZETALU-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZL=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
!
            SIMM=PSMZL-PSMZ+RLOGU
!
            RZ=(ZETAT-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZ=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
            RZ=(ZETALT-ZTMIN1)/DZETA1
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZL=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
            SIMH=(PSHZL-PSHZ+RLOGT)*FH01
!----------------------------------------------------------------------
            USTARK=USTAR*VKARMAN
      if(abs(simm)<1.e-10.or.abs(simh)<1.e-10)then
        write(0,*)' simm=',simm,' simh=',simh,' at i=',i,' j=',j
      endif
            AKMS=MAX(USTARK/SIMM,CXCHS)
            AKHS=MAX(USTARK/SIMH,CXCHS)
!
!----------------------------------------------------------------------
!***  BELJAARS CORRECTION FOR USTAR
!----------------------------------------------------------------------
!
            IF(DTHV<=0.)THEN                                           !zj
              WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)                !zj
            ELSE                                                       !zj
              WSTAR2=0.                                                !zj
            ENDIF                                                      !zj
            USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
!
!----------------------------------------------------------------------
!         ENDIF  !  End of turbulent branch
!----------------------------------------------------------------------
!
        ENDDO  !  End of the iteration loop over sea points
!
!----------------------------------------------------------------------
!
!***  LAND POINTS
!
!----------------------------------------------------------------------
!
      ELSE  
!
!----------------------------------------------------------------------
!nmmb
!!!     IF(NTSD==1)THEN
        IF(NTSD==0)THEN
!nmmb
          QS=QLOW
        ENDIF
!
        ZU=Z0
        UZ0=0.
        VZ0=0.
!
        ZT=ZU*ZTFC
        THZ0=THS
!
        ZQ=ZT
        QZ0=QS
!----------------------------------------------------------------------
        TEM=(TLOW+TZ0)*0.5
        THM=(THELOW+THZ0)*0.5
!
        A=THM*P608
        B=(ELOCP/TEM-1.-P608)*THM
!
        DTHV=((THELOW-THZ0)*((QLOW+QZ0+CWMLOW)*(0.5*P608)+1.)          &
     &        +(QLOW-QZ0+CWMLOW)*A+CWMLOW*B) 
!
        DU2=MAX((ULOW-UZ0)**2+(VLOW-VZ0)**2,EPSU2)
        RIB=BTGX*DTHV*ZSL/DU2
!----------------------------------------------------------------------
!       IF(RIB>=RIC)THEN
!         AKMS=MAX( VISC*RDZ,CXCHL)
!         AKHS=MAX(TVISC*RDZ,CXCHL)
!----------------------------------------------------------------------
!       ELSE  !  Turbulent branch
!----------------------------------------------------------------------
          ZSLU=ZSL+ZU
!
          RZSU=ZSLU/ZU
!
          RLOGU=LOG(RZSU)

          ZSLT=ZSL+ZU ! u,v and t are at the same level
!----------------------------------------------------------------------
!
!mp   Topo modification of ZILFC term
!	
          TOPOTERM=TOPOFAC*ZSFC**2.
          TOPOTERM=MAX(TOPOTERM,3.0)
!
          IF(DTHV>0.)THEN
            ZZIL=ZILFC*TOPOTERM
          ELSE
            ZZIL=ZILFC
          ENDIF
!
!----------------------------------------------------------------------
!
          land_point_iteration: DO ITR=1,ITRMX
!
!----------------------------------------------------------------------
!***  ZILITINKEVITCH FIX FOR ZT
!----------------------------------------------------------------------
!
! oldform   ZT=MAX(EXP(ZZIL*SQRT(USTAR*ZU))*ZU,EPSZT)
            ZT=MAX(EXP(ZZIL*SQRT(USTAR*Z0BASE))*Z0BASE,EPSZT)
!
            RZST=ZSLT/ZT
            RLOGT=LOG(RZST)
!
!----------------------------------------------------------------------
!***  1./MONIN-OBUKHOV LENGTH-SCALE
!----------------------------------------------------------------------
!
            RLMO=ELFC*AKHS*DTHV/USTAR**3
            ZETALU=ZSLU*RLMO
            ZETALT=ZSLT*RLMO
            ZETAU=ZU*RLMO
            ZETAT=ZT*RLMO
!
            ZETALU=MIN(MAX(ZETALU,ZTMIN2),ZTMAX2)
            ZETALT=MIN(MAX(ZETALT,ZTMIN2),ZTMAX2)
            ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)
            ZETAT=MIN(MAX(ZETAT,ZTMIN2/RZST),ZTMAX2/RZST)
!
!----------------------------------------------------------------------
!***  LAND FUNCTIONS
!----------------------------------------------------------------------
!
            RZ=(ZETAU-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZ=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
            RZ=(ZETALU-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSMZL=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
            SIMM=PSMZL-PSMZ+RLOGU
!
            RZ=(ZETAT-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZ=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
            RZ=(ZETALT-ZTMIN2)/DZETA2
            K=INT(RZ)
            RDZT=RZ-REAL(K)
            K=MIN(K,KZTM2)
            K=MAX(K,0)
            PSHZL=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
            SIMH=(PSHZL-PSHZ+RLOGT)*FH02
!----------------------------------------------------------------------
            USTARK=USTAR*VKARMAN
            AKMS=MAX(USTARK/SIMM,CXCHL)
            AKHS=MAX(USTARK/SIMH,CXCHL)
!
!----------------------------------------------------------------------
!***  BELJAARS CORRECTION FOR USTAR
!----------------------------------------------------------------------
!
            IF(DTHV<=0.)THEN                                           !zj
              WSTAR2=WWST2*ABS(BTGH*AKHS*DTHV)**(2./3.)                !zj
            ELSE                                                       !zj
              WSTAR2=0.                                                !zj
            ENDIF                                                      !zj
!
            USTAR=MAX(SQRT(AKMS*SQRT(DU2+WSTAR2)),EPSUST)
!
!----------------------------------------------------------------------
          ENDDO land_point_iteration
!----------------------------------------------------------------------
!
!       ENDIF  !  End of turbulant branch over land
!
!----------------------------------------------------------------------
!
      ENDIF  !  End of land/sea branch
!
!----------------------------------------------------------------------
!***  COUNTERGRADIENT FIX
!----------------------------------------------------------------------
!
!     HV=-AKHS*DTHV
!     IF(HV>0.)THEN
!       FCT=-10.*(BTGX)**(-1./3.)
!       CT=FCT*(HV/(PBLH*PBLH))**(2./3.)
!     ELSE
       CT=0.
!     ENDIF
!
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!***  THE FOLLOWING DIAGNOSTIC BLOCK PRODUCES 2-m and 10-m VALUES
!***  FOR TEMPERATURE, MOISTURE, AND WINDS.  IT IS DONE HERE SINCE
!***  THE VARIOUS QUANTITIES NEEDED FOR THE COMPUTATION ARE LOST
!***  UPON EXIT FROM THE ROTUINE.
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!
      WSTAR=SQRT(WSTAR2)/WWST
!
      UMFLX=AKMS*(ULOW -UZ0 )
      VMFLX=AKMS*(VLOW -VZ0 )
      HSFLX=AKHS*(THLOW-THZ0)
      HLFLX=AKHS*(QLOW -QZ0 )
!----------------------------------------------------------------------
!     IF(RIB>=RIC)THEN
!----------------------------------------------------------------------
!       IF(SEAMASK>0.5)THEN
!         AKMS10=MAX( VISC/10.,CXCHS)
!         AKHS02=MAX(TVISC/02.,CXCHS)
!         AKHS10=MAX(TVISC/10.,CXCHS)
!       ELSE
!         AKMS10=MAX( VISC/10.,CXCHL)
!         AKHS02=MAX(TVISC/02.,CXCHL)
!         AKHS10=MAX(TVISC/10.,CXCHL)
!       ENDIF
!----------------------------------------------------------------------
!     ELSE
!----------------------------------------------------------------------
        ZU10=ZU+10.
        ZT02=ZT+02.
        ZT10=ZT+10.
!
        RLNU10=LOG(ZU10/ZU)
        RLNT02=LOG(ZT02/ZT)
        RLNT10=LOG(ZT10/ZT)
!
        ZTAU10=ZU10*RLMO
        ZTAT02=ZT02*RLMO
        ZTAT10=ZT10*RLMO
!
!----------------------------------------------------------------------
!***  SEA
!----------------------------------------------------------------------
!
        IF(SEAMASK>0.5)THEN
!
!----------------------------------------------------------------------
          ZTAU10=MIN(MAX(ZTAU10,ZTMIN1),ZTMAX1)
          ZTAT02=MIN(MAX(ZTAT02,ZTMIN1),ZTMAX1)
          ZTAT10=MIN(MAX(ZTAT10,ZTMIN1),ZTMAX1)
!----------------------------------------------------------------------
          RZ=(ZTAU10-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSM10=(PSIM1(K+2)-PSIM1(K+1))*RDZT+PSIM1(K+1)
!
          SIMM10=PSM10-PSMZ+RLNU10
!
          RZ=(ZTAT02-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH02=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
          SIMH02=(PSH02-PSHZ+RLNT02)*FH01
!
          RZ=(ZTAT10-ZTMIN1)/DZETA1
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH10=(PSIH1(K+2)-PSIH1(K+1))*RDZT+PSIH1(K+1)
!
          SIMH10=(PSH10-PSHZ+RLNT10)*FH01
!
          AKMS10=MAX(USTARK/SIMM10,CXCHS)
          AKHS02=MAX(USTARK/SIMH02,CXCHS)
          AKHS10=MAX(USTARK/SIMH10,CXCHS)
!
!----------------------------------------------------------------------
!***  LAND
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
          ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
          ZTAT02=MIN(MAX(ZTAT02,ZTMIN2),ZTMAX2)
          ZTAT10=MIN(MAX(ZTAT10,ZTMIN2),ZTMAX2)
!----------------------------------------------------------------------
          RZ=(ZTAU10-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
!
          SIMM10=PSM10-PSMZ+RLNU10
!
          RZ=(ZTAT02-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH02=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
          SIMH02=(PSH02-PSHZ+RLNT02)*FH02
!
          RZ=(ZTAT10-ZTMIN2)/DZETA2
          K=INT(RZ)
          RDZT=RZ-REAL(K)
          K=MIN(K,KZTM2)
          K=MAX(K,0)
          PSH10=(PSIH2(K+2)-PSIH2(K+1))*RDZT+PSIH2(K+1)
!
          SIMH10=(PSH10-PSHZ+RLNT10)*FH02
!
          AKMS10=USTARK/SIMM10
          AKHS02=USTARK/SIMH02
          AKHS10=USTARK/SIMH10
!
          IF(AKMS10<=CXCHL) AKMS10=AKMS
          IF(AKHS02<=CXCHL) AKHS02=AKHS
          IF(AKHS10<=CXCHL) AKHS10=AKHS
!
!----------------------------------------------------------------------
        ENDIF
!----------------------------------------------------------------------
!     ENDIF
!----------------------------------------------------------------------
!
      U10 =UMFLX/AKMS10+UZ0
      V10 =VMFLX/AKMS10+VZ0
      TH02=HSFLX/AKHS02+THZ0
      TH10=HSFLX/AKHS10+THZ0
      Q02 =HLFLX/AKHS02+QZ0
      Q10 =HLFLX/AKHS10+QZ0
      TERM1=-0.068283/TLOW
      PSHLTR=PSFC*EXP(TERM1)
!
!----------------------------------------------------------------------
!***  COMPUTE "EQUIVALENT" Z0 TO APPROXIMATE LOCAL SHELTER READINGS.
!----------------------------------------------------------------------
!
      U10E=U10
      V10E=V10
!
      IF(SEAMASK<0.5)THEN

!1st        ZUUZ=MIN(0.5*ZU,0.1)
!1st        ZU=MAX(0.1*ZU,ZUUZ)
!tst        ZUUZ=amin1(ZU*0.50,0.3)
!tst        ZU=amax1(ZU*0.3,ZUUZ)

        ZUUZ=AMIN1(ZU*0.50,0.18)
        ZU=AMAX1(ZU*0.35,ZUUZ)
!
        ZU10=ZU+10.
        RZSU=ZU10/ZU
        RLNU10=LOG(RZSU)

        ZETAU=ZU*RLMO
        ZTAU10=ZU10*RLMO

        ZTAU10=MIN(MAX(ZTAU10,ZTMIN2),ZTMAX2)
        ZETAU=MIN(MAX(ZETAU,ZTMIN2/RZSU),ZTMAX2/RZSU)

        RZ=(ZTAU10-ZTMIN2)/DZETA2
        K=INT(RZ)
        RDZT=RZ-REAL(K)
        K=MIN(K,KZTM2)
        K=MAX(K,0)
        PSM10=(PSIM2(K+2)-PSIM2(K+1))*RDZT+PSIM2(K+1)
        SIMM10=PSM10-PSMZ+RLNU10
        EKMS10=MAX(USTARK/SIMM10,CXCHL)

        U10E=UMFLX/EKMS10+UZ0
        V10E=VMFLX/EKMS10+VZ0

      ENDIF
!
      U10=U10E
      V10=V10E
!
!----------------------------------------------------------------------
!***  SET OTHER WRF DRIVER ARRAYS
!----------------------------------------------------------------------
!
      RLOW=PLOW/(R_D*TLOW)
      CHS=AKHS
      CHS2=AKHS02
      CQS2=AKHS02
      HFX=-RLOW*CP*HSFLX
      QFX=-RLOW*HLFLX*WETM
      FLX_LH=XLV*QFX
      FLHC=RLOW*CP*AKHS
      FLQC=RLOW*AKHS*WETM
!!!   QGH=PQ0/PSHLTR*EXP(A2S*(TSK-A3S)/(TSK-A4S))
      QGH=((1.-SEAMASK)*PQ0+SEAMASK*PQ0SEA)                            &
     &     /PLOW*EXP(A2S*(TLOW-A3S)/(TLOW-A4S))
      QGH=QGH/(1.-QGH)    !Convert to mixing ratio
      CPM=CP*(1.+0.8*QLOW)
!
!***  DO NOT COMPUTE QS OVER LAND POINTS HERE SINCE IT IS
!***  A PROGNOSTIC VARIABLE THERE.  IT IS OKAY TO USE IT 
!***  AS A DIAGNOSTIC OVER WATER SINCE IT WILL CAUSE NO
!***  INTERFERENCE BEFORE BEING RECOMPUTED IN MYJPBL.
!
      IF(SEAMASK>0.5)THEN
        QS=QLOW+QFX/(RLOW*AKHS)
        QS=QS/(1.-QS)
      ENDIF
!----------------------------------------------------------------------
!
      END SUBROUTINE SFCDIF
!
!----------------------------------------------------------------------
      SUBROUTINE MYJSFC_INIT(LOWLYR,USTAR,Z0                           &
     &                     ,SEAMASK,XICE,IVGTYP,RESTART                &
     &                     ,ALLOWED_TO_READ                            &
     &                     ,IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE)
!----------------------------------------------------------------------
      IMPLICIT NONE
!----------------------------------------------------------------------
      LOGICAL,INTENT(IN) :: RESTART,ALLOWED_TO_READ
!
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,KDS,KDE                    &
     &                     ,IMS,IME,JMS,JME,KMS,KME                    &
     &                     ,ITS,ITE,JTS,JTE,KTS,KTE
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: IVGTYP
!
      INTEGER,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: LOWLYR
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(IN) :: SEAMASK,XICE
!
      REAL,DIMENSION(IMS:IME,JMS:JME),INTENT(INOUT) :: USTAR,Z0
!
!----------------------------------------------------------------------
!***  LOCAL VARIABLES
!----------------------------------------------------------------------
!
      REAL,DIMENSION(0:30) :: VZ0TBL
      REAL,DIMENSION(0:30) :: VZ0TBL_24
!
      INTEGER :: I,IDUM,IRECV,J,JDUM,K,ITF,JTF,KTF,MAXGBL_IVGTYP       &
     &,          MAXLOC_IVGTYP
!
!     INTEGER :: MPI_INTEGER,MPI_MAX
!
      REAL :: SM_LOC,X,ZETA1,ZETA2,ZRNG1,ZRNG2
!
      REAL :: PIHF=3.1415926/2.,EPS=1.E-6
!----------------------------------------------------------------------
      VZ0TBL=                                                          &
     &  (/0.,                                                          &
     &    2.653,0.826,0.563,1.089,0.854,0.856,0.035,0.238,0.065,0.076  &
     &   ,0.011,0.035,0.011,0.000,0.000,0.000,0.000,0.000,0.000,0.000  &
     &   ,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000/)

      VZ0TBL_24= (/0.,                                                 &
     &            1.00,  0.07,  0.07,  0.07,  0.07,  0.15,             &
     &            0.08,  0.03,  0.05,  0.86,  0.80,  0.85,             &
     &            2.65,  1.09,  0.80,  0.001, 0.04,  0.05,             &
     &            0.01,  0.04,  0.06,  0.05,  0.03,  0.001,            &
     &            0.000, 0.000, 0.000, 0.000, 0.000, 0.000 /)

!----------------------------------------------------------------------
!
      JTF=MIN0(JTE,JDE-1)
      KTF=MIN0(KTE,KDE-1)
      ITF=MIN0(ITE,IDE-1)
!
!
!***  FOR NOW, ASSUME SIGMA MODE FOR LOWEST MODEL LAYER
!
      DO J=JTS,JTF
      DO I=ITS,ITF
        LOWLYR(I,J)=1
!       USTAR(I,J)=EPSUST
      ENDDO
      ENDDO
!----------------------------------------------------------------------
!nmmb#if (NMM_CORE == 1)
!
      IF(.NOT.RESTART)THEN
!       CALL WRF_GET_DM_COMMUNICATOR(MPI_COMM_COMP)
        MAXLOC_IVGTYP=MAXVAL(IVGTYP)
        CALL MPI_ALLREDUCE(MAXLOC_IVGTYP,MAXGBL_IVGTYP,1,MPI_INTEGER   &
     &,                    MPI_MAX,MPI_COMM_COMP,IRECV)
        MAXGBL_IVGTYP=MAXVAL(IVGTYP)
!       CALL WRF_DM_MAXVAL(MAXGBL_IVGTYP,IDUM,JDUM)
!
        IF (MAXGBL_IVGTYP<13) THEN
          DO J=JTS,JTE
          DO I=ITS,ITE
            SM_LOC=SEAMASK(I,J)-1.
            IF(SM_LOC+XICE(I,J)<0.5)THEN
!zj              Z0(I,J)=VZ0TBL(IVGTYP(I,J))
            ENDIF
          ENDDO
          ENDDO
!
        ELSE
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            SM_LOC=SEAMASK(I,J)-1.
            IF(SM_LOC+XICE(I,J)<0.5)THEN
!zj              Z0(I,J)=VZ0TBL_24(IVGTYP(I,J))
            ENDIF
          ENDDO
          ENDDO
!
        ENDIF ! Vegtype check
!
      ENDIF ! Restart check

!nmmb#endif
!----------------------------------------------------------------------
      IF(.NOT.RESTART)THEN
        DO J=JTS,JTE
        DO I=ITS,ITF
          USTAR(I,J)=0.1
        ENDDO
        ENDDO
      ENDIF

!----------------------------------------------------------------------
!
!***  COMPUTE SURFACE LAYER INTEGRAL FUNCTIONS
!
!----------------------------------------------------------------------
      FH01=1.
      FH02=1.
!
!      ZTMIN1=-10.0
!      ZTMAX1=2.0
!      ZTMIN2=-10.0
!      ZTMAX2=2.0
      ZTMIN1=-5.0
      ZTMAX1=1.0
      ZTMIN2=-5.0
      ZTMAX2=1.0
!
      ZRNG1=ZTMAX1-ZTMIN1
      ZRNG2=ZTMAX2-ZTMIN2
!
      DZETA1=ZRNG1/(KZTM-1)
      DZETA2=ZRNG2/(KZTM-1)
!
!----------------------------------------------------------------------
!***  FUNCTION DEFINITION LOOP
!----------------------------------------------------------------------
!
      ZETA1=ZTMIN1
      ZETA2=ZTMIN2
!
      DO K=1,KZTM
!
!----------------------------------------------------------------------
!***  UNSTABLE RANGE
!----------------------------------------------------------------------
!
        IF(ZETA1<0.)THEN
!
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
          X=SQRT(SQRT(1.-16.*ZETA1))
!
          PSIM1(K)=-2.*LOG((X+1.)/2.)-LOG((X*X+1.)/2.)+2.*ATAN(X)-PIHF
          PSIH1(K)=-2.*LOG((X*X+1.)/2.)
!
!----------------------------------------------------------------------
!***  STABLE RANGE
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
!
!         PSIM1(K)=5.*ZETA1
!         PSIH1(K)=5.*ZETA1
!----------------------------------------------------------------------
!***   HOLTSLAG AND DE BRUIN 1988
!----------------------------------------------------------------------
!
          PSIM1(K)=0.7*ZETA1+0.75*ZETA1*(6.-0.35*ZETA1)*EXP(-0.35*ZETA1)
          PSIH1(K)=0.7*ZETA1+0.75*ZETA1*(6.-0.35*ZETA1)*EXP(-0.35*ZETA1)
!----------------------------------------------------------------------
!
        ENDIF
!
!----------------------------------------------------------------------
!***  UNSTABLE RANGE
!----------------------------------------------------------------------
!
        IF(ZETA2<0.)THEN
!
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
!
          X=SQRT(SQRT(1.-16.*ZETA2))
!
          PSIM2(K)=-2.*LOG((X+1.)/2.)-LOG((X*X+1.)/2.)+2.*ATAN(X)-PIHF
          PSIH2(K)=-2.*LOG((X*X+1.)/2.)
!----------------------------------------------------------------------
!***  STABLE RANGE
!----------------------------------------------------------------------
!
        ELSE
!
!----------------------------------------------------------------------
!***  PAULSON 1970 FUNCTIONS
!----------------------------------------------------------------------
!
!         PSIM2(K)=5.*ZETA2
!         PSIH2(K)=5.*ZETA2
!
!----------------------------------------------------------------------
!***  HOLTSLAG AND DE BRUIN 1988
!----------------------------------------------------------------------
!
          PSIM2(K)=0.7*ZETA2+0.75*ZETA2*(6.-0.35*ZETA2)*EXP(-0.35*ZETA2)
          PSIH2(K)=0.7*ZETA2+0.75*ZETA2*(6.-0.35*ZETA2)*EXP(-0.35*ZETA2)
!----------------------------------------------------------------------
!
        ENDIF
!
!----------------------------------------------------------------------
        IF(K==KZTM)THEN
          ZTMAX1=ZETA1
          ZTMAX2=ZETA2
        ENDIF
!
        ZETA1=ZETA1+DZETA1
        ZETA2=ZETA2+DZETA2
!----------------------------------------------------------------------
      ENDDO
!----------------------------------------------------------------------
      ZTMAX1=ZTMAX1-EPS
      ZTMAX2=ZTMAX2-EPS
!----------------------------------------------------------------------
!
      END SUBROUTINE MYJSFC_INIT
!
!----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      END MODULE MODULE_SURFACE_LAYER
!
!-----------------------------------------------------------------------
