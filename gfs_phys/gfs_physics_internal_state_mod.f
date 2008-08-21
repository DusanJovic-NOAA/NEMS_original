!
! !module: gfs_physics_internal_state_mod 
!                         --- internal state definition of the
!                             esmf gridded component of the gfs physics.
!
! !description:  define the gfs physics internal state used to
!                                             create the esmf internal state.
!---------------------------------------------------------------------------
! !revision history:
!
!  november 2004      weiyu yang initial code.
!  may 2005           weiyu yang for the updated gfs version.
!  february 2006      shrinivas moorthi updated for the new version of gfs
!  january 2007       hann-ming henry juang for gfs dynamics only
!  july    2007       shrinivas moorthi for gfs physics only
!  november 2007       hann-ming henry juang continue for gfs physics
!
! !interface:
!
      module gfs_physics_internal_state_mod

!!uses:
!------
      use esmf_mod
      use gfs_physics_namelist_mod

      use machine, only: kind_grid, kind_io4, kind_phys, kind_rad
      use namelist_physics_def
      USE namelist_soilveg
      use tracer_const

      use resol_def
      use layout1
      use gg_def
      use date_def
      use mpi_def

      use module_ras , only : nrcmax
      use ozne_def
      use d3d_def
      use gfs_physics_sfc_flx_mod

      implicit none

! -----------------------------------------------
      type gfs_physics_internal_state		! start type define
! -----------------------------------------------

      type(nam_gfs_phy_namelist)   :: nam_gfs_phy
      type(gfs_phy_state_namelist) :: esmf_sta_list

      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld

      integer                   :: me, nodes
      INTEGER                   :: llgg_s, lonr_s, latr_s
!     integer                   :: grib_inp

!
      integer ntrac,nxpt,nypt,jintmx,jcap,levs,lonr,latr,lats_node_r_max
      integer ntoz, ntcw, ncld, lsoil, nmtvr, num_p3d, num_p2d,levr
      integer thermodyn_id, sfcpress_id

      character(16)                     ::  cfhour1

      integer                           ::  nblck, kdt
      real                              ::  deltim

      integer              ,allocatable ::      lonsperlar (:)
      integer              ,allocatable ::  lats_nodes_r   (:)
      integer              ,allocatable ::  global_lats_r  (:)
      integer              ,allocatable ::  lats_nodes_ext (:)
      integer              ,allocatable ::  global_lats_ext(:)

      real(kind=kind_evod) ,allocatable ::      grid_gr(:,:)
      integer   g_gz ,g_ps ,g_t ,g_u ,g_v ,g_q ,g_p ,g_dp ,g_dpdt
      integer   lotgr

      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: XLON(:,:),XLAT(:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: COSZDG(:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: sfalb(:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: CLDCOV(:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: HPRIME(:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: SWH(:,:,:,:),HLW(:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: FLUXR(:,:,:)
!!

      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: phy_f3d(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: phy_f2d(:,:,:)
!
! carry fhour and initial date, may not be necessary later
      real(kind=kind_evod) ,allocatable :: fhour_idate(:,:)
      real phour
      INTEGER :: KFHOUR
      real, allocatable  :: poz(:),ozplin(:,:,:,:)
!     FOR OZON INTERPOLATION:
      INTEGER,ALLOCATABLE:: JINDX1(:),JINDX2(:)
!
      REAL,ALLOCATABLE:: DDY(:)
      REAL(KIND=KIND_RAD) SLAG,SDEC,CDEC

!!
! for nasa ozon production and distruction rates:(input throu fixio_r)
      integer 	lev,levmax
!
      integer              init,jcount,jpt,node
      integer              ibmsign
      integer              lon_dim,ilat

      real(kind=kind_evod) colat1
!!
!     real(kind=kind_evod) rone
!     real(kind=kind_evod) rlons_lat
!     real(kind=kind_evod) scale_ibm


!     integer              ibrad,ifges,ihour,ini,j,jdt,ksout,maxstp
!     integer              mdt,idt,timetot,timer,time0
!     integer              mods,n1,n2,n3,n4,ndgf,ndgi,nfiles,nflps
!     integer              nges,ngpken,niter,nnmod,nradf,nradr
!     integer              nsfcf,nsfci,nsfcs,nsigi,nsigs,nstep
!     integer              nznlf,nznli,nznls,id,iret,nsout

      integer              iret, n3, n4

      integer              ierr,iprint,k,l,locl,n
      integer              lan,lat

      real(kind=kind_phys) chour
      real(kind=kind_phys) zhour

!     logical start_step
!     logical end_step
      logical lsout, lssavd, lsoutd

      integer ikey,nrank_all,kcolor

      real(kind=kind_phys) cons0p5,cons1200,cons3600    !constant
      real(kind=kind_phys) cons0                        !constant

!
! -----------------------------------------------------
      end type gfs_physics_internal_state		! end type define
! -----------------------------------------------------

! this state is supported by c pointer not f90 pointer, thus
! need this wrap.
!-----------------------------------------------------------
      type gfs_phy_wrap		! begin type define
          type (gfs_physics_internal_state), pointer :: int_state
      end type gfs_phy_wrap	! end type define

      end module gfs_physics_internal_state_mod
