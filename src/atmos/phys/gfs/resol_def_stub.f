      module resol_def
      
      implicit none
      
      integer,public ::   jcap,jcap1,jcap2,latr,latr2
      integer,public ::   levh,levm1,levp1,levs,lnt,lnt2,lnt22,levr
      integer,public ::   lnte,lnted,lnto,lntod,lnuv
      integer,public ::   lonr,lonrx
!jw      integer   ntrac
      integer,target,public  ::   ntrac
      integer,public ::   nxpt,nypt,jintmx,latrd
      integer,public ::   ntoz,ntcw
      integer,public ::   lsoil,nmtvr,ncld,num_p3d,num_p2d,nrcm
      integer,public ::   ngrids_sfcc, ngrids_flx, nfxr
!jws      integer   ivsupa, ivssfc, ivssfc_restart, ivsinp
      integer,public ::   ivsupa, ivssfc_restart, ivsinp
      integer,target,public :: thermodyn_id, sfcpress_id                      ! hmhj
      integer,target,public :: ivssfc
      integer,target,public :: ngrids_gg
      integer,public ::   ngrids_sfcc2d,ngrids_sfcc3d
!jwe
!
      integer,public ::   nlunit
!jw      integer   thermodyn_id, sfcpress_id			! hmhj
!
      integer,public ::   g_gz, g_ps, g_t, g_u, g_v, g_q, g_p, g_dp, g_dpdt
      integer,public ::   lotgr

      integer,public :: kwq,kwte,kwdz,kwrq

!     For Ensemble concurrency run. Weiyu
!     INTEGER :: Ensemble_Id, Total_member

!     The option to add 2d/3d diag fields to physics export state
      logical :: lgocart

      end module resol_def
!
      module ozne_def
      implicit none
      
      integer, parameter ,public :: kozpl=28, kozc=48
      integer ,public :: latsozp, levozp, timeoz, latsozc, levozc, timeozc &
     ,       PL_Coeff
      real (kind=8) ,public ::blatc, dphiozc
      real (kind=8), allocatable ,public :: PL_LAT(:), PL_Pres(:) &
     ,                                     PL_TIME(:)
      end module ozne_def
