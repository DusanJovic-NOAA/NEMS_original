      module resol_def
      
      implicit none
      
      integer   jcap,jcap1,jcap2,latg,latg2,latr,latr2
      integer   levh,levm1,levp1,levs,lnt,lnt2,lnt22,levr
      integer   lnte,lnted,lnto,lntod,lnuv
      integer   lonf,lonfx,lonr,lonrx
      integer   ntrac
      integer   nxpt,nypt,jintmx,latrd
      integer   ntoz,ntcw
      integer   lsoil,nmtvr,ncld,num_p3d,num_p2d,nrcm
      integer   ngrids_sfcc, ngrids_flx, nfxr, ngrids_gg
      integer   ivsupa, ivssfc, ivssfc_restart, ivsinp
      integer   nlunit
      integer   thermodyn_id, sfcpress_id			! hmhj
!
      integer   g_gz, g_ps, g_t, g_u, g_v, g_q, g_p, g_dp, g_dpdt
      integer   lotgr

      integer kwq,kwte,kwdz,kwrq

!     For Ensemble concurrency run. Weiyu
!     INTEGER :: Ensemble_Id, Total_member

      end module resol_def
!
      module ozne_def
      use machine , only : kind_phys
      implicit none
      
      integer, parameter :: kozpl=28, kozc=48
      integer latsozp, levozp, timeoz, latsozc, levozc, timeozc
     &,       PL_Coeff
      real (kind=kind_phys) blatc, dphiozc
      real (kind=kind_phys), allocatable :: PL_LAT(:), PL_Pres(:)
     &,                                     PL_TIME(:)
      end module ozne_def
