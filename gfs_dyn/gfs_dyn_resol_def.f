      module gfs_dyn_resol_def
      use gfs_dyn_machine
      implicit none
      save
      integer   jcap,jcap1,jcap2,latg,latg2
      integer   levh,levm1,levp1,levs,lnt,lnt2,lnt22,levr
      integer   lnte,lnted,lnto,lntod,lnuv
      integer   lonf,lonfx
      integer   ntrac
      integer   nxpt,nypt,jintmx,latgd
      integer   ntoz,ntcw,ncld
      integer   ngrids_gg
      integer   ivsupa, ivsinp
      integer   nlunit
      integer   thermodyn_id, sfcpress_id			! hmhj
!
      INTEGER   P_GZ,P_ZEM,P_DIM,P_TEM,P_RM,P_QM
      INTEGER   P_zslam,P_zsphi
      INTEGER   P_ZE,P_DI,P_TE,P_RQ,P_Q,P_DLAM,P_DPHI,P_ULN,P_VLN
      INTEGER   P_W,P_X,P_Y,P_RT,P_ZQ
      INTEGER   G_UUM,G_VVM,G_TTM,G_RM,G_QM,G_GZ
      INTEGER   G_UU ,G_VV ,G_TT ,G_RQ,G_Q 
      INTEGER   G_U  ,G_V  ,G_T  ,G_RT,G_ZQ, g_p, g_dp, g_dpdt
      INTEGER   LOTS,LOTD,LOTA,LOTLS,LOTGR

      integer   kwq,kwte,kwdz,kwrq

      integer   ksz, ksd, kst, ksr, ksq, ksplam, kspphi
      integer   ksu, ksv, kzslam, kzsphi
!
      integer   kau, kav, kat, kar, kaps, kazs
!
      integer   kdtphi, kdrphi, kdtlam, kdrlam 
      integer   kdulam, kdvlam, kduphi, kdvphi


      end module gfs_dyn_resol_def
