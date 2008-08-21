!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright 2002, 2003, Atmospheric & Environmental Research, Inc. (AER).
! This software may be used, copied, or redistributed as long as it is
! not sold and this copyright notice is reproduced on each copy made.
! This model is provided as is without any express or implied warranties.
!                      (http://www.rtweb.aer.com/)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  ==========================================================  !!!!!
!!!!!              rrtm radiation package description              !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the rrtm package includes three parts:                            !
!                                                                      !
!       'radlw_rrtm_param.f'                                           !
!       'radlw_rrtm_datatb.f'                                          !
!       'radlw_rrtm_main.f'                                            !
!                                                                      !
!    the 'radlw_rrtm_param.f' contains:                                !
!                                                                      !
!       'module_radlw_parameters'  -- band parameters set up           !
!       'module_radlw_cntr_para'   -- control parameters set up        !
!                                                                      !
!    the 'radlw_rrtm_datatb.f' contains:                               !
!                                                                      !
!       'module_radlw_aerosols'    -- aerosols data tables (not yet!)  !
!       'module_radlw_avplank'     -- plank flux data                  !
!       'module_radlw_cldprlw'     -- cloud property coefficients      !
!       'module_radlw_kgbnn'       -- absorption coeffients for 16     !
!                                     bands, where nn = 01-16          !
!                                                                      !
!    the 'radlw_rrtm_main.f' contains the main module:                 !
!                                                                      !
!       'module_radlw_main'                                            !
!                                                                      !
!    in the main module 'module_radlw_main' there are only two         !
!    externally callable subroutines:                                  !
!                                                                      !
!                                                                      !
!       'lwrad'     -- main rrtm lw radiation routine                  !
!       'rlwinit'   -- to initialize rrtm lw radiation                 !
!                                                                      !
!    all the lw radiation subprograms become contained subprograms     !
!    in module 'module_radlw_rrtm' and many of them are not directly   !
!    accessable from places outside the module.                        !
!                                                                      !
!    exterior modules referenced:                                      !
!                                                                      !
!       'module machine'                    in 'machine.f'             !
!       'module physcons'                   in 'physcons.f'            !
!aer    'module module_aerosols'            in 'rad_aerosols.f'        !
!                                                                      !
!    compilation sequence is:                                          !
!                                                                      !
!       'radlw_rrtm_param,f'                                           !
!       'radlw_rrtm_datatb,f'                                          !
!       'radlw_rrtm_main.f'                                            !
!                                                                      !
!    and all should be put in front of routines that use lw modules    !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!
 
 
 
!========================================!
      module module_radlw_parameters     !
!........................................!
!
!  ---  parameter constants for lw band structures
!
      implicit none
!
      integer  :: NBANDS, NGPT, N5000, N200, MAXGAS, MAXXSEC, NPLNK
      parameter (NBANDS=16, NGPT=140, N5000=5000, N200=200, NPLNK=181)
      parameter (MAXGAS=6, MAXXSEC=4)
 
!  ---  number of g-point in each band
      integer  :: NG01, NG02, NG03, NG04, NG05, NG06, NG07, NG08,       &
     &            NG09, NG10, NG11, NG12, NG13, NG14, NG15, NG16
      parameter (NG01=08, NG02=14, NG03=16, NG04=14, NG05=16, NG06=08,  &
     &           NG07=12, NG08=08, NG09=12, NG10=06, NG11=08, NG12=08,  &
     &           NG13=04, NG14=02, NG15=02, NG16=02)
 
!  ---  begining index of each band
      integer  :: NS01, NS02, NS03, NS04, NS05, NS06, NS07, NS08,       &
     &            NS09, NS10, NS11, NS12, NS13, NS14, NS15, NS16
      parameter (NS01=00, NS02=08, NS03=22, NS04=38, NS05=52, NS06=68,  &
     &           NS07=76, NS08=88, NS09=96, NS10=108, NS11=114,         &
     &           NS12=122, NS13=130, NS14=134, NS15=136, NS16=138)
 
!........................................!
      end module module_radlw_parameters !
!========================================!
 
 
 
!========================================!
      module module_radlw_cntr_para      !
!........................................!
!
        implicit   none
!
        integer :: ilwrate, iaerlw, icfclw, iflagliq, iflagice
 
!
!  ---  set up control parameters for lw radiation
!
        parameter ( ilwrate=1 )     !===> ... lw heating rate unit selection
                        !(default)  ! =1: output in k/day
                                    ! =2: output in k/second
 
        parameter ( iaerlw=0 )      !===> ... control flag for aerosols (not yet)
                        !(default)  ! =0: do not include aerosol effect
                                    ! =1: include aerosol effect
 
        parameter ( icfclw=0  )     !===> ... control flag for cfc gases
                                    ! =0: do not include cfc gases
                        !(default)  ! =1: include all cfc gases
 
        parameter ( iflagliq=3 )    !===> ... liq-cloud optical properties contrl flag
                                    ! =0: input cloud opt depth, ignor iflagice setting
                                    ! =1: input cwp,cip, (ccm2 method) ignor iflagice setting
                                    ! =2: input cwp rew, ccm3 method for liquid clouds
                        !(default)  ! =3: input cwp rew, hu and stamnes(1993) method for liq cld
 
        parameter ( iflagice=1 )    !===> ... ice-cloud optical properties contrl flag
                                    !         only used when iflagliq .ge. 2, else is ignored
                                    ! =0: input cip rei, ccm3 method for ice clouds
                        !(default)  ! =1: input cip rei, ebert and curry(1997) for ice clouds
                                    ! =2: input cip rei, streamer (1996) for ice clouds
 
!
!........................................!
      end module module_radlw_cntr_para  !
!========================================!
