!-----------------------------------------------------------------------
!
      MODULE MODULE_CLOCKTIMES
!
!-----------------------------------------------------------------------
!
!***  List the clocktime counters for the various parts of
!***  the integration and print them as desired.
!
!-----------------------------------------------------------------------
!
      USE MODULE_INCLUDE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: INTEGRATION_TIMERS                                      &
               ,PRINT_CLOCKTIMES                                        &
               ,TIMERS                                                  &
               ,cbcst_tim,pbcst_tim
!
!-----------------------------------------------------------------------
!
      TYPE INTEGRATION_TIMERS
!
        REAL(kind=KDBL) :: total_integ_tim,totalsum_tim
!
        REAL(kind=KDBL) :: adv1_tim,adv2_tim,bocoh_tim,bocov_tim        &
                          ,cdwdt_tim,cdzdt_tim,consts_tim               &
                          ,ddamp_tim,dht_tim                            &
                          ,exch_tim                                     &
                          ,fftfhn_tim,fftfwn_tim,hadv2_tim              &
                          ,hdiff_tim,mono_tim                           &
                          ,pdtsdt_tim,pgforce_tim,poavhn_tim            &
                          ,polehn_tim,polewn_tim                        &
                          ,prefft_tim,presmud_tim                       &
                          ,solver_init_tim,solver_run_tim               &
                          ,swaphn_tim,swapwn_tim                        &
                          ,updatet_tim                                  &
                          ,vadv2_tim,vsound_tim,vtoa_tim
!
        REAL(kind=KDBL) :: adjppt_tim,cucnvc_tim                        &
                          ,gsmdrive_tim,h_to_v_tim,gfs_phy_tim          &
                          ,phy_sum_tim                                  &
                          ,pole_swap_tim,radiation_tim,rdtemp_tim       &
                          ,turbl_tim
!
        REAL(kind=KDBL) :: domain_run_1                                 &
                          ,domain_run_2                                 &
                          ,domain_run_3                                 &
                          ,pc_cpl_run_cpl1                              &
                          ,pc_cpl_run_cpl2                              &
                          ,pc_cpl_run_cpl3                              &
                          ,pc_cpl_run_cpl4                              &
                          ,cpl1_recv_tim                                &
                          ,cpl2_send_tim                                &
                          ,cpl2_comp_tim                                &
                          ,cpl2_wait_tim                                &
                          ,parent_bookkeep_moving_tim                   &
                          ,parent_update_moving_tim                     &
                          ,t0_recv_move_tim
!
!-----------------------------------------------------------------------
!***  Associated with moving nests
!-----------------------------------------------------------------------
!
        REAL(kind=KDBL) :: update_interior_from_nest_tim                &
                          ,update_interior_from_parent_tim
!
      END TYPE INTEGRATION_TIMERS
!
!-----------------------------------------------------------------------
!
      TYPE(INTEGRATION_TIMERS),DIMENSION(:),ALLOCATABLE,TARGET :: TIMERS   !<-- Timers for each domain
!
      REAL(kind=KDBL),DIMENSION(5) :: cbcst_tim,pbcst_tim
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE PRINT_CLOCKTIMES(NTIMESTEP                             &
                                 ,MY_DOMAIN_ID                          &
                                 ,MYPE                                  &
                                 ,NPE_PRINT                             &
                                 ,TIMERS_DOMAIN)
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER(kind=KINT),INTENT(IN) :: NTIMESTEP                        &  !<-- Forecast timestep
                                      ,MY_DOMAIN_ID                     &  !<-- The domain's ID
                                      ,MYPE                             &  !<-- The task ID
                                      ,NPE_PRINT                           !<-- ID of task providing clocktime diagnostics
!
      TYPE(INTEGRATION_TIMERS),TARGET,INTENT(INOUT) :: TIMERS_DOMAIN       !<-- Assorted clocktime timers for current domain
!
!---------------------
!***  Local variables
!---------------------
!
      REAL(kind=KDBL) :: FACTOR
!
      TYPE(INTEGRATION_TIMERS),POINTER :: TD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      TD=>TIMERS_DOMAIN                                                    !<-- Abbreviate name of this domain's timers
!
!-----------------------------------------------------------------------
!
#ifdef IBM
      FACTOR=1.
#else
      FACTOR=1.0E-3
#endif
!
      td%totalsum_tim=td%adv1_tim                                       &
                  +td%bocoh_tim                                         &
                  +td%bocov_tim                                         &
                  +td%cdwdt_tim                                         &
                  +td%cdzdt_tim                                         &
!jaa                  +td%consts_tim                                    &
                  +td%dht_tim                                           &
                  +td%ddamp_tim                                         &
                  +td%exch_tim                                          &
!jaa                  +td%init_tim                                      &
                  +td%fftfhn_tim                                        &
                  +td%fftfwn_tim                                        &
                  +td%hadv2_tim                                         &
                  +td%hdiff_tim                                         &
                  +td%pdtsdt_tim                                        &
                  +td%pgforce_tim                                       &
                  +td%poavhn_tim                                        &
                  +td%polehn_tim                                        &
                  +td%polewn_tim                                        &
!jaa                  +td%prefft_tim                                    &
!jaa                  +td%presmud_tim                                   &
                  +td%swaphn_tim                                        &
                  +td%swapwn_tim                                        &
                  +td%updatet_tim                                       &
                  +td%vadv2_tim                                         &
                  +td%vsound_tim                                        &
                  +td%vtoa_tim
!
      td%totalsum_tim=td%totalsum_tim                                   &
                  +td%cucnvc_tim                                        &
                  +td%gsmdrive_tim                                      &
                  +td%h_to_v_tim                                        &
                  +td%pole_swap_tim                                     &
                  +td%radiation_tim                                     &
                  +td%rdtemp_tim                                        &
                  +td%turbl_tim
!
      td%totalsum_tim=td%totalsum_tim                                   &
!xxx              +td%dyn_init_tim                                      &
!xxx              +td%phy_init_tim                                      &
!xxx              +td%cpl_dyn_phy_tim
                  +td%solver_init_tim 
!
!-----------------------------------------------------------------------
!***  The designated MPI task writes clocktimes for its work.
!-----------------------------------------------------------------------
!
      IF(MYPE==NPE_PRINT)THEN
!
        write(0,*)' '
        write(0,FMT='(" Clocktimes for domain #",I2.2)') my_domain_id
!
        write(0,FMT='(" ntsd= ",I6," total_integration_tim=  ",g10.5)') ntimestep,td%total_integ_tim*factor
! &
!                 ,' totalsum_tim=',td%totalsum_tim*factor
!
        write(0,*)' DYNAMICS'
!
        write(0,FMT='("  solver_run=            ",g10.5," pct= ",f7.2)') td%solver_run_tim*factor &
                 ,td%solver_run_tim/td%total_integ_tim*100.
        write(0,FMT='("   solver_init=          ",g10.5," pct= ",f7.2)') td%solver_init_tim*factor &
                 ,td%solver_init_tim/td%total_integ_tim*100.
!xxx    write(0,FMT='("   update_dyn_int_state= ",g10.5," pct= ",f7.2)') td%update_dyn_int_state_tim*factor &
!xxx             ,td%update_dyn_int_state_tim/td%total_integ_tim*100.
        write(0,FMT='("   consts=               ",g10.5," pct= ",f7.2)') td%consts_tim*factor &
                 ,td%consts_tim/td%total_integ_tim*100.
!       write(0,FMT='("   init=                 ",g10.5," pct= ",f7.2)') td%init_tim*factor &
!                ,td%init_tim/td%total_integ_tim*100.
        write(0,FMT='("   pgforce=              ",g10.5," pct= ",f7.2)') td%pgforce_tim*factor &
                 ,td%pgforce_tim/td%total_integ_tim*100.
        write(0,FMT='("   dht=                  ",g10.5," pct= ",f7.2)') td%dht_tim*factor &
                 ,td%dht_tim/td%total_integ_tim*100.
        write(0,FMT='("   ddamp=                ",g10.5," pct= ",f7.2)') td%ddamp_tim*factor &
                ,td%ddamp_tim/td%total_integ_tim*100.
        write(0,FMT='("   pdtsdt=               ",g10.5," pct= ",f7.2)') td%pdtsdt_tim*factor &
                ,td%pdtsdt_tim/td%total_integ_tim*100.
        write(0,FMT='("   vtoa=                 ",g10.5," pct= ",f7.2)') td%vtoa_tim*factor &
                ,td%vtoa_tim/td%total_integ_tim*100.
        write(0,FMT='("   adv1=                 ",g10.5," pct= ",f7.2)') td%adv1_tim*factor &
                ,td%adv1_tim/td%total_integ_tim*100.
!
        if(td%vadv2_tim/=0.)then
          write(0,FMT='("   vadv2=                ",g10.5," pct= ",f7.2)') td%vadv2_tim*factor &
                  ,td%vadv2_tim/td%total_integ_tim*100.
        endif
!
        if(td%hadv2_tim/=0.)then
          write(0,FMT='("   hadv2=                ",g10.5," pct= ",f7.2)') td%hadv2_tim*factor &
                ,td%hadv2_tim/td%total_integ_tim*100.
        endif
!
        if(td%adv2_tim/=0.)then
          write(0,FMT='("   adv2=                 ",g10.5," pct= ",f7.2)') td%adv2_tim*factor &
               ,td%adv2_tim/td%total_integ_tim*100.
        endif
!
        if(td%mono_tim/=0.)then
          write(0,FMT='("   mono=                 ",g10.5," pct= ",f7.2)') td%mono_tim*factor &
               ,td%mono_tim/td%total_integ_tim*100.
        endif
!
        write(0,FMT='("   cdzdt=                ",g10.5," pct= ",f7.2)') td%cdzdt_tim*factor &
                ,td%cdzdt_tim/td%total_integ_tim*100.
        write(0,FMT='("   cdwdt=                ",g10.5," pct= ",f7.2)') td%cdwdt_tim*factor &
                ,td%cdwdt_tim/td%total_integ_tim*100.
        write(0,FMT='("   vsound=               ",g10.5," pct= ",f7.2)') td%vsound_tim*factor &
                ,td%vsound_tim/td%total_integ_tim*100. 
        write(0,FMT='("   hdiff=                ",g10.5," pct= ",f7.2)') td%hdiff_tim*factor &
                ,td%hdiff_tim/td%total_integ_tim*100.
        write(0,FMT='("   bocoh=                ",g10.5," pct= ",f7.2)') td%bocoh_tim*factor &
                ,td%bocoh_tim/td%total_integ_tim*100.
        write(0,FMT='("   bocov=                ",g10.5," pct= ",f7.2)') td%bocov_tim*factor &
                ,td%bocov_tim/td%total_integ_tim*100.
        write(0,FMT='("   updatet=              ",g10.5," pct= ",f7.2)') td%updatet_tim*factor &
                ,td%updatet_tim/td%total_integ_tim*100.

        if(td%prefft_tim/=0.)then
          write(0,FMT='("   prefft=               ",g10.5," pct= ",f7.2)') td%prefft_tim*factor &
                ,td%prefft_tim/td%total_integ_tim*100.
          write(0,FMT='("   fftfhn=               ",g10.5," pct= ",f7.2)') td%fftfhn_tim*factor &
                ,td%fftfhn_tim/td%total_integ_tim*100.
          write(0,FMT='("   fftfwn=               ",g10.5," pct= ",f7.2)') td%fftfwn_tim*factor &
                ,td%fftfwn_tim/td%total_integ_tim*100.
          write(0,FMT='("   polewn=               ",g10.5," pct= ",f7.2)') td%polewn_tim*factor &
                ,td%polewn_tim/td%total_integ_tim*100.
          write(0,FMT='("   poavhn=               ",g10.5," pct= ",f7.2)') td%poavhn_tim*factor &
                ,td%poavhn_tim/td%total_integ_tim*100.
        endif
!
        if(td%presmud_tim/=0.)then
          write(0,FMT='("  presmud=               ",g10.5," pct= ",f7.2)') td%presmud_tim*factor &
                ,td%presmud_tim/td%total_integ_tim*100.
        endif
        write(0,*)' PHYSICS '
!
!xxx    write(0,FMT='("  phy_run=               ",g10.5," pct= ",f7.2)') td%phy_run_tim*factor &
!xxx            ,td%phy_run_tim/td%total_integ_tim*100.
!xxx    write(0,FMT='("   phy_init=             ",g10.5," pct= ",f7.2)') td%phy_init_tim*factor &
!xxx            ,td%phy_init_tim/td%total_integ_tim*100.
!xxx    write(0,FMT='("   update_phy_int_state= ",g10.5," pct= ",f7.2)') td%update_phy_int_state_tim*factor &
!xxx            ,td%update_phy_int_state_tim/td%total_integ_tim*100.
        write(0,FMT='("   cucnvc=               ",g10.5," pct= ",f7.2)') td%cucnvc_tim*factor &
                ,td%cucnvc_tim/td%total_integ_tim*100.
        write(0,FMT='("   gsmdrive=             ",g10.5," pct= ",f7.2)') td%gsmdrive_tim*factor &
                ,td%gsmdrive_tim/td%total_integ_tim*100.
        write(0,FMT='("   radiation=            ",g10.5," pct= ",f7.2)') td%radiation_tim*factor &
                ,td%radiation_tim/td%total_integ_tim*100.
        write(0,FMT='("   rdtemp=               ",g10.5," pct= ",f7.2)') td%rdtemp_tim*factor &
                ,td%rdtemp_tim/td%total_integ_tim*100.
        write(0,FMT='("   turbl=                ",g10.5," pct= ",f7.2)') td%turbl_tim*factor &
                ,td%turbl_tim/td%total_integ_tim*100.
        write(0,FMT='("   h_to_v=               ",g10.5," pct= ",f7.2)') td%h_to_v_tim*factor &
                ,td%h_to_v_tim/td%total_integ_tim*100.
!
        if(td%pole_swap_tim/=0.)then
          write(0,FMT='(" pole_swap=            ",g10.5," pct= ",f7.2)') td%pole_swap_tim*factor &
                ,td%pole_swap_tim/td%total_integ_tim*100.
        endif
!
        write(0,*)' EXCHANGE TIMES '
!
        write(0,FMT='("   exch_dyn=             ",g10.5," pct= ",f7.2)') td%exch_tim*factor &
                ,td%exch_tim/td%total_integ_tim*100.
!
        if(td%swaphn_tim/=0.)then
          write(0,FMT='("   swaphn=               ",g10.5," pct= ",f7.2)') td%swaphn_tim*factor &
                ,td%swaphn_tim/td%total_integ_tim*100.
        endif
!
        if(td%swapwn_tim/=0.)then
          write(0,FMT='("   swapwn=               ",g10.5," pct= ",f7.2)') td%swapwn_tim*factor &
                ,td%swapwn_tim/td%total_integ_tim*100.
        endif
!
        if(td%polehn_tim/=0.)then
          write(0,FMT='("   polehn=               ",g10.5," pct= ",f7.2)') td%polehn_tim*factor &
                ,td%polehn_tim/td%total_integ_tim*100.
        endif
!
      ENDIF
!
!-----------------------------------------------------------------------
!
!     total_integ_tim=0.
!     totalsum_tim=0
!     adv1_tim=0.
!     adv2_tim=0.
!     bocoh_tim=0.
!     bocov_tim=0.
!     cdwdt_tim=0.
!     cdzdt_tim=0.
!     consts_tim=0.
!     ddamp_tim=0.
!     dht_tim=0.
!     exch_tim=0.
!     fftfhn_tim=0.
!     fftfwn_tim=0.
!     hadv2_tim=0.
!     hdiff_tim=0.
!     mono_tim=0.
!     pdtsdt_tim=0.
!     pgforce_tim=0.
!     poavhn_tim=0.
!     polehn_tim=0.
!     polewn_tim=0.
!     prefft_tim=0.
!     presmud_tim=0.
!     solver_init_tim=0.
!     solver_run_tim=0.
!     swaphn_tim=0.
!     swapwn_tim=0.
!     updatet_tim=0.
!     vadv2_tim=0.
!     vsound_tim=0.
!     vtoa_tim=0.
!     adjppt_tim=0.
!     cucnvc_tim=0.
!     gsmdrive_tim=0.
!     h_to_v_tim=0.
!     gfs_phy_tim=0.
!     phy_sum_tim=0.
!     pole_swap_tim=0.
!     radiation_tim=0.
!     rdtemp_tim=0.
!     turbl_tim=0.
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PRINT_CLOCKTIMES
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_CLOCKTIMES
!
!-----------------------------------------------------------------------
