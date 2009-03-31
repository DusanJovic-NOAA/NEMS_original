!-----------------------------------------------------------------------
!
      MODULE MODULE_CLOCKTIMES
!
!-----------------------------------------------------------------------
!
!***  LIST THE CLOCKTIME COUNTERS FOR THE VARIOUS PARTS OF
!***  THE INTEGRATION AND PRINT THEM AS DESIRED.
!
!-----------------------------------------------------------------------
!
      REAL :: total_tim,totalsum_tim
!
!-----------------------------------------------------------------------
!***  ASSOCIATED WITH DYNAMICS
!-----------------------------------------------------------------------
!
      REAL :: adv1_tim,adv2_tim,bocoh_tim,bocov_tim                     &
             ,cdwdt_tim,cdzdt_tim,consts_tim                            &
             ,ddamp_tim,dht_tim                                         &
             ,dyn_init_tim,dyn_run_tim                                  &
             ,exch_dyn_tim                                              &
             ,fftfhn_tim,fftfwn_tim,hadv2_tim                           &
             ,hdiff_tim,init_tim,mono_tim                               &
             ,pdtsdt_tim,pgforce_tim,poavhn_tim                         &
             ,polehn_tim,polewn_tim                                     &
             ,prefft_tim,presmud_tim                                    &
             ,swaphn_tim,swapwn_tim                                     &
             ,update_dyn_int_state_tim,updatet_tim                      &
             ,vadv2_tim,vsound_tim,vtoa_tim
!
!-----------------------------------------------------------------------
!***  ASSOCIATED WITH PHYSICS
!-----------------------------------------------------------------------
!
      REAL :: cucnvc_tim,exch_phy_tim,gsmdrive_tim,h_to_v_tim           &
             ,phy_init_tim,phy_run_tim,phy_sum_tim                      &
             ,pole_swap_phy_tim,radiation_tim,rdtemp_tim                &
             ,turbl_tim,update_phy_int_state_tim
!
!-----------------------------------------------------------------------
!***  ASSOCIATED WITH DYNAMICS-PHYSICS COUPLER
!-----------------------------------------------------------------------
!
      REAL :: add_fld_tim,cpl_dyn_phy_tim,get_fld_tim
!
!-----------------------------------------------------------------------
      CHARACTER(LEN=10) :: timer_name(50)
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!
      SUBROUTINE PRINT_CLOCKTIMES(NTIMESTEP,MYPE,NPE_PRINT)
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: NTIMESTEP                                   &  !<-- Forecast timestep
                           ,MYPE                                        &  !<-- My task ID
                           ,NPE_PRINT                                      !<-- ID of task providing clocktime diagnostics
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!

      totalsum_tim=adv1_tim                                         &
                  +bocoh_tim                                        &
                  +bocov_tim                                        &
                  +cdwdt_tim                                        &
                  +cdzdt_tim                                        &
                  +consts_tim                                       &
                  +dht_tim                                          &
                  +ddamp_tim                                        &
                  +exch_dyn_tim                                     &
                  +init_tim                                         &
                  +fftfhn_tim                                       &
                  +fftfwn_tim                                       &
                  +hadv2_tim                                        &
                  +hdiff_tim                                        &
                  +pdtsdt_tim                                       &
                  +pgforce_tim                                      &
                  +poavhn_tim                                       &
                  +polehn_tim                                       &
                  +polewn_tim                                       &
                  +prefft_tim                                       &
                  +presmud_tim                                      &
                  +swaphn_tim                                       &
                  +swapwn_tim                                       &
                  +update_dyn_int_state_tim                         &
                  +updatet_tim                                      &
                  +vadv2_tim                                        &
                  +vsound_tim                                       &
                  +vtoa_tim
!
      totalsum_tim=totalsum_tim                                     &
                  +cucnvc_tim                                       &
                  +exch_phy_tim                                     &
                  +gsmdrive_tim                                     &
                  +h_to_v_tim                                       &
                  +pole_swap_phy_tim                                &
                  +radiation_tim                                    &
                  +rdtemp_tim                                       &
                  +turbl_tim                                        &
                  +update_phy_int_state_tim
!
      totalsum_tim=totalsum_tim                                     &
                  +dyn_init_tim                                     &
                  +phy_init_tim                                     &
                  +cpl_dyn_phy_tim
!
!-----------------------------------------------------------------------
!***  THE DESIGNATED MPI TASK WRITES CLOCKTIMES FOR ITS WORK.
!-----------------------------------------------------------------------
!
      IF(MYPE==NPE_PRINT)THEN
!
        write(0,*)' ntsd=',ntimestep,' total_tim=',total_tim*1.e-3 &
                 ,' totalsum_tim=',totalsum_tim*1.e-3
!
        write(0,*)' DYNAMICS'
!
        write(0,*)'  dyn_run=',dyn_run_tim*1.e-3 &
                 ,' pct=',dyn_run_tim/total_tim*100.
        write(0,*)'   dyn_init=',dyn_init_tim*1.e-3 &
                 ,' pct=',dyn_init_tim/total_tim*100.
        write(0,*)'   update_dyn_int_state=',update_dyn_int_state_tim*1.e-3 &
                 ,' pct=',update_dyn_int_state_tim/total_tim*100.
!jaa        write(0,*)' cpl_dyn_phy_tim=',cpl_dyn_phy_tim*1.e-3 &
!jaa                  ,' pct=',cpl_dyn_phy_tim/total_tim*100.
!jaa         write(0,*)' get_fld_tim=',get_fld_tim*1.e-3 &
!jaa                  ,' pct=',get_fld_tim/total_tim*100.
!jaa         write(0,*)' add_fld_tim=',add_fld_tim*1.e-3 &
!jaa                  ,' pct=',add_fld_tim/total_tim*100.
        write(0,*)'   consts=',consts_tim*1.e-3 &
                 ,' pct=',consts_tim/total_tim*100.
        write(0,*)'   init=',init_tim*1.e-3 &
                 ,' pct=',init_tim/total_tim*100.
        write(0,*)'   pgforce=',pgforce_tim*1.e-3 &
                 ,' pct=',pgforce_tim/total_tim*100.
        write(0,*)'   dht=',dht_tim*1.e-3 &
                 ,' pct=',dht_tim/total_tim*100.
        write(0,*)'   ddamp=',ddamp_tim*1.e-3 &
                 ,' pct=',ddamp_tim/total_tim*100.
        write(0,*)'   pdtsdt=',pdtsdt_tim*1.e-3 &
                 ,' pct=',pdtsdt_tim/total_tim*100.
        write(0,*)'   vtoa=',vtoa_tim*1.e-3 &
                 ,' pct=',vtoa_tim/total_tim*100.
        write(0,*)'   adv1=',adv1_tim*1.e-3 &
                 ,' pct=',adv1_tim/total_tim*100.
!
        if(vadv2_tim/=0.)then
         write(0,*)'   vadv2=',vadv2_tim*1.e-3 &
                  ,' pct=',vadv2_tim/total_tim*100.
        endif
!
        if(hadv2_tim/=0.)then
        write(0,*)'   hadv2=',hadv2_tim*1.e-3 &
                 ,' pct=',hadv2_tim/total_tim*100.
        endif
!
        if(adv2_tim/=0.)then
        write(0,*)'   adv2=',adv2_tim*1.e-3 &
                 ,' pct=',adv2_tim/total_tim*100.
        endif
!
        if(mono_tim/=0.)then
        write(0,*)'   mono=',mono_tim*1.e-3 &
                 ,' pct=',mono_tim/total_tim*100.
        endif
!
        write(0,*)'   cdzdt=',cdzdt_tim*1.e-3 &
                 ,' pct=',cdzdt_tim/total_tim*100.
        write(0,*)'   cdwdt=',cdwdt_tim*1.e-3 &
                 ,' pct=',cdwdt_tim/total_tim*100.
        write(0,*)'   vsound=',vsound_tim*1.e-3 &
                 ,' pct=',vsound_tim/total_tim*100.
        write(0,*)'   hdiff=',hdiff_tim*1.e-3 &
                 ,' pct=',hdiff_tim/total_tim*100.
        write(0,*)'   bocoh=',bocoh_tim*1.e-3 &
                 ,' pct=',bocoh_tim/total_tim*100.
        write(0,*)'   bocov=',bocov_tim*1.e-3 &
                 ,' pct=',bocov_tim/total_tim*100.
        write(0,*)'   updatet=',updatet_tim*1.e-3 &
                 ,' pct=',updatet_tim/total_tim*100.
        write(0,*)'   prefft=',prefft_tim*1.e-3 &
                 ,' pct=',prefft_tim/total_tim*100.
        write(0,*)'   presmud=',presmud_tim*1.e-3 &
                 ,' pct=',presmud_tim/total_tim*100.
        write(0,*)'   fftfhn=',fftfhn_tim*1.e-3 &
                 ,' pct=',fftfhn_tim/total_tim*100.
        write(0,*)'   fftfwn=',fftfwn_tim*1.e-3 &
                 ,' pct=',fftfwn_tim/total_tim*100.
        write(0,*)'   polewn=',polewn_tim*1.e-3 &
                 ,' pct=',polewn_tim/total_tim*100.
        write(0,*)'   poavhn=',poavhn_tim*1.e-3 &
                 ,' pct=',poavhn_tim/total_tim*100.
!
        write(0,*)' PHYSICS '
!
        write(0,*)'  phy_run =',phy_run_tim*1.e-3 &
                 ,' pct=',phy_run_tim/total_tim*100.
        write(0,*)'   phy_init=',phy_init_tim*1.e-3 &
                 ,' pct=',phy_init_tim/total_tim*100.
        write(0,*)'   update_phy_int_state=',update_phy_int_state_tim*1.e-3 &
                 ,' pct=',update_phy_int_state_tim/total_tim*100.
        write(0,*)'   exch_phy=',exch_phy_tim*1.e-3 &
                 ,' pct=',exch_phy_tim/total_tim*100.
        write(0,*)'   cucnvc=',cucnvc_tim*1.e-3 &
                 ,' pct=',cucnvc_tim/total_tim*100.
        write(0,*)'   gsmdrive=',gsmdrive_tim*1.e-3 &
                 ,' pct=',gsmdrive_tim/total_tim*100.
        write(0,*)'   radiation=',radiation_tim*1.e-3 &
                 ,' pct=',radiation_tim/total_tim*100.
        write(0,*)'   rdtemp=',rdtemp_tim*1.e-3 &
                 ,' pct=',rdtemp_tim/total_tim*100.
        write(0,*)'   turbl=',turbl_tim*1.e-3 &
                 ,' pct=',turbl_tim/total_tim*100.
        write(0,*)'   h_to_v=',h_to_v_tim*1.e-3 &
                 ,' pct=',h_to_v_tim/total_tim*100.
        write(0,*)'   pole_swap_phy=',pole_swap_phy_tim*1.e-3 &
                 ,' pct=',pole_swap_phy_tim/total_tim*100.
!
        write(0,*)' EXCHANGE TIMES '
!
        write(0,*)'  swaphn=',swaphn_tim*1.e-3 &
                 ,' pct=',swaphn_tim/total_tim*100.
        write(0,*)'  swapwn=',swapwn_tim*1.e-3 &
                 ,' pct=',swapwn_tim/total_tim*100.
        write(0,*)'  polehn=',polehn_tim*1.e-3 &
                 ,' pct=',polehn_tim/total_tim*100.
        write(0,*)'  exch_dyn=',exch_dyn_tim*1.e-3 &
                  ,' pct=',exch_dyn_tim/total_tim*100.
!
      ENDIF
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
