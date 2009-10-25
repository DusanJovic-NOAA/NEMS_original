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
      INCLUDE '../../inc/kind.inc'
!
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: total_integ_tim,totalsum_tim
!
!-----------------------------------------------------------------------
!***  ASSOCIATED WITH DYNAMICS
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: adv1_tim,adv2_tim,bocoh_tim,bocov_tim          &
                        ,cdwdt_tim,cdzdt_tim,consts_tim                 &
                        ,ddamp_tim,dht_tim                              &
                        ,dyn_init_tim,dyn_run_tim                       &
                        ,exch_dyn_tim                                   &
                        ,fftfhn_tim,fftfwn_tim,hadv2_tim                &
                        ,hdiff_tim,init_tim,mono_tim                    &
                        ,pdtsdt_tim,pgforce_tim,poavhn_tim              &
                        ,polehn_tim,polewn_tim                          &
                        ,prefft_tim,presmud_tim                         &
                        ,swaphn_tim,swapwn_tim                          &
                        ,update_dyn_int_state_tim,updatet_tim           &
                        ,vadv2_tim,vsound_tim,vtoa_tim
!
!-----------------------------------------------------------------------
!***  ASSOCIATED WITH PHYSICS
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: adjppt_tim,cucnvc_tim,exch_phy_tim             &
                        ,gsmdrive_tim,h_to_v_tim,gfs_phy_tim            &
                        ,phy_init_tim,phy_run_tim,phy_sum_tim           &
                        ,pole_swap_phy_tim,radiation_tim,rdtemp_tim     &
                        ,turbl_tim,update_phy_int_state_tim    
!
!-----------------------------------------------------------------------
!***  ASSOCIATED WITH DYNAMICS-PHYSICS COUPLER
!-----------------------------------------------------------------------
!
      REAL(kind=KDBL) :: add_fld_tim,cpl_dyn_phy_tim,get_fld_tim
!
!-----------------------------------------------------------------------
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
      REAL(kind=KDBL) :: FACTOR
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
#ifdef IBM
      FACTOR=1.
#else
      FACTOR=1.0E-3
#endif
!
      totalsum_tim=adv1_tim                                             &
                  +bocoh_tim                                            &
                  +bocov_tim                                            &
                  +cdwdt_tim                                            &
                  +cdzdt_tim                                            &
!jaa                  +consts_tim                                       &
                  +dht_tim                                              &
                  +ddamp_tim                                            &
                  +exch_dyn_tim                                         &
!jaa                  +init_tim                                         &
                  +fftfhn_tim                                           &
                  +fftfwn_tim                                           &
                  +hadv2_tim                                            &
                  +hdiff_tim                                            &
                  +pdtsdt_tim                                           &
                  +pgforce_tim                                          &
                  +poavhn_tim                                           &
                  +polehn_tim                                           &
                  +polewn_tim                                           &
!jaa                  +prefft_tim                                       &
!jaa                  +presmud_tim                                      &
                  +swaphn_tim                                           &
                  +swapwn_tim                                           &
                  +update_dyn_int_state_tim                             &
                  +updatet_tim                                          &
                  +vadv2_tim                                            &
                  +vsound_tim                                           &
                  +vtoa_tim
!
      totalsum_tim=totalsum_tim                                         &
                  +cucnvc_tim                                           &
                  +exch_phy_tim                                         &
                  +gsmdrive_tim                                         &
                  +h_to_v_tim                                           &
                  +pole_swap_phy_tim                                    &
                  +radiation_tim                                        &
                  +rdtemp_tim                                           &
                  +turbl_tim                                            &
                  +update_phy_int_state_tim
!
      totalsum_tim=totalsum_tim                                         &
                  +dyn_init_tim                                         &
                  +phy_init_tim                                         &
                  +cpl_dyn_phy_tim
!
!-----------------------------------------------------------------------
!***  The designated MPI task writes clocktimes for its work.
!-----------------------------------------------------------------------
!
      IF(MYPE==NPE_PRINT)THEN
!
        write(0,FMT='(" ntsd= ",I6," total_integration_tim=  ",g10.5)') ntimestep,total_integ_tim*factor
! &
!                 ,' totalsum_tim=',totalsum_tim*factor
!
        write(0,*)' DYNAMICS'
!
        write(0,FMT='("  dyn_run=               ",g10.5," pct= ",f7.2)') dyn_run_tim*factor &
                 ,dyn_run_tim/total_integ_tim*100.
        write(0,FMT='("   dyn_init=             ",g10.5," pct= ",f7.2)') dyn_init_tim*factor &
                 ,dyn_init_tim/total_integ_tim*100.
        write(0,FMT='("   update_dyn_int_state= ",g10.5," pct= ",f7.2)') update_dyn_int_state_tim*factor &
                 ,update_dyn_int_state_tim/total_integ_tim*100.
        write(0,FMT='("   consts=               ",g10.5," pct= ",f7.2)') consts_tim*factor &
                 ,consts_tim/total_integ_tim*100.
        write(0,FMT='("   init=                 ",g10.5," pct= ",f7.2)') init_tim*factor &
                 ,init_tim/total_integ_tim*100.
        write(0,FMT='("   pgforce=              ",g10.5," pct= ",f7.2)') pgforce_tim*factor &
                 ,pgforce_tim/total_integ_tim*100.
        write(0,FMT='("   dht=                  ",g10.5," pct= ",f7.2)') dht_tim*factor &
                 ,dht_tim/total_integ_tim*100.
        write(0,FMT='("   ddamp=                ",g10.5," pct= ",f7.2)') ddamp_tim*factor &
                ,ddamp_tim/total_integ_tim*100.
        write(0,FMT='("   pdtsdt=               ",g10.5," pct= ",f7.2)') pdtsdt_tim*factor &
                ,pdtsdt_tim/total_integ_tim*100.
        write(0,FMT='("   vtoa=                 ",g10.5," pct= ",f7.2)') vtoa_tim*factor &
                ,vtoa_tim/total_integ_tim*100.
        write(0,FMT='("   adv1=                 ",g10.5," pct= ",f7.2)') adv1_tim*factor &
                ,adv1_tim/total_integ_tim*100.
!
        if(vadv2_tim/=0.)then
          write(0,FMT='("   vadv2=                ",g10.5," pct= ",f7.2)') vadv2_tim*factor &
                  ,vadv2_tim/total_integ_tim*100.
        endif
!
        if(hadv2_tim/=0.)then
          write(0,FMT='("   hadv2=                ",g10.5," pct= ",f7.2)') hadv2_tim*factor &
                ,hadv2_tim/total_integ_tim*100.
        endif
!
        if(adv2_tim/=0.)then
          write(0,FMT='("   adv2=                 ",g10.5," pct= ",f7.2)') adv2_tim*factor &
               ,adv2_tim/total_integ_tim*100.
        endif
!
        if(mono_tim/=0.)then
          write(0,FMT='("   mono=                 ",g10.5," pct= ",f7.2)') mono_tim*factor &
               ,mono_tim/total_integ_tim*100.
        endif
!
        write(0,FMT='("   cdzdt=                ",g10.5," pct= ",f7.2)') cdzdt_tim*factor &
                ,cdzdt_tim/total_integ_tim*100.
        write(0,FMT='("   cdwdt=                ",g10.5," pct= ",f7.2)') cdwdt_tim*factor &
                ,cdwdt_tim/total_integ_tim*100.
        write(0,FMT='("   vsound=               ",g10.5," pct= ",f7.2)') vsound_tim*factor &
                ,vsound_tim/total_integ_tim*100. 
        write(0,FMT='("   hdiff=                ",g10.5," pct= ",f7.2)') hdiff_tim*factor &
                ,hdiff_tim/total_integ_tim*100.
        write(0,FMT='("   bocoh=                ",g10.5," pct= ",f7.2)') bocoh_tim*factor &
                ,bocoh_tim/total_integ_tim*100.
        write(0,FMT='("   bocov=                ",g10.5," pct= ",f7.2)') bocov_tim*factor &
                ,bocov_tim/total_integ_tim*100.
        write(0,FMT='("   updatet=              ",g10.5," pct= ",f7.2)') updatet_tim*factor &
                ,updatet_tim/total_integ_tim*100.

        if(prefft_tim/=0.)then
          write(0,FMT='("   prefft=               ",g10.5," pct= ",f7.2)') prefft_tim*factor &
                ,prefft_tim/total_integ_tim*100.
          write(0,FMT='("   fftfhn=               ",g10.5," pct= ",f7.2)') fftfhn_tim*factor &
                ,fftfhn_tim/total_integ_tim*100.
          write(0,FMT='("   fftfwn=               ",g10.5," pct= ",f7.2)') fftfwn_tim*factor &
                ,fftfwn_tim/total_integ_tim*100.
          write(0,FMT='("   polewn=               ",g10.5," pct= ",f7.2)') polewn_tim*factor &
                ,polewn_tim/total_integ_tim*100.
          write(0,FMT='("   poavhn=               ",g10.5," pct= ",f7.2)') poavhn_tim*factor &
                ,poavhn_tim/total_integ_tim*100.
        endif
!
        if(presmud_tim/=0.)then
          write(0,FMT='("  presmud=               ",g10.5," pct= ",f7.2)') presmud_tim*factor &
                ,presmud_tim/total_integ_tim*100.
        endif
!jaa        write(0,FMT='("  cpl_dyn_phy_tim=       ",g10.5," pct= ",f7.2)') cpl_dyn_phy_tim*factor &
!jaa                 ,cpl_dyn_phy_tim/total_integ_tim*100.
!jaa         write(0,FMT='("  get_fld_tim=           ",g10.5," pct= ",f7.2)') get_fld_tim*factor &
!jaa                 ,get_fld_tim/total_integ_tim*100.
!jaa         write(0,FMT='("  add_fld_tim=           ",g10.5," pct= ",f7.2)') add_fld_tim*factor &
!jaa                 ,add_fld_tim/total_integ_tim*100.
        write(0,*)' PHYSICS '
!
        write(0,FMT='("  phy_run=               ",g10.5," pct= ",f7.2)') phy_run_tim*factor &
                ,phy_run_tim/total_integ_tim*100.
        write(0,FMT='("   phy_init=             ",g10.5," pct= ",f7.2)') phy_init_tim*factor &
                ,phy_init_tim/total_integ_tim*100.
        write(0,FMT='("   update_phy_int_state= ",g10.5," pct= ",f7.2)') update_phy_int_state_tim*factor &
                ,update_phy_int_state_tim/total_integ_tim*100.
        write(0,FMT='("   cucnvc=               ",g10.5," pct= ",f7.2)') cucnvc_tim*factor &
                ,cucnvc_tim/total_integ_tim*100.
        write(0,FMT='("   gsmdrive=             ",g10.5," pct= ",f7.2)') gsmdrive_tim*factor &
                ,gsmdrive_tim/total_integ_tim*100.
        write(0,FMT='("   radiation=            ",g10.5," pct= ",f7.2)') radiation_tim*factor &
                ,radiation_tim/total_integ_tim*100.
        write(0,FMT='("   rdtemp=               ",g10.5," pct= ",f7.2)') rdtemp_tim*factor &
                ,rdtemp_tim/total_integ_tim*100.
        write(0,FMT='("   turbl=                ",g10.5," pct= ",f7.2)') turbl_tim*factor &
                ,turbl_tim/total_integ_tim*100.
        write(0,FMT='("   h_to_v=               ",g10.5," pct= ",f7.2)') h_to_v_tim*factor &
                ,h_to_v_tim/total_integ_tim*100.
!
        if(pole_swap_phy_tim/=0.)then
          write(0,FMT='("   pole_swap_phy=        ",g10.5," pct= ",f7.2)') pole_swap_phy_tim*factor &
                ,pole_swap_phy_tim/total_integ_tim*100.
        endif
!
        write(0,*)' EXCHANGE TIMES '
!
        write(0,FMT='("   exch_dyn=             ",g10.5," pct= ",f7.2)') exch_dyn_tim*factor &
                ,exch_dyn_tim/total_integ_tim*100.
        write(0,FMT='("   exch_phy=             ",g10.5," pct= ",f7.2)') exch_phy_tim*factor &
                ,exch_phy_tim/total_integ_tim*100.
!
        if(swaphn_tim/=0.)then
          write(0,FMT='("   swaphn=               ",g10.5," pct= ",f7.2)') swaphn_tim*factor &
                ,swaphn_tim/total_integ_tim*100.
        endif
!
        if(swapwn_tim/=0.)then
          write(0,FMT='("   swapwn=               ",g10.5," pct= ",f7.2)') swapwn_tim*factor &
                ,swapwn_tim/total_integ_tim*100.
        endif
!
        if(polehn_tim/=0.)then
          write(0,FMT='("   polehn=               ",g10.5," pct= ",f7.2)') polehn_tim*factor &
                ,polehn_tim/total_integ_tim*100.
        endif
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
