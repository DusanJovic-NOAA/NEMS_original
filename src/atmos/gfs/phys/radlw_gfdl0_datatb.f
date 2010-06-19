!!!!!  ==========================================================  !!!!!
!!!!!            gfdl0 radiation package description               !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!    the gfdl0 package includes these parts:                           !
!                                                                      !
!       'radlw_gfdl0_param.f'                                          !
!       'radlw_gfdl0_datatb.f'                                         !
!       'radlw_gfdl0_main.f'                                           !
!                                                                      !
!    the 'radlw_gfdl0_param.f' contains:                               !
!                                                                      !
!       'module_radlw_cntr_para'   -- control parameters set up        !
!       'module_radlw_parameters'  -- band parameters set up           !
!                                                                      !
!    the 'radlw_gfdl0_datatb.f' contains:                              !
!                                                                      !
!       'module_radlw_banddata'    -- data for each lw spectral band   !
!       'module_radlw_cldprlw'     -- cloud property coefficients      !
!                                                                      !
!    the 'radlw_gfdl0_main.f' contains:                                !
!                                                                      !
!       'module_radlw_main'        -- main lw radiation transfer       !
!                                                                      !
!    in the main module 'module_radlw_main' there are only two         !
!    externally callable subroutines:                                  !
!                                                                      !
!                                                                      !
!       'lwrad'     -- main rrtm lw radiation routine                  !
!       'rlwinit'   -- initialization routine                          !
!                                                                      !
!    all the lw radiation subprograms become contained subprograms     !
!    in module 'module_radlw_main' and many of them are not directly   !
!    accessable from places outside the module.                        !
!                                                                      !
!                                                                      !
!    compilation sequence is:                                          !
!                                                                      !
!       'radlw_gfdl0_param.f'                                          !
!       'radlw_gfdl0_datatb.f'                                         !
!       'radlw_gfdl0_main.f'                                           !
!                                                                      !
!    and all should be put in front of routines that use lw modules    !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radlw_banddata       !
!........................................!

      use machine,                 only : kind_phys
      use module_radlw_parameters, only : NBLW, NBLY
!
      implicit none
!
      private
!  ---  random band parameters for the lw calculations using 10 cm-1 wide
!       bands.the 15 um co2 complex is 2 bands,560-670 and 670-800 cm-1.
!       ozone coefficients are in 3 bands,670-800 (14.1 um),990-1070 and
!       1070-1200 (9.6 um). the  (nblw=163) bands now include:
!
!                56 bands, 10  cm-1 wide    0  -   560  cm-1
!                 2 bands, 15 um complex  560  -   670  cm-1
!                                         670  -   800  cm-1
!                 3 "continuum" bands     800  -   900  cm-1
!                                         900  -   990  cm-1
!                                        1070  -   1200 cm-1
!                 1 band for 9.6 um band  990  -   1070 cm-1
!               100 bands, 10 cm-1 wide  1200  -   2200 cm-1
!                 1 band for 4.3 um src  2270  -   2380 cm-1

!  ---  all bands are arranged in order of increasing wavenumber, data for
!       arndm, brndm, ao3rnd, bo3rnd are obtained by using the afgl 1982
!       catalog. continuum coefficients are from roberts (1976).

!  ---  arndm,brndm   =  random "a" and "b" parameter for (NBLW) bands

      real (kind=kind_phys), public, dimension(NBLW) :: arndm, brndm

      data arndm (  1: 93)  /    0.354693e+00,0.269857e+03,0.167062e+03,&
     & 0.201314e+04,0.964533e+03,0.547971e+04,0.152933e+04,0.599429e+04,&
     & 0.699329e+04,0.856721e+04,0.962489e+04,0.233348e+04,0.127091e+05,&
     & 0.104383e+05,0.504249e+04,0.181227e+05,0.856480e+03,0.136354e+05,&
     & 0.288635e+04,0.170200e+04,0.209761e+05,0.126797e+04,0.110096e+05,&
     & 0.336436e+03,0.491663e+04,0.863701e+04,0.540389e+03,0.439786e+04,&
     & 0.347836e+04,0.130557e+03,0.465332e+04,0.253086e+03,0.257387e+04,&
     & 0.488041e+03,0.892991e+03,0.117148e+04,0.125880e+03,0.458852e+03,&
     & 0.142975e+03,0.446355e+03,0.302887e+02,0.394451e+03,0.438112e+02,&
     & 0.348811e+02,0.615503e+02,0.143165e+03,0.103958e+02,0.725108e+02,&
     & 0.316628e+02,0.946456e+01,0.542675e+02,0.351557e+02,0.301797e+02,&
     & 0.381010e+01,0.126319e+02,0.548010e+01,0.600199e+01,0.640803e+00,&
     & 0.501549e-01,0.167961e-01,0.178110e-01,0.170166e+00,0.273514e-01,&
     & 0.983767e+00,0.753946e+00,0.941763e-01,0.970547e+00,0.268862e+00,&
     & 0.564373e+01,0.389794e+01,0.310955e+01,0.128235e+01,0.196414e+01,&
     & 0.247113e+02,0.593435e+01,0.377552e+02,0.305173e+02,0.852479e+01,&
     & 0.116780e+03,0.101490e+03,0.138939e+03,0.324228e+03,0.683729e+02,&
     & 0.471304e+03,0.159684e+03,0.427101e+03,0.114716e+03,0.106190e+04,&
     & 0.294607e+03,0.762948e+03,0.333199e+03,0.830645e+03,0.162512e+04/

      data arndm ( 94:163)  /    0.525676e+03,0.137739e+04,0.136252e+04,&
     & 0.147164e+04,0.187196e+04,0.131118e+04,0.103975e+04,0.621637e+01,&
     & 0.399459e+02,0.950648e+02,0.943161e+03,0.526821e+03,0.104150e+04,&
     & 0.905610e+03,0.228142e+04,0.806270e+03,0.691845e+03,0.155237e+04,&
     & 0.192241e+04,0.991871e+03,0.123907e+04,0.457289e+02,0.146146e+04,&
     & 0.319382e+03,0.436074e+03,0.374214e+03,0.778217e+03,0.140227e+03,&
     & 0.562540e+03,0.682685e+02,0.820292e+02,0.178779e+03,0.186150e+03,&
     & 0.383864e+03,0.567416e+01,0.225129e+03,0.473099e+01,0.753149e+02,&
     & 0.233689e+02,0.339802e+02,0.108855e+03,0.380016e+02,0.151039e+01,&
     & 0.660346e+02,0.370165e+01,0.234169e+02,0.440206e+00,0.615283e+01,&
     & 0.304077e+02,0.117769e+01,0.125248e+02,0.142652e+01,0.241831e+00,&
     & 0.483721e+01,0.226357e-01,0.549835e+01,0.597067e+00,0.404553e+00,&
     & 0.143584e+01,0.294291e+00,0.466273e+00,0.156048e+00,0.656185e+00,&
     & 0.172727e+00,0.118349e+00,0.141598e+00,0.588581e-01,0.919409e-01,&
     & 0.155521e-01,0.537083e-02 /

      data brndm (  1: 93)  /    0.789571e-01,0.920256e-01,0.696960e-01,&
     & 0.245544e+00,0.188503e+00,0.266127e+00,0.271371e+00,0.330917e+00,&
     & 0.190424e+00,0.224498e+00,0.282517e+00,0.130675e+00,0.212579e+00,&
     & 0.227298e+00,0.138585e+00,0.187106e+00,0.194527e+00,0.177034e+00,&
     & 0.115902e+00,0.118499e+00,0.142848e+00,0.216869e+00,0.149848e+00,&
     & 0.971585e-01,0.151532e+00,0.865628e-01,0.764246e-01,0.100035e+00,&
     & 0.171133e+00,0.134737e+00,0.105173e+00,0.860832e-01,0.148921e+00,&
     & 0.869234e-01,0.106018e+00,0.184865e+00,0.767454e-01,0.108981e+00,&
     & 0.123094e+00,0.177287e+00,0.848146e-01,0.119356e+00,0.133829e+00,&
     & 0.954505e-01,0.155405e+00,0.164167e+00,0.161390e+00,0.113287e+00,&
     & 0.714720e-01,0.741598e-01,0.719590e-01,0.140616e+00,0.355356e-01,&
     & 0.832779e-01,0.128680e+00,0.983013e-01,0.629660e-01,0.643346e-01,&
     & 0.717082e-01,0.629730e-01,0.875182e-01,0.857907e-01,0.358808e+00,&
     & 0.178840e+00,0.254265e+00,0.297901e+00,0.153916e+00,0.537774e+00,&
     & 0.267906e+00,0.104254e+00,0.400723e+00,0.389670e+00,0.263701e+00,&
     & 0.338116e+00,0.351528e+00,0.267764e+00,0.186419e+00,0.238237e+00,&
     & 0.210408e+00,0.176869e+00,0.114715e+00,0.173299e+00,0.967770e-01,&
     & 0.172565e+00,0.162085e+00,0.157782e+00,0.886832e-01,0.242999e+00,&
     & 0.760298e-01,0.164248e+00,0.221428e+00,0.166799e+00,0.312514e+00/

      data brndm ( 94:163)  /    0.380600e+00,0.353828e+00,0.269500e+00,&
     & 0.254759e+00,0.285408e+00,0.159764e+00,0.721058e-01,0.170528e+00,&
     & 0.231595e+00,0.307184e+00,0.564136e-01,0.159884e+00,0.147907e+00,&
     & 0.185666e+00,0.183567e+00,0.182482e+00,0.230650e+00,0.175348e+00,&
     & 0.195978e+00,0.255323e+00,0.198517e+00,0.195500e+00,0.208356e+00,&
     & 0.309603e+00,0.112011e+00,0.102570e+00,0.128276e+00,0.168100e+00,&
     & 0.177836e+00,0.105533e+00,0.903330e-01,0.126036e+00,0.101430e+00,&
     & 0.124546e+00,0.221406e+00,0.137509e+00,0.911365e-01,0.724508e-01,&
     & 0.795788e-01,0.137411e+00,0.549175e-01,0.787714e-01,0.165544e+00,&
     & 0.136484e+00,0.146729e+00,0.820496e-01,0.846211e-01,0.785821e-01,&
     & 0.122527e+00,0.125359e+00,0.101589e+00,0.155756e+00,0.189239e+00,&
     & 0.999086e-01,0.480993e+00,0.100233e+00,0.153754e+00,0.130780e+00,&
     & 0.136136e+00,0.159353e+00,0.156634e+00,0.272265e+00,0.186874e+00,&
     & 0.192090e+00,0.135397e+00,0.131497e+00,0.127463e+00,0.227233e+00,&
     & 0.190562e+00,0.214005e+00 /

!  ---  bandlo  =  lowest frequency in each of (NBLW) freq. bands
!       bandhi  =  highest frequency in each of (NBLW) freq. bands
      real (kind=kind_phys), public, dimension(NBLW) :: bandlo, bandhi

      data bandlo(  1: 93)  /     .000000e+00, .100000e+02, .200000e+02,&
     &  .300000e+02, .400000e+02, .500000e+02, .600000e+02, .700000e+02,&
     &  .800000e+02, .900000e+02, .100000e+03, .110000e+03, .120000e+03,&
     &  .130000e+03, .140000e+03, .150000e+03, .160000e+03, .170000e+03,&
     &  .180000e+03, .190000e+03, .200000e+03, .210000e+03, .220000e+03,&
     &  .230000e+03, .240000e+03, .250000e+03, .260000e+03, .270000e+03,&
     &  .280000e+03, .290000e+03, .300000e+03, .310000e+03, .320000e+03,&
     &  .330000e+03, .340000e+03, .350000e+03, .360000e+03, .370000e+03,&
     &  .380000e+03, .390000e+03, .400000e+03, .410000e+03, .420000e+03,&
     &  .430000e+03, .440000e+03, .450000e+03, .460000e+03, .470000e+03,&
     &  .480000e+03, .490000e+03, .500000e+03, .510000e+03, .520000e+03,&
     &  .530000e+03, .540000e+03, .550000e+03, .560000e+03, .670000e+03,&
     &  .800000e+03, .900000e+03, .990000e+03, .107000e+04, .120000e+04,&
     &  .121000e+04, .122000e+04, .123000e+04, .124000e+04, .125000e+04,&
     &  .126000e+04, .127000e+04, .128000e+04, .129000e+04, .130000e+04,&
     &  .131000e+04, .132000e+04, .133000e+04, .134000e+04, .135000e+04,&
     &  .136000e+04, .137000e+04, .138000e+04, .139000e+04, .140000e+04,&
     &  .141000e+04, .142000e+04, .143000e+04, .144000e+04, .145000e+04,&
     &  .146000e+04, .147000e+04, .148000e+04, .149000e+04, .150000e+04/

      data bandlo( 94:163)  /     .151000e+04, .152000e+04, .153000e+04,&
     &  .154000e+04, .155000e+04, .156000e+04, .157000e+04, .158000e+04,&
     &  .159000e+04, .160000e+04, .161000e+04, .162000e+04, .163000e+04,&
     &  .164000e+04, .165000e+04, .166000e+04, .167000e+04, .168000e+04,&
     &  .169000e+04, .170000e+04, .171000e+04, .172000e+04, .173000e+04,&
     &  .174000e+04, .175000e+04, .176000e+04, .177000e+04, .178000e+04,&
     &  .179000e+04, .180000e+04, .181000e+04, .182000e+04, .183000e+04,&
     &  .184000e+04, .185000e+04, .186000e+04, .187000e+04, .188000e+04,&
     &  .189000e+04, .190000e+04, .191000e+04, .192000e+04, .193000e+04,&
     &  .194000e+04, .195000e+04, .196000e+04, .197000e+04, .198000e+04,&
     &  .199000e+04, .200000e+04, .201000e+04, .202000e+04, .203000e+04,&
     &  .204000e+04, .205000e+04, .206000e+04, .207000e+04, .208000e+04,&
     &  .209000e+04, .210000e+04, .211000e+04, .212000e+04, .213000e+04,&
     &  .214000e+04, .215000e+04, .216000e+04, .217000e+04, .218000e+04,&
     &  .219000e+04, .227000e+04 /

      data bandhi(  1: 93)  /     .100000e+02, .200000e+02, .300000e+02,&
     &  .400000e+02, .500000e+02, .600000e+02, .700000e+02, .800000e+02,&
     &  .900000e+02, .100000e+03, .110000e+03, .120000e+03, .130000e+03,&
     &  .140000e+03, .150000e+03, .160000e+03, .170000e+03, .180000e+03,&
     &  .190000e+03, .200000e+03, .210000e+03, .220000e+03, .230000e+03,&
     &  .240000e+03, .250000e+03, .260000e+03, .270000e+03, .280000e+03,&
     &  .290000e+03, .300000e+03, .310000e+03, .320000e+03, .330000e+03,&
     &  .340000e+03, .350000e+03, .360000e+03, .370000e+03, .380000e+03,&
     &  .390000e+03, .400000e+03, .410000e+03, .420000e+03, .430000e+03,&
     &  .440000e+03, .450000e+03, .460000e+03, .470000e+03, .480000e+03,&
     &  .490000e+03, .500000e+03, .510000e+03, .520000e+03, .530000e+03,&
     &  .540000e+03, .550000e+03, .560000e+03, .670000e+03, .800000e+03,&
     &  .900000e+03, .990000e+03, .107000e+04, .120000e+04, .121000e+04,&
     &  .122000e+04, .123000e+04, .124000e+04, .125000e+04, .126000e+04,&
     &  .127000e+04, .128000e+04, .129000e+04, .130000e+04, .131000e+04,&
     &  .132000e+04, .133000e+04, .134000e+04, .135000e+04, .136000e+04,&
     &  .137000e+04, .138000e+04, .139000e+04, .140000e+04, .141000e+04,&
     &  .142000e+04, .143000e+04, .144000e+04, .145000e+04, .146000e+04,&
     &  .147000e+04, .148000e+04, .149000e+04, .150000e+04, .151000e+04/

      data bandhi( 94:163)  /     .152000e+04, .153000e+04, .154000e+04,&
     &  .155000e+04, .156000e+04, .157000e+04, .158000e+04, .159000e+04,&
     &  .160000e+04, .161000e+04, .162000e+04, .163000e+04, .164000e+04,&
     &  .165000e+04, .166000e+04, .167000e+04, .168000e+04, .169000e+04,&
     &  .170000e+04, .171000e+04, .172000e+04, .173000e+04, .174000e+04,&
     &  .175000e+04, .176000e+04, .177000e+04, .178000e+04, .179000e+04,&
     &  .180000e+04, .181000e+04, .182000e+04, .183000e+04, .184000e+04,&
     &  .185000e+04, .186000e+04, .187000e+04, .188000e+04, .189000e+04,&
     &  .190000e+04, .191000e+04, .192000e+04, .193000e+04, .194000e+04,&
     &  .195000e+04, .196000e+04, .197000e+04, .198000e+04, .199000e+04,&
     &  .200000e+04, .201000e+04, .202000e+04, .203000e+04, .204000e+04,&
     &  .205000e+04, .206000e+04, .207000e+04, .208000e+04, .209000e+04,&
     &  .210000e+04, .211000e+04, .212000e+04, .213000e+04, .214000e+04,&
     &  .215000e+04, .216000e+04, .217000e+04, .218000e+04, .219000e+04,&
     &  .220000e+04, .238000e+04 /

!  ---  ao3rnd, bo3rnd =  random "a" and "b" parameter for ozone band 2.
      real (kind=kind_phys), public :: ao3rnd, bo3rnd

      data ao3rnd, bo3rnd / 0.234676e+04,  0.922424e+01 /

!  ---  random band parameters using combined wide frequency bands between 160
!       and 1200 cm-1,as well as the 2270-2380 band for source calc.
!          bands 1-8: combined wide frequency bands for 160-560 cm-1
!          bands 9-14: freq bands, as in bandta (narrow bands) for 560-1200 cm-1
!          band  15:  freq band 2270-2380 cm-1,used for source calculation only
!       bands are arranged in order of increasing wavenumber

!  ---  data for acomb, bcomb, apcm, bpcm, atpcm, btpcm are obtained by using
!       the afgl 1982 catalog. continuum coefficients are from roberts (1976).
!       iband index values are obtained by experimentation.

!  ---  acomb, bcomb   = random "a" an "b" parameters
      real (kind=kind_phys), public, dimension(NBLY) :: acomb, bcomb

      data acomb  /                                                     &
     &  .152070e+05, .332194e+04, .527177e+03, .163124e+03, .268808e+03,&
     &  .534591e+02, .268071e+02, .123133e+02, .600199e+01, .640803e+00,&
     &  .501549e-01, .167961e-01, .178110e-01, .170166e+00, .537083e-02/

      data bcomb  /                                                     &
     &  .152538e+00, .118677e+00, .103660e+00, .100119e+00, .127518e+00,&
     &  .118409e+00, .904061e-01, .642011e-01, .629660e-01, .643346e-01,&
     &  .717082e-01, .629730e-01, .875182e-01, .857907e-01, .214005e+00/

!  ---  apcm, bpcm   = capphi coefficients
      real (kind=kind_phys), public, dimension(NBLY) :: apcm, bpcm

      data apcm   /                                                     &
     & -.671879e-03, .654345e-02, .143657e-01, .923593e-02, .117022e-01,&
     &  .159596e-01, .181600e-01, .145013e-01, .170062e-01, .233303e-01,&
     &  .256735e-01, .274745e-01, .279259e-01, .197002e-01, .349782e-01/

      data bpcm   /                                                     &
     & -.113520e-04,-.323965e-04,-.448417e-04,-.230779e-04,-.361981e-04,&
     & -.145117e-04, .198349e-04,-.486529e-04,-.550050e-04,-.684057e-04,&
     & -.447093e-04,-.778390e-04,-.982953e-04,-.772497e-04,-.748263e-04/

!  ---  atpcm, btpcm   = cappsi coefficients
      real (kind=kind_phys), public, dimension(NBLY) :: atpcm, btpcm

      data atpcm  /                                                     &
     & -.106346e-02, .641531e-02, .137362e-01, .922513e-02, .136162e-01,&
     &  .169791e-01, .206959e-01, .166223e-01, .171776e-01, .229724e-01,&
     &  .275530e-01, .302731e-01, .281662e-01, .199525e-01, .370962e-01/

      data btpcm  /                                                     &
     & -.735731e-05,-.294149e-04,-.505592e-04,-.280894e-04,-.492972e-04,&
     & -.341508e-04,-.362947e-04,-.250487e-04,-.521369e-04,-.746260e-04,&
     & -.744124e-04,-.881905e-04,-.933645e-04,-.664045e-04,-.115290e-03/

!  ---  betacm   = continuum coeff
      real (kind=kind_phys), public, dimension(NBLY) :: betacm

      data betacm /                                                     &
     &  .000000e+00, .000000e+00, .000000e+00, .000000e+00, .188625e+03,&
     &  .144293e+03, .174098e+03, .909366e+02, .497489e+02, .221212e+02,&
     &  .113124e+02, .754174e+01, .589554e+01, .495227e+01, .000000e+00/

!  ---  random band parameters for specific wide bands. at present, the
!       information consists of
!          1) random model parameters for the 15 um band,560-800 cm-1;
!          2) the continuum coefficient for the 800-990,1070-1200 cm-1 band
!       data for awide, bwide, are obtained by using the afgl 1982 catalog.
!       continuum coefficients are from roberts (1976).

      real (kind=kind_phys),public,parameter :: awide  = 0.309801e+01  !random "a" parameter
      real (kind=kind_phys),public,parameter :: bwide  = 0.495357e-01  !random "b" parameter
      real (kind=kind_phys),public,parameter :: betawd = 0.347839e+02  !continuum coefficients
      real (kind=kind_phys),public,parameter :: betinw = 0.766811e+01  !continuum coeff. for a
                                 ! specified wide freq.band (800-990 and 1070-1200 cm-1)
      real (kind=kind_phys),public,parameter :: ab15wd = awide*bwide   !15 um band complex of co2
      real (kind=kind_phys),public,parameter :: skc1r  = betawd/betinw !for 15 um band cont. coeff
      real (kind=kind_phys),public,parameter :: sko3r  = 5.89554/betinw!for 9.6 um cont. coeff
      real (kind=kind_phys),public,parameter :: sko2d  = 1.0/betinw    !cont. coeff

!

!........................................!
      end module module_radlw_banddata   !
!========================================!
