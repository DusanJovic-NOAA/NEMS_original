!!!!!  ==========================================================  !!!!!
!!!!!            sw-gfdl1 radiation package description            !!!!!
!!!!!  ==========================================================  !!!!!
!                                                                      !
!   the sw-gfdl1 package includes these parts:                         !
!                                                                      !
!      'radsw_gfdl1_param.f'                                           !
!      'radsw_gfdl1_datatb.f'                                          !
!      'radsw_gfdl1_main.f'                                            !
!                                                                      !
!   the 'radsw_gfdl1_param.f' contains:                                !
!                                                                      !
!      'module_radsw_parameters'  -- band parameters set up            !
!      'module_radsw_cntr_para'   -- control parameters set up         !
!                                                                      !
!   the 'radsw_gfdl1_datatb.f' contains:                               !
!                                                                      !
!      'module_radsw_bandtbl'     -- band structure and data           !
!      'module_radsw_cldprtb'     -- cloud property coefficients table !
!                                                                      !
!   the 'radsw_gfdl1_main.f' contains the main module:                 !
!                                                                      !
!      'module_radsw_main'        -- main sw radiation transfer        !
!                                                                      !
!   in the main module 'module_radsw_main' there are only two          !
!   externally callable subroutines:                                   !
!                                                                      !
!      'swrad'      -- main gfdl1 sw radiation routine                 !
!      'rswinit'    -- initialization routine                          !
!                                                                      !
!   all the sw radiation subprograms become contained subprograms      !
!   in module 'module_radsw_main' and many of them are not directly    !
!   accessable from places outside the module.                         !
!                                                                      !
!   compilation sequence is:                                           !
!                                                                      !
!      'radsw_gfdl1_param.f'                                           !
!      'radsw_gfdl1_datatb.f'                                          !
!      'radsw_gfdl1_main.f'                                            !
!                                                                      !
!   and all should be put in front of routines that use sw modules     !
!                                                                      !
!!!!!  ==========================================================  !!!!!
!!!!!                       end descriptions                       !!!!!
!!!!!  ==========================================================  !!!!!



!========================================!
      module module_radsw_bandtbl        !
!........................................!
!
      use machine,                 only : kind_phys
      use module_radsw_parameters, only : NBANDS, NSTREAMS, NFRQPTS,    &
     &                                    NINTSOLAR, NH2OBANDS
!
      implicit   none

!  ---   gaussian weights for two, four, and eight-point quadratures
      real (kind=kind_phys) :: gwt2(2), gwt4(4), gwt8(8)

      data gwt2 / 0.5, 0.5 /
      data gwt4 / .1739274226, .3260725774, .3260725774, .1739274226 /
      data gwt8 / .0506142681, .1111905172, .1568533229, .1813418917,   &
     &            .1813418917, .1568533229, .1111905172, .0506142681 /

!  ---   solar flux in each parameterization band
      real (kind=kind_phys), dimension(NBANDS) :: solflxband

      data solflxband   /          12.1587,  6.5070, 10.7300, 23.8226,  &
     &  19.2689, 43.7116, 35.7886,135.0955,239.2806,222.9263,138.7890,  &
     & 182.3105,101.2186, 72.2298, 48.5104, 28.2587, 15.4827,  6.0424,  &
     &   3.7148,  3.0384,  1.7734,  1.9695,  3.1789,  1.0869,  1.0672 /

!  ---  wavenumber value for the band limits
      integer, dimension(0:NBANDS) :: endwvnbands

      data endwvnbands  /          0,   2500,   2900,   3400,   4200,   &
     &          4700,   5600,   6200,   8200,  11500,  14600,  16700,   &
     &         20000,  22300,  24600,  27500,  30000,  31900,  33000,   &
     &         33800,  34500,  35300,  36500,  40000,  43300,  57600 /

!  ---  number of pseudo-monochromatic frequencies in each band
      integer, dimension(NBANDS) :: nfreqpts

      data nfreqpts     /      6,  1,  5,  6,  1,  9,  5,  9,  7,  8,   &
     &     1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1  /

!  ---  
      real (kind=kind_phys), dimension(NSTREAMS) :: wtstr

      data wtstr / 0.347854845, 0.652145155, 0.347854845, 0.652145155 /

!  ---  weight associated with each exponential term
      real (kind=kind_phys), dimension(NFRQPTS) :: wtfreq

      data wtfreq / 9.312504E-2, 1.626168E-1, 1.433068E-1, 3.071999E-1, &
     & 2.737514E-1, 2.000000E-2, 1.000000E+0, 7.874545E-2, 2.476534E-1, &
     & 3.794250E-1, 2.755277E-1, 1.864845E-2, 1.921193E-1, 2.751345E-1, &
     & 2.480785E-1, 2.262776E-1, 1.839013E-2, 4.000000E-2, 1.000000E+0, &
     & 3.837498E-2, 8.791451E-2, 1.637795E-1, 1.189705E-1, 2.020078E-1, &
     & 1.509591E-1, 1.693291E-1, 6.466451E-2, 4.000000E-3, 1.605849E-2, &
     & 8.165883E-2, 6.004498E-2, 3.943264E-1, 4.479113E-1, 4.476276E-2, &
     & 1.119772E-1, 1.180453E-1, 6.352831E-2, 8.756711E-2, 5.027652E-2, &
     & 1.652990E-1, 3.485438E-1, 1.000000E-2, 5.315826E-3, 2.793918E-2, &
     & 4.138399E-2, 1.803244E-1, 1.955244E-1, 1.832192E-1, 3.662930E-1, &
     & 3.601121E-3, 1.291813E-2, 4.997123E-2, 8.941965E-2, 1.016408E-1, &
     & 1.407492E-1, 2.304239E-1, 3.712796E-1, 1.000000E+0, 1.000000E+0, &
     & 1.000000E+0, 1.000000E+0, 1.000000E+0, 1.000000E+0, 1.000000E+0, &
     & 1.000000E+0, 1.000000E+0, 1.000000E+0, 1.000000E+0, 1.000000E+0, &
     & 1.000000E+0, 1.000000E+0, 1.000000E+0 /

!  ---  number of wavenumbers in each region where the solar flux is constant
      integer, dimension(NINTSOLAR) :: nwvnsolar

      data nwvnsolar/ 100,  11,  14,  18,  24,  33,  50,  83,  12,  12, &
     &  13,  15,  15,  17,  18,  20,  21,  24,  26,  30,  32,  37,  42, &
     &  47,  55,  64,  76,  91, 111, 139, 179, 238, 333,  41,  42,  45, &
     &  46,  48,  51,  53,  55,  58,  61,  64,  68,  71,  75,  79,  84, &
     &  89,  95, 101, 107, 115, 123, 133, 142, 154, 167, 181, 197, 217, &
     & 238, 263, 293, 326, 368, 417, 476, 549, 641, 758, 909, 101, 103, &
     & 105, 108, 109, 112, 115, 117, 119, 122, 125, 128, 130, 134, 137, &
     & 140, 143, 147, 151, 154, 158, 163, 166, 171, 175, 181, 185, 190, &
     & 196, 201, 207, 213, 219, 227, 233, 240, 248, 256, 264, 274, 282, &
     & 292, 303, 313, 325, 337, 349, 363, 377, 392, 408, 425, 444, 462, &
     & 483, 505, 529, 554, 580, 610, 641, 675, 711, 751, 793, 841, 891, &
     & 947,1008,1075,1150,1231,1323,1425,1538,1667,1633,14300 /

!  ---  solar flux in watts per meter**2 in each wavenumber region where it is constant
      real (kind=kind_phys), dimension(NINTSOLAR) :: solint

      data  solint(  1: 50)       /                                     &
     &     1.60000E-6, 2.88000E-5, 3.60000E-5, 4.59200E-5, 6.13200E-5,  &
     &     8.55000E-5, 1.28600E-4, 2.16000E-4, 2.90580E-4, 3.10184E-4,  &
     &     3.34152E-4, 3.58722E-4, 3.88050E-4, 4.20000E-4, 4.57056E-4,  &
     &     4.96892E-4, 5.45160E-4, 6.00600E-4, 6.53600E-4, 7.25040E-4,  &
     &     7.98660E-4, 9.11200E-4, 1.03680E-3, 1.18440E-3, 1.36682E-3,  &
     &     1.57560E-3, 1.87440E-3, 2.25500E-3, 2.74500E-3, 3.39840E-3,  &
     &     4.34000E-3, 5.75400E-3, 7.74000E-3, 9.53050E-3, 9.90192E-3,  &
     &     1.02874E-2, 1.06803E-2, 1.11366E-2, 1.15830E-2, 1.21088E-2,  &
     &     1.26420E-2, 1.32250E-2, 1.38088E-2, 1.44612E-2, 1.51164E-2,  &
     &     1.58878E-2, 1.66500E-2, 1.75140E-2, 1.84450E-2, 1.94106E-2 /
      data  solint( 51:100)       /                                     &
     &     2.04864E-2, 2.17248E-2, 2.30640E-2, 2.44470E-2, 2.59840E-2,  &
     &     2.75940E-2, 2.94138E-2, 3.13950E-2, 3.34800E-2, 3.57696E-2,  &
     &     3.84054E-2, 4.13490E-2, 4.46880E-2, 4.82220E-2, 5.22918E-2,  &
     &     5.70078E-2, 6.19888E-2, 6.54720E-2, 6.69060E-2, 6.81226E-2,  &
     &     6.97788E-2, 7.12668E-2, 7.27100E-2, 7.31610E-2, 7.33471E-2,  &
     &     7.34814E-2, 7.34717E-2, 7.35072E-2, 7.34939E-2, 7.35202E-2,  &
     &     7.33249E-2, 7.31713E-2, 7.35462E-2, 7.36920E-2, 7.23677E-2,  &
     &     7.25023E-2, 7.24258E-2, 7.20766E-2, 7.18284E-2, 7.32757E-2,  &
     &     7.31645E-2, 7.33277E-2, 7.36128E-2, 7.33752E-2, 7.28965E-2,  &
     &     7.24924E-2, 7.23307E-2, 7.21050E-2, 7.12620E-2, 7.10903E-2 /
      data  solint(101:151)       /                        7.12714E-2,  &
     &     7.08012E-2, 7.03752E-2, 7.00350E-2, 6.98639E-2, 6.90690E-2,  &
     &     6.87621E-2, 6.52080E-2, 6.65184E-2, 6.60038E-2, 6.47615E-2,  &
     &     6.44831E-2, 6.37206E-2, 6.24102E-2, 6.18698E-2, 6.06320E-2,  &
     &     5.83498E-2, 5.67028E-2, 5.51232E-2, 5.48645E-2, 5.12340E-2,  &
     &     4.85581E-2, 4.85010E-2, 4.79220E-2, 4.44058E-2, 4.48718E-2,  &
     &     4.29373E-2, 4.15242E-2, 3.81744E-2, 3.16342E-2, 2.99615E-2,  &
     &     2.92740E-2, 2.67484E-2, 1.76904E-2, 1.40049E-2, 1.46224E-2,  &
     &     1.39993E-2, 1.19574E-2, 1.06386E-2, 1.00980E-2, 8.63808E-3,  &
     &     6.52736E-3, 4.99410E-3, 4.39350E-3, 2.21676E-3, 1.33812E-3,  &
     &     1.12320E-3, 5.59000E-4, 3.60000E-4, 2.98080E-4, 7.46294E-5  /

!  ---  gaussian points and weights for evaluation of the diffuse beam
      real (kind=kind_phys), dimension(NSTREAMS) :: ptstr

      data ptstr  / -0.861136312,-0.339981044,0.861136312,0.339981044 /

!  ---  scaling factor used in the fit of the h2o transmission function
      real (kind=kind_phys), dimension(NH2OBANDS) :: powph2o

      data powph2o    / 0.84, 0.00, 0.96, 0.60, 0.00, 0.88, 0.91,       &
     &                  0.80, 0.49, 0.38, 0.00, 0.00, 0.00, 0.00 /

!  ---  reference pressure (Pa) used in the fit of the h2o transmission function
      real (kind=kind_phys), dimension(NH2OBANDS) :: p0h2o

      data p0h2o      /                                                 &
     &     300., 101325.,  50000.,   5000., 101325.,   4000.,  70000.,  &
     &  101325., 101325.,  70000., 101325., 101325., 101325., 101325. /

!  ---  c(n)co2(str) = coefficients for the absorptivity expression for co2
!                      for the pressure-scaled and non-scaled, respectively,
!                      portions of the fit (=1.0e-99 if no absorption)
      real (kind=kind_phys), dimension(NH2OBANDS) :: c1co2, c2co2,      &
     &       c3co2, c1co2str, c2co2str, c3co2str

      data c1co2      /                                                 &
     &   5.4E+02, 1.0E-99, 5.0E-02, 1.6E-01, 1.0E-99, 1.3E-02, 2.4E-06, &
     &   3.3E-04, 2.8E+00, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /
      data c2co2      /                                                 &
     &   7.2E-04, 1.0E-99, 1.0E+03, 1.0E-02, 1.0E-99, 4.2E-01, 9.9E+03, &
     &   6.6E+00, 5.8E+05, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /
      data c3co2      /                                                 &
     &   1.8E-05, 1.0E-99, 6.8E-02, 1.0E-01, 1.0E-99, 4.1E-01, 9.9E-01, &
     &   5.6E-01, 2.8E-02, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /

      data c1co2str   /                                                 &
     &   7.2E-03, 1.0E-99, 1.3E-04, 7.1E-03, 1.0E-99, 3.2E-02, 6.6E-01, &
     &   2.8E-04, 7.0E-02, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /
      data c2co2str   /                                                 &
     &   1.0E-99, 1.0E-99, 1.0E+02, 4.4E-03, 1.0E-99, 6.3E-01, 7.2E+02, &
     &   9.5E+00, 1.0E+03, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /
      data c3co2str   /                                                 &
     &   4.0E-01, 1.0E-99, 5.0E-01, 9.7E-02, 1.0E-99, 3.3E-02, 4.9E-03, &
     &   4.6E-01, 9.4E-03, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99 /

!  --- c(n)o2(str) = coefficients for the absorptivity expression for o2
!                    for the pressure-scaled and non-scaled, respectively,
!                    portions of the fit (=1.0e-99 if no absorption)
      real (kind=kind_phys), dimension(NH2OBANDS) :: c1o2, c2o2,        &
     &       c3o2, c1o2str, c2o2str, c3o2str

      data c1o2       /                                                 &
     &   1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, &
     &   3.4E-02, 1.0E-99, 3.9E-03, 1.4E-07, 1.0E-99, 1.0E-99, 1.0E-99 /
      data c2o2       /                                                 &
     &   1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, &
     &   4.5E+05, 1.0E-99, 3.3E+03, 5.2E+08, 1.0E-99, 1.0E-99, 1.0E-99 /
      data c3o2       /                                                 &
     &   1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, &
     &   8.4E-02, 1.0E-99, 2.1E-01, 7.8E-01, 1.0E-99, 1.0E-99, 1.0E-99 /

      data c1o2str    /                                                 &
     &   1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, &
     &   7.9E-06, 1.0E-99, 1.9E-03, 4.9E-06, 1.0E-99, 1.0E-99, 1.0E-99 /
      data c2o2str    /                                                 &
     &   1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, &
     &   9.8E+03, 1.0E-99, 1.5E+02, 7.0E+04, 1.0E-99, 1.0E-99, 1.0E-99 /
      data c3o2str    /                                                 &
     &   1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, 1.0E-99, &
     &   4.7E-01, 1.0E-99, 9.7E-02, 4.6E-01, 1.0E-99, 1.0E-99, 1.0E-99 /

!  ---  c(n)o2strschrun = coefficients for the absorptivity expression for the
!                         Schuman-Runge o2 band (non-scaled only)
      real (kind=kind_phys) :: c1o2strschrun,c2o2strschrun,c3o2strschrun

      data  c1o2strschrun, c2o2strschrun, c3o2strschrun                 &
     &      / 2.7E-02,       1.2E-01,       1.5E-01 /

!  ---  kh2o =  the psuedo-absorption coefficients in m2/kg for h2o
      real (kind=kind_phys), dimension(NFRQPTS) :: kh2o

      data kh2o   /              4.000000E-1, 2.546100E-2, 1.489350E-3, &
     & 5.228619E-5, 0.000000E+0, 1.000000E+2, 3.524000E-3, 4.000000E+0, &
     & 2.923541E-1, 3.820742E-2, 4.629296E-3, 0.000000E+0, 5.000000E+0, &
     & 4.756542E-1, 3.740251E-2, 2.046388E-3, 0.000000E+0, 1.000000E+2, &
     & 1.443000E-3, 1.000000E+1, 6.131070E-1, 8.325077E-2, 2.336348E-2, &
     & 5.893526E-3, 4.474983E-4, 3.560266E-5, 0.000000E+0, 4.000000E+2, &
     & 4.000000E-1, 2.855524E-2, 5.188692E-3, 6.285804E-4, 0.000000E+0, &
     & 1.000000E+1, 1.124617E+0, 2.061909E-1, 6.405800E-2, 2.300286E-2, &
     & 8.792879E-3, 1.722225E-3, 0.000000E+0, 5.000000E+1, 5.000000E+0, &
     & 7.355553E-1, 2.057103E-1, 4.088761E-2, 6.036039E-3, 1.437820E-3, &
     & 0.000000E+0, 3.000000E-1, 7.057908E-2, 1.566582E-2, 3.665539E-3, &
     & 1.463247E-3, 5.343198E-4, 1.353690E-4, 0.000000E+0, 1.925000E-4, &
     & 2.025000E-4, 3.664000E-6, 6.239000E-6, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0 /

!  ---  ko3  = the absorption coefficients in m2/kg for o3
      real (kind=kind_phys), dimension(NFRQPTS) :: ko3

      data ko3    /              0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, 0.000000E+0, &
     & 0.000000E+0, 5.041992E-1, 5.041992E-1, 5.041992E-1, 5.041992E-1, &
     & 5.041992E-1, 5.041992E-1, 5.041992E-1, 5.041992E-1, 3.745368E+0, &
     & 4.021207E+0, 7.062558E-1, 9.620166E-2, 0.000000E+0, 9.738058E-1, &
     & 2.378498E+1, 1.567498E+2, 5.394976E+2, 1.207743E+3, 2.706951E+3, &
     & 5.277161E+3, 1.177364E+4, 1.035882E+4, 2.475921E+3 /

!  ---  strterm = logical flag to indicate whether or not a h2o pseudo-
!                 absorption coefficient is assigned a non-scaled
!                 (true) or pressure-scaled (false) gas amount
      logical, dimension(NFRQPTS) :: strterm

      data strterm  /                  .FALSE.,.FALSE.,.FALSE.,.FALSE., &
     & .FALSE., .TRUE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., &
     & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., .TRUE.,.FALSE.,.FALSE., &
     & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., .TRUE., &
     & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., &
     & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., .TRUE.,.FALSE.,.FALSE., &
     & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., &
     & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., &
     & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE., &
     & .FALSE.,.FALSE.,.FALSE.,.FALSE.  / 

!........................................!
      end module module_radsw_bandtbl    !
!========================================!



!========================================!
      module module_radsw_cldprtb        !
!........................................!
!
!     use machine,                 only : kind_phys
      use module_radsw_parameters, only : NLIQCLDV,  NICECLDV,          &
     &                                    NRAINCLDV, NSNOWCLDV
!
      implicit none

!  ---  define the spectral limits for drop, ice, rain, and snow single
!       scattering properties in shortwave frequency ranges. 

      integer, dimension(NLIQCLDV)   :: endliqwvn
      integer, dimension(NICECLDV)   :: endicewvn
      integer, dimension(NRAINCLDV)  :: endrainwvn
      integer, dimension(NSNOWCLDV)  :: endsnowwvn

!  ---  wavenumber limits for slingo cloud drop intervals (=24)
      data endliqwvn  /                       2924, 3437, 4202, 4695,   &
     &    6098, 6536, 7813, 8404, 9091,10000,11494,12821,13333,14493,   &
     &   15625,17544,19231,20833,22727,25000,27778,30303,33333,57600 /

!  ---  wavenumber limits for fu ice crystal intervals (=25)
      data endicewvn  /                 2000, 2924, 3437, 4202, 4695,   &
     &    6098, 6536, 7092, 8404, 9091,10000,11494,12821,13333,14493,   &
     &   15625,17544,19231,20833,22727,25000,27778,30303,33333,57600 /

!  ---  wavenumber limits for Savijarvi rain drop intervals (=4)
      data endrainwvn /  4202, 8403, 14493, 57600 /

!  ---  wavenumber limits for the Fu snow intervals (=6)
      data endsnowwvn /  2857, 4000, 5263, 7692, 14493, 57600 /


!........................................!
      end module module_radsw_cldprtb    !
!========================================!

