      MODULE module_n_microphysics_gfs
!
      USE n_MACHINE , ONLY : kind_phys
      USE n_FUNCPHYS
      USE n_PHYSCONS, CP => con_CP, RD => con_RD, RV => con_RV            &
     &,             T0C => con_T0C, HVAP => con_HVAP, HFUS => con_HFUS  &
     &,             EPS => con_EPS, EPSM1 => con_EPSM1                  &
     &,             EPS1 => con_FVirt, pi => con_pi, grav => con_g 
!!!CARLOS
      !!WE USE MASSI AND SDENS FROM NMMB MICROPHYSICS.... FAST FIX
      USE MODULE_MP_ETANEW, ONLY : MASSI,SDENS

      implicit none
!
!--- Common block of constants used in column microphysics
!
      real,private ::  ABFR, CBFR, CIACW, CIACR, C_N0r0,                 &
     &CN0r0, CN0r_DMRmin, CN0r_DMRmax, CRACW, CRAUT, ESW0,               &
     &QAUTx, RFmax,        RQR_DR1, RQR_DR2, RQR_DR3, RQR_DRmin,         &
     &RQR_DRmax, RR_DRmin, RR_DR1, RR_DR2, RR_DR3, RR_DRmax
!
      real,private :: mic_step
!
!--- Common block for lookup table used in calculating growth rates of
!    nucleated ice crystals growing in water saturated conditions
!--- Discretized growth rates of small ice crystals after their nucleation
!     at 1 C intervals from -1 C to -35 C, based on calculations by Miller
!     and Young (1979, JAS) after 600 s of growth.  Resultant growth rates
!     are multiplied by physics time step in GSMCONST.
!
      INTEGER, PRIVATE,PARAMETER :: MY_T1=1, MY_T2=35
      REAL,PRIVATE,DIMENSION(MY_T1:MY_T2) :: MY_GROWTH
!
!--- Parameters for ice lookup tables, which establish the range of mean ice particle
!      diameters; from a minimum mean diameter of 0.05 mm (DMImin) to a
!      maximum mean diameter of 1.00 mm (DMImax).  The tables store solutions
!      at 1 micron intervals (DelDMI) of mean ice particle diameter.
!
      REAL, PRIVATE,PARAMETER :: DMImin=.05e-3, DMImax=1.e-3,             &
     &      DelDMI=1.e-6,XMImin=1.e6*DMImin, XMImax=1.e6*DMImax
      INTEGER, PRIVATE,PARAMETER :: MDImin=XMImin, MDImax=XMImax
!
!--- Various ice lookup tables
!
      REAL, PRIVATE,DIMENSION(MDImin:MDImax) ::                           &
!     &      ACCRI,MASSI,SDENS,VSNOWI,VENTI1,VENTI2
!!!CARLOS
     &      ACCRI,VSNOWI,VENTI1,VENTI2
!
!--- Mean rain drop diameters varying from 50 microns (0.05 mm) to 450 microns
!      (0.45 mm), assuming an exponential size distribution.
!
      REAL, PRIVATE,PARAMETER :: DMRmin=.05e-3, DMRmax=.45e-3,            &
     &      DelDMR=1.e-6,XMRmin=1.e6*DMRmin, XMRmax=1.e6*DMRmax           &
     &,     NLImin=100.
      INTEGER, PRIVATE,PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax
!
    !--- Factor of 1.5 for RECImin, RESNOWmin, & RERAINmin accounts for
    !    integrating exponential distributions for effective radius
    !    (i.e., the r**3/r**2 moments).
    !
      INTEGER, PRIVATE, PARAMETER :: INDEXSmin=100
      REAL, PRIVATE, PARAMETER :: RERAINmin=1.5*XMRmin                  &
     &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=10.

!
!--- Various rain lookup tables
!--- Rain lookup tables for mean rain drop diameters from DMRmin to DMRmax,
!      assuming exponential size distributions for the rain drops
!
      REAL, PRIVATE,DIMENSION(MDRmin:MDRmax)::                          &
     &      ACCRR,MASSR,RRATE,VRAIN,VENTR1,VENTR2
!
!--- Common block for riming tables
!--- VEL_RF - velocity increase of rimed particles as functions of crude
!      particle size categories (at 0.1 mm intervals of mean ice particle
!      sizes) and rime factor (different values of Rime Factor of 1.1**N,
!      where N=0 to Nrime).
!
      INTEGER, PRIVATE,PARAMETER :: Nrime=40
      REAL, DIMENSION(2:9,0:Nrime),PRIVATE :: VEL_RF
!
!--- The following variables are for microphysical statistics
!
      INTEGER, PARAMETER :: ITLO=-60, ITHI=40
      INTEGER  NSTATS(ITLO:ITHI,4)
      REAL     QMAX(ITLO:ITHI,5),  QTOT(ITLO:ITHI,22)
!
      REAL, PRIVATE,  PARAMETER ::                                      &
     &  T_ICE=-40., T_ICE_init=-15.     !- Ver2
!
!     Some other miscellaneous parameters
!
      REAL, PRIVATE, PARAMETER :: Thom=T_ICE, TNW=50., TOLER=1.0E-20    &
! Assume fixed cloud ice effective radius
     &, RECICE=RECImin                                                  &
     &, EPSQ=1.0E-20                                                    &
     &, FLG0P1=0.1, FLG0P2=0.2, FLG1P0=1.0                              
!
      CONTAINS
!
!#######################################################################
!------- Initialize constants & lookup tables for microphysics ---------
!#######################################################################
!



!-----------------------------------
      subroutine n_rsipath2                                               &
!...................................

!  ---  inputs:
     &     ( plyr, plvl, tlyr, qlyr, qcwat, qcice, qrain, rrime,        &
     &       IM, LEVS, iflip, flgmin,                                   &
!  ---  outputs:
     &       cwatp, cicep, rainp, snowp, recwat, rerain, resnow, snden  &
     &     )

! =================   subprogram documentation block   ================ !
!                                                                       !
! abstract:  this program is a modified version of ferrier's original   !
!   "rsipath" subprogram.  it computes layer's cloud liquid, ice, rain, !
!   and snow water condensate path and the partical effective radius    !
!   for liquid droplet, rain drop, and snow flake.                      !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input variables:                                                      !
!   plyr  (IM,LEVS) : model layer mean pressure in mb (100Pa)           !
!   plvl  (IM,LEVS+1):model level pressure in mb (100Pa)                !
!   tlyr  (IM,LEVS) : model layer mean temperature in k                 !
!   qlyr  (IM,LEVS) : layer specific humidity in gm/gm                  !
!   qcwat (IM,LEVS) : layer cloud liquid water condensate amount        !
!   qcice (IM,LEVS) : layer cloud ice water condensate amount           !
!   qrain (IM,LEVS) : layer rain drop water amount                      !
!   rrime (IM,LEVS) : mass ratio of total to unrimed ice ( >= 1 )       !
!   IM              : horizontal dimention                              !
!   LEVS            : vertical layer dimensions                         !
!   iflip           : control flag for in/out vertical indexing         !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   flgmin          : Minimum large ice fraction                        !
!   lprnt           : logical check print control flag                  !
!                                                                       !
! output variables:                                                     !
!   cwatp (IM,LEVS) : layer cloud liquid water path                     !
!   cicep (IM,LEVS) : layer cloud ice water path                        !
!   rainp (IM,LEVS) : layer rain water path                             !
!   snowp (IM,LEVS) : layer snow water path                             !
!   recwat(IM,LEVS) : layer cloud eff radius for liqid water (micron)   !
!   rerain(IM,LEVS) : layer rain water effective radius      (micron)   !
!   resnow(IM,LEVS) : layer snow flake effective radius      (micron)   !
!   snden (IM,LEVS) : 1/snow density                                    !
!                                                                       !
!                                                                       !
! usage:     call rsipath2                                              !
!                                                                       !
! subroutines called:  none                                             !
!                                                                       !
! program history log:                                                  !
!      xx-xx-2001   b. ferrier     - original program                   !
!      xx-xx-2004   s. moorthi     - modified for use in gfs model      !
!      05-20-2004   y. hou         - modified, added vertical index flag!
!                     to reduce data flipping, and rearrange code to    !
!                     be comformable with radiation part programs.      !
!                                                                       !
!  ====================    end of description    =====================  !
!

      implicit none

!  ---  constant parameter:
      real, parameter :: CEXP= 1.0/3.0

!  ---  inputs:
      real, dimension(:,:), intent(in) ::                               &
     &       plyr, plvl, tlyr, qlyr, qcwat, qcice, qrain, rrime

      integer, intent(in) :: IM, LEVS, iflip
      real, dimension(:),   intent(in) :: flgmin
!     logical, intent(in) :: lprnt

!  ---  output:
      real, dimension(:,:), intent(out) ::                              &
     &       cwatp, cicep, rainp, snowp, recwat, rerain, resnow, snden

!  ---  locals:
!     real,    dimension(IM,LEVS) :: delp, pp1, pp2

      real    :: recw1, dsnow, qsnow, qqcice, flarge, xsimass, pfac,    &
     &           nlice, xli, nlimax, dum, tem,                          &
     &           rho, cpath, rc, totcnd, tc

      integer :: i, k, indexs, ksfc, k1
!
!===>  ...  begin here
!
      recw1 = 620.3505 / TNW**CEXP         ! cloud droplet effective radius

      do k = 1, LEVS
        do i = 1, IM
                                           !--- hydrometeor's optical path
           cwatp(i,k) = 0.0
           cicep(i,k) = 0.0
           rainp(i,k) = 0.0
           snowp(i,k) = 0.0
           snden(i,k) = 0.0
                                           !--- hydrometeor's effective radius
           recwat(i,k) = RECWmin
           rerain(i,k) = RERAINmin
           resnow(i,k) = RESNOWmin
        enddo
      enddo

!  ---  set up pressure related arrays, convert unit from mb to cb (10Pa)
!       cause the rest part uses cb in computation

      if (iflip == 0) then        ! data from toa to sfc
        ksfc = levs + 1
        k1   = 0
      else                        ! data from sfc to top
        ksfc = 1
        k1   = 1
      endif                       ! end_if_iflip
!
      do k = 1, LEVS
        do i = 1, IM
          totcnd = qcwat(i,k) + qcice(i,k) + qrain(i,k)
          qsnow = 0.0
          if(totcnd > EPSQ) then

!  ---  air density (rho), model mass thickness (cpath), temperature in c (tc)

            rho   = 0.1 * plyr(i,k)                                     &
     &            / (RD* tlyr(i,k) * (1.0 + EPS1*qlyr(i,k)))
            cpath = abs(plvl(i,k+1) - plvl(i,k)) * (100000.0 / GRAV)
            tc    = tlyr(i,k) - T0C

!! cloud water
!
!  ---  effective radius (recwat) & total water path (cwatp):
!       assume monodisperse distribution of droplets (no factor of 1.5)

            if (qcwat(i,k) > 0.0) then
              recwat(i,k) = max(RECWmin,recw1*(rho*qcwat(i,k))**CEXP)
              cwatp (i,k) = cpath * qcwat(i,k)           ! cloud water path
            endif

!! rain
!
!  ---  effective radius (rerain) & total water path (rainp):
!       factor of 1.5 accounts for r**3/r**2 moments for exponentially
!       distributed drops in effective radius calculations
!       (from m.d. chou's code provided to y.-t. hou)

            if (qrain(i,k) > 0.0) then
              tem         = CN0r0 * sqrt(sqrt(rho*qrain(i,k)))
              rerain(i,k) = 1.5 * max(XMRmin, min(XMRmax, tem))
              rainp (i,k) = cpath * qrain(i,k)           ! rain water path
            endif

!! snow (large ice) & cloud ice
!
!  ---  effective radius (resnow) & total ice path (snowp) for snow, and
!       total ice path (cicep) for cloud ice:
!       factor of 1.5 accounts for r**3/r**2 moments for exponentially
!       distributed ice particles in effective radius calculations
!       separation of cloud ice & "snow" uses algorithm from subroutine gsmcolumn

            pfac = 1.0

            if (qcice(i,k) > 0.0) then

!  ---  mean particle size following houze et al. (jas, 1979, p. 160),
!       converted from fig. 5 plot of lamdas.  an analogous set of
!       relationships also shown by fig. 8 of ryan (bams, 1996, p. 66),
!       but with a variety of different relationships that parallel
!       the houze curves.

              dum = max(0.05, min(1.0, exp(0.0564*tc) ))
              indexs = min(MDImax, max(MDImin, int(XMImax*dum) ))
              DUM=MAX(FLGmin(i)*pfac, DUM)

!  ---  assumed number fraction of large ice to total (large & small) ice
!       particles, which is based on a general impression of the literature.
!       small ice are assumed to have a mean diameter of 50 microns.

              if (tc >= 0.0) then
                flarge = FLG1P0
              else
                flarge = dum
              endif
!------------------------commented by moorthi -----------------------------
!             elseif (tc >= -25.0) then
!
!  ---  note that absence of cloud water (qcwat) is used as a quick
!       substitute for calculating water subsaturation as in gsmcolumn
!
!               if (qcwat(i,k) <= 0.0 .or. tc < -8.0                 &
!    &                                .or. tc > -3.0) then
!                 flarge = FLG0P2
!               else
!
!  ---  parameterize effects of rime splintering by increasing
!       number of small ice particles
!
!                 flarge = FLG0P1
!               endif
!             elseif (tc <= -50.0) then
!               flarge = 0.01
!             else
!               flarge = 0.2 * exp(0.1198*(tc+25.0))
!             endif
!____________________________________________________________________________

              xsimass = MASSI(MDImin) * (1.0 - flarge) / flarge
              NLImax=10.E3/sqrt(DUM)       !- Ver3

              tem = rho * qcice(i,k)
              nlice = tem / (xsimass +rrime(i,k)*MASSI(indexs))

!  ---  from subroutine gsmcolumn:
!       minimum number concentration for large ice of NLImin=10/m**3
!       at t>=0c.  done in order to prevent unrealistically small
!       melting rates and tiny amounts of snow from falling to
!       unrealistically warm temperatures.

              if (tc >= 0.0) then

                nlice = max(NLImin, nlice)

              elseif (nlice > nlimax) then

!  ---  ferrier 6/13/01:  prevent excess accumulation of ice

                xli = (tem/nlimax - xsimass) / rrime(i,k)

                if (xli <= MASSI(450) ) then
                  dsnow = 9.5885e5 * xli**0.42066
                else
                  dsnow = 3.9751e6 * xli** 0.49870
                endif

                indexs = min(MDImax, max(indexs, int(dsnow)))
                nlice = tem / (xsimass + rrime(i,k)*MASSI(indexs))

              endif                               ! end if_tc block

              if (plvl(i,ksfc) > 850.0 .and.                            &
     &            plvl(i,k+k1) > 700.0 .and. indexs >= indexsmin) then ! 20060516
                qsnow = min( qcice(i,k),                                &
     &                       nlice*rrime(i,k)*MASSI(indexs)/rho )
              endif

              qqcice      = max(0.0, qcice(i,k)-qsnow)
              cicep (i,k) = cpath * qqcice          ! cloud ice path
              resnow(i,k) = 1.5 * float(indexs)
              snden (i,k) = SDENS(indexs) / rrime(i,k)   ! 1/snow density
              snowp (i,k) = cpath*qsnow             ! snow path

            endif                                 ! end if_qcice block
          endif                                   ! end if_totcnd block

        enddo
      enddo
!
!...................................
      end subroutine n_rsipath2
!-----------------------------------

      end MODULE module_n_microphysics_gfs
