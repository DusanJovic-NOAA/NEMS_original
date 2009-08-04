!-------------------------------------------------------------------------------
      module funcphys
!$$$  Module Documentation Block
!
! Module:    funcphys        API for basic thermodynamic physics
!   Author: Iredell          Org: W/NX23     Date: 1999-03-01
!
! Abstract: This module provides an Application Program Interface
!   for computing basic thermodynamic physics functions, in particular
!     (1) saturation vapor pressure as a function of temperature,
!     (2) dewpoint temperature as a function of vapor pressure,
!     (3) equivalent potential temperature as a function of temperature
!         and scaled pressure to the kappa power,
!     (4) temperature and specific humidity along a moist adiabat
!         as functions of equivalent potential temperature and
!         scaled pressure to the kappa power,
!     (5) scaled pressure to the kappa power as a function of pressure, and
!     (6) temperature at the lifting condensation level as a function
!         of temperature and dewpoint depression.
!   The entry points required to set up lookup tables start with a "g".
!   All the other entry points are functions starting with an "f" or
!   are subroutines starting with an "s".  These other functions and
!   subroutines are elemental; that is, they return a scalar if they
!   are passed only scalars, but they return an array if they are passed
!   an array.  These other functions and subroutines can be inlined, too.
!
! Program History Log:
!   1999-03-01  Mark Iredell
!   1999-10-15  Mark Iredell  SI unit for pressure (Pascals)
!   2001-02-26  Mark Iredell  Ice phase changes of Hong and Moorthi
!
! Public Variables:
!   krealfp         Integer parameter kind or length of reals (=kind_phys)
!
! Public Subprograms:
!   gpvsl            Compute saturation vapor pressure over liquid table
!
!   fpvsl           Elementally compute saturation vapor pressure over liquid
!     function result Real(krealfp) saturation vapor pressure in Pascals
!     t               Real(krealfp) temperature in Kelvin
!
!   fpvslq          Elementally compute saturation vapor pressure over liquid
!     function result Real(krealfp) saturation vapor pressure in Pascals
!     t               Real(krealfp) temperature in Kelvin
!
!   fpvslx          Elementally compute saturation vapor pressure over liquid
!     function result Real(krealfp) saturation vapor pressure in Pascals
!     t               Real(krealfp) temperature in Kelvin
!
!   gpvsi            Compute saturation vapor pressure over ice table
!
!   fpvsi           Elementally compute saturation vapor pressure over ice
!     function result Real(krealfp) saturation vapor pressure in Pascals
!     t               Real(krealfp) temperature in Kelvin
!
!   fpvsiq          Elementally compute saturation vapor pressure over ice
!     function result Real(krealfp) saturation vapor pressure in Pascals
!     t               Real(krealfp) temperature in Kelvin
!
!   fpvsix          Elementally compute saturation vapor pressure over ice
!     function result Real(krealfp) saturation vapor pressure in Pascals
!     t               Real(krealfp) temperature in Kelvin
!
!   gpvs            Compute saturation vapor pressure table
!
!   fpvs            Elementally compute saturation vapor pressure
!     function result Real(krealfp) saturation vapor pressure in Pascals
!     t               Real(krealfp) temperature in Kelvin
!
!   fpvsq           Elementally compute saturation vapor pressure
!     function result Real(krealfp) saturation vapor pressure in Pascals
!     t               Real(krealfp) temperature in Kelvin
!
!   fpvsx           Elementally compute saturation vapor pressure
!     function result Real(krealfp) saturation vapor pressure in Pascals
!     t               Real(krealfp) temperature in Kelvin
!
!   gtdpl           Compute dewpoint temperature over liquid table
!
!   ftdpl           Elementally compute dewpoint temperature over liquid
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   ftdplq          Elementally compute dewpoint temperature over liquid
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   ftdplx          Elementally compute dewpoint temperature over liquid
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   ftdplxg         Elementally compute dewpoint temperature over liquid
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     t               Real(krealfp) guess dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   gtdpi           Compute dewpoint temperature table over ice
!
!   ftdpi           Elementally compute dewpoint temperature over ice
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   ftdpiq          Elementally compute dewpoint temperature over ice
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   ftdpix          Elementally compute dewpoint temperature over ice
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   ftdpixg         Elementally compute dewpoint temperature over ice
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     t               Real(krealfp) guess dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   gtdp            Compute dewpoint temperature table
!
!   ftdp            Elementally compute dewpoint temperature
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   ftdpq           Elementally compute dewpoint temperature
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   ftdpx           Elementally compute dewpoint temperature
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   ftdpxg          Elementally compute dewpoint temperature
!     function result Real(krealfp) dewpoint temperature in Kelvin
!     t               Real(krealfp) guess dewpoint temperature in Kelvin
!     pv              Real(krealfp) vapor pressure in Pascals
!
!   gthe            Compute equivalent potential temperature table
!
!   fthe            Elementally compute equivalent potential temperature
!     function result Real(krealfp) equivalent potential temperature in Kelvin
!     t               Real(krealfp) LCL temperature in Kelvin
!     pk              Real(krealfp) LCL pressure over 1e5 Pa to the kappa power
!
!   ftheq           Elementally compute equivalent potential temperature
!     function result Real(krealfp) equivalent potential temperature in Kelvin
!     t               Real(krealfp) LCL temperature in Kelvin
!     pk              Real(krealfp) LCL pressure over 1e5 Pa to the kappa power
!
!   fthex           Elementally compute equivalent potential temperature
!     function result Real(krealfp) equivalent potential temperature in Kelvin
!     t               Real(krealfp) LCL temperature in Kelvin
!     pk              Real(krealfp) LCL pressure over 1e5 Pa to the kappa power
!
!   gtma            Compute moist adiabat tables
!
!   stma            Elementally compute moist adiabat temperature and moisture
!     the             Real(krealfp) equivalent potential temperature in Kelvin
!     pk              Real(krealfp) pressure over 1e5 Pa to the kappa power
!     tma             Real(krealfp) parcel temperature in Kelvin
!     qma             Real(krealfp) parcel specific humidity in kg/kg
!
!   stmaq           Elementally compute moist adiabat temperature and moisture
!     the             Real(krealfp) equivalent potential temperature in Kelvin
!     pk              Real(krealfp) pressure over 1e5 Pa to the kappa power
!     tma             Real(krealfp) parcel temperature in Kelvin
!     qma             Real(krealfp) parcel specific humidity in kg/kg
!
!   stmax           Elementally compute moist adiabat temperature and moisture
!     the             Real(krealfp) equivalent potential temperature in Kelvin
!     pk              Real(krealfp) pressure over 1e5 Pa to the kappa power
!     tma             Real(krealfp) parcel temperature in Kelvin
!     qma             Real(krealfp) parcel specific humidity in kg/kg
!
!   stmaxg          Elementally compute moist adiabat temperature and moisture
!     tg              Real(krealfp) guess parcel temperature in Kelvin
!     the             Real(krealfp) equivalent potential temperature in Kelvin
!     pk              Real(krealfp) pressure over 1e5 Pa to the kappa power
!     tma             Real(krealfp) parcel temperature in Kelvin
!     qma             Real(krealfp) parcel specific humidity in kg/kg
!
!   gpkap           Compute pressure to the kappa table
!
!   fpkap           Elementally raise pressure to the kappa power.
!     function result Real(krealfp) p over 1e5 Pa to the kappa power
!     p               Real(krealfp) pressure in Pascals
!
!   fpkapq          Elementally raise pressure to the kappa power.
!     function result Real(krealfp) p over 1e5 Pa to the kappa power
!     p               Real(krealfp) pressure in Pascals
!
!   fpkapo          Elementally raise pressure to the kappa power.
!     function result Real(krealfp) p over 1e5 Pa to the kappa power
!     p               Real(krealfp) surface pressure in Pascals
!
!   fpkapx          Elementally raise pressure to the kappa power.
!     function result Real(krealfp) p over 1e5 Pa to the kappa power
!     p               Real(krealfp) pressure in Pascals
!
!   grkap           Compute pressure to the 1/kappa table
!
!   frkap           Elementally raise pressure to the 1/kappa power.
!     function result Real(krealfp) pressure in Pascals
!     pkap            Real(krealfp) p over 1e5 Pa to the 1/kappa power
!
!   frkapq          Elementally raise pressure to the kappa power.
!     function result Real(krealfp) pressure in Pascals
!     pkap            Real(krealfp) p over 1e5 Pa to the kappa power
!
!   frkapx          Elementally raise pressure to the kappa power.
!     function result Real(krealfp) pressure in Pascals
!     pkap            Real(krealfp) p over 1e5 Pa to the kappa power
!
!   gtlcl           Compute LCL temperature table
!
!   ftlcl           Elementally compute LCL temperature.
!     function result Real(krealfp) temperature at the LCL in Kelvin
!     t               Real(krealfp) temperature in Kelvin
!     tdpd            Real(krealfp) dewpoint depression in Kelvin
!
!   ftlclq          Elementally compute LCL temperature.
!     function result Real(krealfp) temperature at the LCL in Kelvin
!     t               Real(krealfp) temperature in Kelvin
!     tdpd            Real(krealfp) dewpoint depression in Kelvin
!
!   ftlclo          Elementally compute LCL temperature.
!     function result Real(krealfp) temperature at the LCL in Kelvin
!     t               Real(krealfp) temperature in Kelvin
!     tdpd            Real(krealfp) dewpoint depression in Kelvin
!
!   ftlclx          Elementally compute LCL temperature.
!     function result Real(krealfp) temperature at the LCL in Kelvin
!     t               Real(krealfp) temperature in Kelvin
!     tdpd            Real(krealfp) dewpoint depression in Kelvin
!
!   gfuncphys       Compute all physics function tables
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  implicit none
!-------------------------------------------------------------------------------
  end module funcphys
