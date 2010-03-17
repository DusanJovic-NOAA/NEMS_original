!-------------------------------------------------------------------------------
  subroutine gfuncphys
!$$$     Subprogram Documentation Block
!
! Subprogram: gfuncphys    Compute all physics function tables
!   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
!
! Abstract: Compute all physics function tables.  Lookup tables are
!   set up for computing saturation vapor pressure, dewpoint temperature,
!   equivalent potential temperature, moist adiabatic temperature and humidity,
!   pressure to the kappa, and lifting condensation level temperature.
!
! Program History Log:
! 1999-03-01  Iredell             f90 module
!
! Usage:  call gfuncphys
!
! Subprograms called:
!   gpvsl       compute saturation vapor pressure over liquid table
!   gpvsi       compute saturation vapor pressure over ice table
!   gpvs        compute saturation vapor pressure table
!   gtdpl       compute dewpoint temperature over liquid table
!   gtdpi       compute dewpoint temperature over ice table
!   gtdp        compute dewpoint temperature table
!   gthe        compute equivalent potential temperature table
!   gtma        compute moist adiabat tables
!   gpkap       compute pressure to the kappa table
!   grkap       compute pressure to the 1/kappa table
!   gtlcl       compute LCL temperature table
!
! Attributes:
!   Language: Fortran 90.
!
!$$$
    implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine
!-------------------------------------------------------------------------------

