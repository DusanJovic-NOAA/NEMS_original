!-----------------------------------------------------------------------
      subroutine compns_dynamics (deltim,iret,
     &  ntrac,nxpt,nypt,jintmx,jcap,
     &  levs,levr,lonf,latg, ntoz,
     &  ntcw,ncld, spectral_loop, me,
     &  nlunit, gfs_dyn_namelist)
!$$$  Subprogram Documentation Block
!
! Subprogram:  compns     Check and compute namelist frequencies
!   Prgmmr: Iredell       Org: NP23          Date: 1999-01-26
!
! Abstract: This subprogram checks global spectral model namelist
!           frequencies in hour units for validity.  If they are valid,
!           then the frequencies are computed in timestep units.
!           The following rules are applied:
!             1. the timestep must be positive;
!             2. the output frequency must be positive and
!                a multiple of the timestep to within tolerance;
!             3. the shortwave frequency must be positive and
!                a multiple of the timestep to within tolerance;
!             4. the longwave frequency must be positive and
!                a multiple of the timestep to within tolerance and
!                a multiple of the shortwave frequency;
!             5. the zeroing frequency must be positive and
!                a multiple of the timestep to within tolerance and
!                a multiple of the output frequency;
!             6. the restart frequency must be positive and
!                a multiple of the timestep to within tolerance and
!                a multiple of the longwave frequency and
!                a multiple of the zeroing frequency;
!             7. the initialization window must be non-negative and
!                a multiple of the timestep to within tolerance and
!                a multiple of the longwave frequency and
!                no longer than the restart frequency;
!             8. the cycling frequency must be non-negative and
!                a multiple of the timestep to within tolerance and
!                a multiple of the longwave frequency.
!
! Program History Log:
!   1999-01-26  Iredell
!   2007-02-01  H.-M. H. Juang modify to be used for dynamics only
!
! Usage:    call compns(deltim,
!    &                  fhout,fhres,
!    &                  nsout,nsres,
!    &                  iret)
!   Input Arguments:
!     tol      - real error tolerance allowed for input frequencies
!                (e.g. 0.01 for 1% of timestep maximum error allowed)
!     deltim   - real timestep in seconds
!     fhout    - real output frequency in hours
!     fhres    - real restart frequency in hours
!   Output Arguments:
!     nsout    - integer output frequency in timesteps
!     nsres    - integer restart frequency in timesteps
!     iret     - integer return code (0 if successful or
!                between 1 and 8 for which rule above was broken)
!
! Attributes:
!   Language: Fortran 90
!
!$$$

      
      use namelist_dynamics_def
      use gfs_dyn_mpi_def, only : liope
      implicit none

      real tol
 
      character (len=*), intent(in) :: gfs_dyn_namelist
      integer, intent(in)           :: me, nlunit
      real,intent(inout)            :: deltim
      integer,intent(out)           :: iret
      integer ntrac,nxpt,nypt,jintmx,jcap,levs,lonf,latg
      integer levr,lsoil,nmtvr,lonr,latr
      integer ntoz,ntcw,ncld,spectral_loop,member_num
      real    tfiltc

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      namelist /nam_dyn/FHMAX,FHOUT,FHRES,FHROT,DELTIM,IGEN,
     & NGPTC,spectral_loop,nislfv,
     & shuff_lats_a,reshuff_lats_a,reduced_grid,
     & explicit,hybrid,gen_coord_hybrid,liope, 
     & ntrac,nxpt,nypt,jintmx,jcap,levs,lonf,latg,levr,
     & ntoz,ntcw,ncld,nsout,tfiltc,
     & gfsio_in,gfsio_out,ref_temp,
     & zflxtvd

!
      fhmax=0
      fhout=0
      fhres=0
      fhrot=0
      deltim=0
      igen=0
      tfiltc  = 0.92
      NGPTC=lonf
!
      shuff_lats_a     = .true.
      reshuff_lats_a   = .false.
!
      reduced_grid     = .true.
!
      explicit         = .false.
      hybrid           = .false.
      gen_coord_hybrid = .false.                                     !hmhj
      liope            = .true.
!
      zflxtvd          = .false.
      nislfv           = 0        ! non_iteration semi_Lagrangian finite volume
!
      gfsio_in         = .true.
      gfsio_out        = .true.
!
      ref_temp         = 300.0
!
      nsout            = 0
      levr             = 0
      spectral_loop    = 2	! 1 for one-loop or 2 for two-loop
!
      print *,' nlunit=',nlunit,' gfs_dyn_namelist=',gfs_dyn_namelist
c$$$      read(5,nam_dyn)
      open(unit=nlunit,file=gfs_dyn_namelist)
      rewind (nlunit)
      read(nlunit,nam_dyn)
      write(*,nam_dyn)
c
!     if (me.eq.0) write(6,nam_dyn)
      filta = tfiltc
!
      if (levr == 0) then
        levr = levs
      endif
c
csela - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tol=0.01
!  Check rule 1.
      if(deltim.le.0) then
        iret=1
        return
      endif
!      write(0,*)'lver=',levr,'deltim=',deltim,'nsout=',nsout,'fhout=',
!     & fhout,'fhres=',fhres,'gen_coord_hybrid=',gen_coord_hybrid, 
!     & 'ntoz=',ntoz
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsout and check rule 2.
      if(nsout.gt.0) fhout=nsout*deltim/3600.
      nsout=nint(fhout*3600./deltim)
      if(nsout.le.0.or.abs(nsout-fhout*3600./deltim).gt.tol) then
        iret=2
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsres and check rule 6.
      nsres=nint(fhres*3600./deltim)
      if(nsres.le.0.or.abs(nsres-fhres*3600./deltim).gt.tol) then
        iret=6
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      IF (NGPTC.GT.lonf) THEN
         NGPTC=lonf
         WRITE(0,*) "NGPTC IS TOO BIG, RESET NGPTC TO lonf",NGPTC
      ENDIF
      IF (ME.EQ.0)   WRITE(0,*) "NGPTC IS SET TO NGPTC :",NGPTC
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  All checks are successful.
!     print *,' done compns_dynamics '
      iret=0
c
      end subroutine compns_dynamics
