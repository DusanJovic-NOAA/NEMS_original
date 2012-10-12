!-----------------------------------------------------------------------
      subroutine compns_dynamics (deltim,iret,
     &  ntrac,nxpt,nypt,jintmx,jcap,jcapg,
     &  levs,levr,lonf,latg, ntoz,
     &  ntcw,ncld, spectral_loop, me,
     &  thermodyn_id,sfcpress_id,
     &  nlunit, gfs_dyn_namelist,ndfi)
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
!             9. the difgital filter  must be non-negative and
!
! Program History Log:
!   1999-01-26  Iredell
!   2007-02-01  H.-M. H. Juang modify to be used for dynamics only
!   2009-11-09  Jun Wang       added ndfi 
!   2010-09-08  Jun Wang       change gfsio to nemsio
!   2011-02-11  Henry Juang    add codes to fit mass_dp and ndslfv
!   2011-02-28  Sarah Lu       add thermodyn_id and sfcpress_id
!   2011-03-15  Henry Juang    add jcapg for usual truncation in grid
!                              for Eulerian jcap=jcapg<lonf/3
!                              for NDSL jcap<lonf/2, jcapg<lonf/3
!   2012-04-06  Henry Juang    add idea for lsidea
!   2012-10-05  Jun Wang       add sigio_out 
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
      integer ntrac,nxpt,nypt,jintmx,jcap,jcapg,levs,lonf,latg
      integer levr,lsoil,nmtvr,lonr,latr,ndfi
      integer ntoz,ntcw,ncld,spectral_loop,member_num
      integer thermodyn_id, sfcpress_id
      real    tfiltc

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      namelist /nam_dyn/FHMAX,FHOUT,FHRES,FHROT,FHDFI,DELTIM,IGEN,
     & NGPTC,shuff_lats_a,reshuff_lats_a,
     & nxpt,nypt,jintmx,jcap,jcapg,levs,lonf,latg,levr,
     & ntrac,ntoz,ntcw,ncld,nsout,tfiltc,
     & nemsio_in,nemsio_out,liope,ref_temp,lsidea,
     & explicit,hybrid,gen_coord_hybrid,process_split, 
     & spectral_loop,ndslfv,mass_dp,semi_implicit_temp_profile,
     & reduced_grid,
     & thermodyn_id, sfcpress_id,
     & zflxtvd,num_reduce,sigio_out

!
      num_reduce=-4
      fhmax=0
      fhout=0
      fhres=0
      fhrot=0
      fhdfi=0
      deltim=0
      igen=0
      tfiltc  = 0.85
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
      mass_dp          = .false.                                     !hmhj
      semi_implicit_temp_profile    = .false.                        !hmhj
      process_split    = .false.                                     !hmhj
      liope            = .true.
!
      thermodyn_id     = 1
      sfcpress_id      = 1
!
      zflxtvd          = .false.
      ndslfv           = .false. ! non_iteration semi_Lagrangian finite volume
!
      nemsio_in         = .true.
      nemsio_out        = .true.
      sigio_out         = .false.
! idea add
      lsidea           = .false. 
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

! idea add
      if( lsidea ) then
        if (levs > 100 .and. ref_temp < 400.0) ref_temp = 2500.0
      endif

      if (me.eq.0) write(6,nam_dyn)
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
!  Compute ndfi and check rule 7.
      if(fhdfi.eq.0.) then
        ndfi=0
      else
        ndfi=nint(2*fhdfi*3600./deltim)
        if(ndfi.le.0.or.abs(ndfi-2*fhdfi*3600./deltim).gt.tol.or.
     &     ndfi.gt.nsres) then
           print *,'ndfi=',ndfi,'is not equal to',2*fhdfi*3600./deltim
          iret=7
          return
        endif
      endif
!
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
      return
      end subroutine compns_dynamics
