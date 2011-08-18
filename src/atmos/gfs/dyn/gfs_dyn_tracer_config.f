!
!! ! Module: gfs_dyn_tracer_config
!
! ! Description: gfs dynamics tracer configuration module
!
! ! Revision history:
!   Aug 08 2009   Sarah Lu, initial code
!   Aug 09 2009   Sarah Lu, add ntrac_chem, ntrac_met
!   Aug 10 2009   Sarah Lu, gfs_dyn_tracer is determined from ChemRegistry
!   Oct 16 2009   Sarah Lu, remove ChemRegistry; hardwire tracer specification 
!                           for testing; port to the latest trunk
!   Nov 13, 2009  Weiyu Yang, modified for the ensemble GEFS code.
!   Nov 19 2009   Sarah Lu, chem tracer specified from ChemRegistry
!   Feb 09 2009   Sarah Lu, ri/cpi added to gfs_dyn_tracer_type
!   Aug 17 2010   Sarah Lu, remove debug print
!   Aug 30 2010   Sarah Lu, set glbsum default as F
!   Aug 08 2011   Jun Wang, compile gocart only when running GOCART
! -------------------------------------------------------------------------
!
      module gfs_dyn_tracer_config
      use gfs_dyn_machine , only : kind_grid
      use gfs_dyn_tracer_const, only : cpi,ri

      implicit none
      SAVE
!
! tracer specification
!
      type    gfs_dyn_tracer_type
        character*20,         pointer      :: vname(:, :)    ! variable name
        real(kind=kind_grid), pointer      :: ri(:)
        real(kind=kind_grid), pointer      :: cpi(:)
        integer                  :: ntrac
        integer                  :: ntrac_met
        integer                  :: ntrac_chem
        logical                  :: doing_DU
        logical                  :: doing_SU
        logical                  :: doing_SS
        logical                  :: doing_OC
        logical                  :: doing_BC
        logical                  :: doing_GOCART
      endtype gfs_dyn_tracer_type

      type (gfs_dyn_tracer_type), save     ::  gfs_dyn_tracer
!
! misc tracer options
!
      logical, save                  :: glbsum  = .false.
!

! --- public interface
      public     tracer_config_init

      contains

! -------------------------------------------------------------------   
      subroutine tracer_config_init (gfs_dyn_tracer,ntrac,
     &                               ntoz,ntcw,ncld,me)

c  
c  This subprogram sets up gfs_dyn_tracer
c 
      implicit none
! input
      integer, intent(in)    ::  me, ntoz,ntcw,ncld
! output
      type (gfs_dyn_tracer_type), intent(out)    ::  gfs_dyn_tracer
!
      integer, intent(inout)  :: ntrac
! local
      integer                 :: i, status, ierr

! ntrac_chem = number of chem tracers
      gfs_dyn_tracer%ntrac_chem = 0
      gfs_dyn_tracer%doing_OC = .false.
      gfs_dyn_tracer%doing_BC = .false.
      gfs_dyn_tracer%doing_DU = .false.
      gfs_dyn_tracer%doing_SS = .false.
      gfs_dyn_tracer%doing_SU = .false.
      gfs_dyn_tracer%doing_GOCART = .false.

! ntrac_met = number of met tracers
      if ( ntoz < ntcw ) then                       
        gfs_dyn_tracer%ntrac_met = ntcw + ncld - 1   
      else                                                           
        gfs_dyn_tracer%ntrac_met = ntoz                              
      endif                                          
      if ( gfs_dyn_tracer%ntrac_met /= ntrac ) then
        print *,'LU_TRC: ERROR ! inconsistency in ntrac:',
     &           ntrac, gfs_dyn_tracer%ntrac_met
        stop     
      endif

! update ntrac = total number of tracers
      gfs_dyn_tracer%ntrac = gfs_dyn_tracer%ntrac_met +     
     &                       gfs_dyn_tracer%ntrac_chem
      ntrac = gfs_dyn_tracer%ntrac

! Set up tracer name, cpi, and ri
      if ( gfs_dyn_tracer%ntrac > 0 ) then      
       allocate(gfs_dyn_tracer%vname(ntrac, 5), stat=status)
           if( status .ne. 0 ) go to 999         
       allocate(gfs_dyn_tracer%ri(0:ntrac),  stat=status)
           if( status .ne. 0 ) go to 999
       allocate(gfs_dyn_tracer%cpi(0:ntrac), stat=status)
           if( status .ne. 0 ) go to 999
!--- fill in met tracers
      gfs_dyn_tracer%vname(1,    1) = 'spfh'   
      gfs_dyn_tracer%vname(1,    2) = 'spfh_q'   
      gfs_dyn_tracer%vname(1,    3) = 'spfh_m'   
      gfs_dyn_tracer%vname(1,    4) = 'spfh_q6'   
      gfs_dyn_tracer%vname(1,    5) = 'spfh_m6'   
      if(ntoz>1) then
        gfs_dyn_tracer%vname(ntoz, 1) = 'o3mr'  
        gfs_dyn_tracer%vname(ntoz, 2) = 'o3mr_q'  
        gfs_dyn_tracer%vname(ntoz, 3) = 'o3mr_m'  
        gfs_dyn_tracer%vname(ntoz, 4) = 'o3mr_q6'  
        gfs_dyn_tracer%vname(ntoz, 5) = 'o3mr_m6'  
      endif
      if(ntcw>1) then
        gfs_dyn_tracer%vname(ntcw, 1) = 'clwmr' 
        gfs_dyn_tracer%vname(ntcw, 2) = 'clwmr_q' 
        gfs_dyn_tracer%vname(ntcw, 3) = 'clwmr_m' 
        gfs_dyn_tracer%vname(ntcw, 4) = 'clwmr_q6' 
        gfs_dyn_tracer%vname(ntcw, 5) = 'clwmr_m6' 
      endif

      gfs_dyn_tracer%cpi(0:gfs_dyn_tracer%ntrac_met) =
     &               cpi(0:gfs_dyn_tracer%ntrac_met)
      gfs_dyn_tracer%ri(0:gfs_dyn_tracer%ntrac_met) =
     &               ri(0:gfs_dyn_tracer%ntrac_met)

      endif
!      write(0,*)'in trac_config,vname=',
!     &  gfs_dyn_tracer%vname(1:3,1),'ntoz=',ntoz,'ntcw=',ntcw
!
!-- call chem tracer if gocart is running
      call dyn_gocart_tracer_config(gfs_dyn_tracer,ntrac,
     &                               ntoz,ntcw,ncld,me)
!      write(0,*)'in trac_config,af gocart,vname=',
!     &  gfs_dyn_tracer%vname(1:3,1),'ntoz=',ntoz,'ntcw=',ntcw

      return

999   print *,'LU_TRC: error in allocate gfs_dyn_tracer :',status,me

      end subroutine tracer_config_init

! ========================================================================= 

      end module gfs_dyn_tracer_config
