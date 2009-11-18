!
!! ! Module: gfs_phy_tracer_config
!
! ! Description: gfs physics tracer configuration module
!
! ! Revision history:
!   Oct 16 2009   Sarah Lu, adopted from dyn fc
! -------------------------------------------------------------------------
!
      module gfs_phy_tracer_config
      use machine , only : kind_phys
      implicit none
      SAVE
!
! tracer specification
!
      type    gfs_phy_tracer_type
        character*20,    pointer     :: vname(:)    ! variable name
        integer                  :: ntrac
        integer                  :: ntrac_met
        integer                  :: ntrac_chem
      endtype gfs_phy_tracer_type

      type (gfs_phy_tracer_type), save     ::  gfs_phy_tracer
!
! misc tracer options
!
      logical, save                  :: glbsum  = .true.
!

! --- public interface
      public     tracer_config_init

      contains

! -------------------------------------------------------------------   
      subroutine tracer_config_init (gfs_phy_tracer,ntrac,
     &                               ntoz,ntcw,ncld,me)

c  
c  This subprogram sets up gfs_phy_tracer
c 
      implicit none
! input
      integer, intent(in)    ::  me, ntoz,ntcw,ncld
! output
      type (gfs_phy_tracer_type), intent(out)    ::  gfs_phy_tracer
!
      integer, intent(inout)  :: ntrac
! local
      integer                 :: i, status, ierr

! ntrac_chem = number of chem tracers

      gfs_phy_tracer%ntrac_chem = 0                        
! ntrac_met = number of met tracers
      if ( ntoz < ntcw ) then                       
        gfs_phy_tracer%ntrac_met = ntcw + ncld - 1   
      else                                                           
        gfs_phy_tracer%ntrac_met = ntoz                              
      endif                                          
      if ( gfs_phy_tracer%ntrac_met /= ntrac ) then
        print *,'LU_TRC: ERROR ! inconsistency in ntrac:',
     &           ntrac, gfs_phy_tracer%ntrac_met
        stop     
      endif

! update ntrac = total number of tracers
      gfs_phy_tracer%ntrac = gfs_phy_tracer%ntrac_met +     
     &                       gfs_phy_tracer%ntrac_chem
      ntrac = gfs_phy_tracer%ntrac

! Set up tracer name
      if ( gfs_phy_tracer%ntrac > 0 ) then      
       allocate(gfs_phy_tracer%vname(ntrac), stat=status)
           if( status .ne. 0 ) go to 999         
      gfs_phy_tracer%vname(1) = 'spfh'   
      gfs_phy_tracer%vname(ntoz) = 'o3mr'  
      gfs_phy_tracer%vname(ntcw) = 'clwmr' 
      endif

      print *,'LU_TRC: exit tracer_config_init'
      return

999   print *,'LU_TRC: error in allocate gfs_phy_tracer :',status,me

      end subroutine tracer_config_init

! ========================================================================= 

      end module gfs_phy_tracer_config