!
!! ! Module: gfs_phy_tracer_config
!
! ! Description: gfs physics tracer configuration module
!
! ! Revision history:
!   Oct 16 2009   Sarah Lu, adopted from dyn fc
!   Nov 21 2009   Sarah Lu, chem tracer specified from ChemRegistry
!   Dec 10 2009   Sarah Lu, add doing_GOCART
! -------------------------------------------------------------------------
!
      module gfs_phy_tracer_config
      use machine , only : kind_phys
      use Chem_RegistryMod
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
        logical                  :: doing_DU
        logical                  :: doing_SU
        logical                  :: doing_SS
        logical                  :: doing_OC
        logical                  :: doing_BC
        logical                  :: doing_GOCART
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
! -------------------------------------------------------------------   
       subroutine tracer_config_init (gfs_phy_tracer,ntrac,          
     &                                     ntoz,ntcw,ncld,me)

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
      type(Chem_Registry)     :: reg


! Read Chem_Registry
      reg = Chem_RegistryCreate ( ierr )
      call Chem_RegistryPrint (reg)

! ntrac_chem = number of chem tracers
      gfs_phy_tracer%ntrac_chem = reg%nq
      gfs_phy_tracer%doing_OC = reg%doing_OC
      gfs_phy_tracer%doing_BC = reg%doing_BC
      gfs_phy_tracer%doing_DU = reg%doing_DU
      gfs_phy_tracer%doing_SS = reg%doing_SS
      gfs_phy_tracer%doing_SU = reg%doing_SU
      gfs_phy_tracer%doing_GOCART = reg%doing_GOCART

! ntrac_met = number of met tracers
      if ( ntoz < ntcw ) then                       
        gfs_phy_tracer%ntrac_met = ntcw + ncld - 1   
      else                                                           
        gfs_phy_tracer%ntrac_met = ntoz                              
      endif                                          
      if ( gfs_phy_tracer%ntrac_met /= ntrac ) then
        print *,'LU_TRC: ERROR ! inconsistency in ntrac:',    
     &                 ntrac, gfs_phy_tracer%ntrac_met
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
!--- fill in met tracers
      gfs_phy_tracer%vname(1) = 'spfh'   
      gfs_phy_tracer%vname(ntoz) = 'o3mr'  
      gfs_phy_tracer%vname(ntcw) = 'clwmr' 
!--- fill in chem tracers
      do i = 1,gfs_phy_tracer%ntrac_chem
       gfs_phy_tracer%vname(i+gfs_phy_tracer%ntrac_met)=reg%vname(i)
      enddo

      endif

! Destroy Chem_Registry
      call Chem_RegistryDestroy ( reg, ierr )

      print *,'LU_TRC: exit tracer_config_init'

      return

999   print *,'LU_TRC: error in allocate gfs_phy_tracer :',status,me

      end subroutine tracer_config_init

! ========================================================================= 

      end module gfs_phy_tracer_config
