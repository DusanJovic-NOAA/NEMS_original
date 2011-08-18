!
!! ! Subroutine : gocart_tracer_config
!
! ! Description: reset gfs phys when gocart is running
!
! ! Revision history:
!   Aug 09 2011   Jun Wang, initial code from tracer_config_init
! -------------------------------------------------------------------------
!
      subroutine gocart_tracer_config(gfs_phy_tracers,ntrac,
     &                                     ntoz,ntcw,ncld,me)
!
      use gfs_phy_tracer_config
      use Chem_RegistryMod
!
      implicit none
!
! input
      integer, intent(in)    ::  me, ntoz,ntcw,ncld
! output
      type (gfs_phy_tracer_type), intent(out)    ::  gfs_phy_tracers
!
      integer, intent(inout)  :: ntrac
! local
      integer                 :: i, status, ierr
      type(Chem_Registry)     :: reg


! Read Chem_Registry
      reg = Chem_RegistryCreate ( ierr )

!!    if ( me == 0) call Chem_RegistryPrint (reg)

! ntrac_chem = number of chem tracers
      gfs_phy_tracers%ntrac_chem = reg%nq
      gfs_phy_tracers%doing_OC = reg%doing_OC
      gfs_phy_tracers%doing_BC = reg%doing_BC
      gfs_phy_tracers%doing_DU = reg%doing_DU
      gfs_phy_tracers%doing_SS = reg%doing_SS
      gfs_phy_tracers%doing_SU = reg%doing_SU
      gfs_phy_tracers%doing_GOCART = reg%doing_GOCART

! ntrac_met = number of met tracers
      if ( ntoz < ntcw ) then
        gfs_phy_tracers%ntrac_met = ntcw + ncld - 1
      else
        gfs_phy_tracers%ntrac_met = ntoz
      endif
      if ( gfs_phy_tracers%ntrac_met /= ntrac ) then
        print *,'TRAC_CONFIG: ERROR ! inconsistency in ntrac:',
     &                 ntrac, gfs_phy_tracers%ntrac_met
        stop
      endif
!
! update ntrac = total number of tracers
      gfs_phy_tracers%ntrac = gfs_phy_tracers%ntrac_met +       
     &                       gfs_phy_tracers%ntrac_chem
      ntrac = gfs_phy_tracers%ntrac

! Set up tracer name and allocate ri, cpi, fscav
      if(associated(gfs_phy_tracers%vname))
     &   deallocate(gfs_phy_tracers%vname)
      if(associated(gfs_phy_tracers%ri))deallocate(gfs_phy_tracers%ri)
      if(associated(gfs_phy_tracers%cpi))deallocate(gfs_phy_tracers%cpi)
      if(associated(gfs_phy_tracers%fscav))
     &   deallocate(gfs_phy_tracers%fscav)
      if ( gfs_phy_tracers%ntrac > 0 ) then      
       allocate(gfs_phy_tracers%vname(ntrac), stat=status)
           if( status .ne. 0 ) go to 999         
       allocate(gfs_phy_tracers%ri(0:ntrac), stat=status)
           if( status .ne. 0 ) go to 999         
       allocate(gfs_phy_tracers%cpi(0:ntrac), stat=status)
           if( status .ne. 0 ) go to 999         
       allocate(gfs_phy_tracers%fscav(ntrac), stat=status)
           if( status .ne. 0 ) go to 999         

!--- fill in met tracers
      gfs_phy_tracers%vname(1) = 'spfh'   
      gfs_phy_tracers%vname(ntoz) = 'o3mr'  
      gfs_phy_tracers%vname(ntcw) = 'clwmr' 
!--- fill in chem tracers
      do i = 1,gfs_phy_tracers%ntrac_chem
       gfs_phy_tracers%vname(i+gfs_phy_tracers%ntrac_met)=reg%vname(i)
      enddo
!--- fill in default values for fscav
      gfs_phy_tracers%fscav(:) = 0.

      endif

! Destroy Chem_Registry
      call Chem_RegistryDestroy ( reg, ierr )

      return

999   print *,'TRAC_CONFIG: error in allocate gfs_phy_tracers :',
     &    status,me

      end subroutine gocart_tracer_config
