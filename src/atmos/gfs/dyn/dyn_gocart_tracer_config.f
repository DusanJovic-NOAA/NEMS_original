!
!
!! ! Subroutine : dyn_gocart_tracer_config
!
! ! Description: reset gfs dyn when gocart is running
!
! ! Revision history:
!   Aug 09 2011   Jun Wang, initial code from tracer_config_init
! -------------------------------------------------------------------------
!
      subroutine dyn_gocart_tracer_config (gfs_dyn_tracers,ntrac,
     &                               ntoz,ntcw,ncld,me)
!
      use gfs_dyn_tracer_config
      use Chem_RegistryMod
!
      implicit none
! input
      integer, intent(in)    ::  me, ntoz,ntcw,ncld
! output
      type (gfs_dyn_tracer_type), intent(out)    ::  gfs_dyn_tracers
!
      integer, intent(inout)  :: ntrac
! local
      integer                 :: i, status, ierr
      type(Chem_Registry)     :: reg

! Read Chem_Registry
      reg = Chem_RegistryCreate ( ierr )

! ntrac_chem = number of chem tracers
      gfs_dyn_tracers%ntrac_chem = reg%nq
      gfs_dyn_tracers%doing_OC = reg%doing_OC
      gfs_dyn_tracers%doing_BC = reg%doing_BC
      gfs_dyn_tracers%doing_DU = reg%doing_DU
      gfs_dyn_tracers%doing_SS = reg%doing_SS
      gfs_dyn_tracers%doing_SU = reg%doing_SU
      gfs_dyn_tracers%doing_GOCART = reg%doing_GOCART

! ntrac_met = number of met tracers
      if ( ntoz < ntcw ) then                       
        gfs_dyn_tracers%ntrac_met = ntcw + ncld - 1   
      else                                                           
        gfs_dyn_tracers%ntrac_met = ntoz                              
      endif                                          
      if ( gfs_dyn_tracers%ntrac_met /= ntrac ) then
        print *,'LU_TRC: ERROR ! inconsistency in ntrac:',
     &           ntrac, gfs_dyn_tracers%ntrac_met
        stop     
      endif

! update ntrac = total number of tracers
      gfs_dyn_tracers%ntrac = gfs_dyn_tracers%ntrac_met +     
     &                       gfs_dyn_tracers%ntrac_chem
      ntrac = gfs_dyn_tracers%ntrac

! Set up tracer name, cpi, and ri
      if(associated(gfs_dyn_tracers%vname)) 
     &  deallocate(gfs_dyn_tracers%vname)
      if(associated(gfs_dyn_tracers%ri)) 
     &  deallocate(gfs_dyn_tracers%ri)
      if(associated(gfs_dyn_tracers%cpi)) 
     &  deallocate(gfs_dyn_tracers%cpi)
      if ( gfs_dyn_tracers%ntrac > 0 ) then      
       allocate(gfs_dyn_tracers%vname(ntrac, 5), stat=status)
           if( status .ne. 0 ) go to 999         
       allocate(gfs_dyn_tracers%ri(0:ntrac),  stat=status)
           if( status .ne. 0 ) go to 999
       allocate(gfs_dyn_tracers%cpi(0:ntrac), stat=status)
           if( status .ne. 0 ) go to 999
!--- fill in met tracers
      gfs_dyn_tracers%vname(1,    1) = 'spfh'   
      gfs_dyn_tracers%vname(1,    2) = 'spfh_q'   
      gfs_dyn_tracers%vname(1,    3) = 'spfh_m'   
      gfs_dyn_tracers%vname(1,    4) = 'spfh_q6'   
      gfs_dyn_tracers%vname(1,    5) = 'spfh_m6'   
      if(ntoz>1) then
        gfs_dyn_tracers%vname(ntoz, 1) = 'o3mr'  
        gfs_dyn_tracers%vname(ntoz, 2) = 'o3mr_q'  
        gfs_dyn_tracers%vname(ntoz, 3) = 'o3mr_m'  
        gfs_dyn_tracers%vname(ntoz, 4) = 'o3mr_q6'  
        gfs_dyn_tracers%vname(ntoz, 5) = 'o3mr_m6'  
      endif
      if(ntcw>1) then
        gfs_dyn_tracers%vname(ntcw, 1) = 'clwmr' 
        gfs_dyn_tracers%vname(ntcw, 2) = 'clwmr_q' 
        gfs_dyn_tracers%vname(ntcw, 3) = 'clwmr_m' 
        gfs_dyn_tracers%vname(ntcw, 4) = 'clwmr_q6' 
        gfs_dyn_tracers%vname(ntcw, 5) = 'clwmr_m6' 
      endif

      gfs_dyn_tracers%cpi(0:gfs_dyn_tracers%ntrac_met) =
     &               cpi(0:gfs_dyn_tracers%ntrac_met)
      gfs_dyn_tracers%ri(0:gfs_dyn_tracers%ntrac_met) =
     &               ri(0:gfs_dyn_tracers%ntrac_met)

!--- fill in chem tracers
      do i = 1,gfs_dyn_tracers%ntrac_chem
       gfs_dyn_tracers%vname(i+gfs_dyn_tracers%ntrac_met, 1)=
     &    reg%vname(i)
       gfs_dyn_tracers%vname(i+gfs_dyn_tracers%ntrac_met, 2)=
     &    reg%vname(i) // '_q'
       gfs_dyn_tracers%vname(i+gfs_dyn_tracers%ntrac_met, 3)=
     &    reg%vname(i) // '_m'
       gfs_dyn_tracers%vname(i+gfs_dyn_tracers%ntrac_met, 4)=
     &    reg%vname(i) // '_q6'
       gfs_dyn_tracers%vname(i+gfs_dyn_tracers%ntrac_met, 5)=
     &    reg%vname(i) // '_m6'
       gfs_dyn_tracers%cpi(i+gfs_dyn_tracers%ntrac_met) = 0.
       gfs_dyn_tracers%ri(i+gfs_dyn_tracers%ntrac_met) = 0.
      enddo

      endif

! Destroy Chem_Registry
      call Chem_RegistryDestroy ( reg, ierr )

      return

999   print *,'LU_TRC: error in allocate gfs_dyn_tracers :',status,me

      end subroutine dyn_gocart_tracer_config

