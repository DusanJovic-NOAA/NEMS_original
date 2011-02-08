!
!! ! Module: gfs_phy_tracer_config
!
! ! Description: gfs physics tracer configuration module
!
! ! Revision history:
!   Oct 16 2009   Sarah Lu, adopted from dyn fc
!   Nov 21 2009   Sarah Lu, chem tracer specified from ChemRegistry
!   Dec 10 2009   Sarah Lu, add doing_GOCART
!   Jan 12 2010   Sarah Lu, add trcindx
!   Feb 08 2009   Sarah Lu, ri/cpi added to gfs_phy_tracer_type
!   Aug 17 2010   Sarah Lu, remove debug print
!   Oct 16 2010   Sarah Lu, add fscav
! -------------------------------------------------------------------------
!
      module gfs_phy_tracer_config
      use machine , only : kind_phys
      use tracer_const, only : cpi,ri

      use Chem_RegistryMod
      implicit none
      SAVE
!
! tracer specification
!
      type    gfs_phy_tracer_type
        character*20        , pointer      :: vname(:)    ! variable name
        real(kind=kind_phys), pointer      :: ri(:)
        real(kind=kind_phys), pointer      :: cpi(:)
        real(kind=kind_phys), pointer      :: fscav(:)
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
      logical                              :: glbsum  = .true.
!

! --- public interface
      public     tracer_config_init, trcindx

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

!!    if ( me == 0) call Chem_RegistryPrint (reg)

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

! Set up tracer name and allocate ri, cpi, fscav
      if ( gfs_phy_tracer%ntrac > 0 ) then      
       allocate(gfs_phy_tracer%vname(ntrac), stat=status)
           if( status .ne. 0 ) go to 999         
       allocate(gfs_phy_tracer%ri(0:ntrac), stat=status)
           if( status .ne. 0 ) go to 999         
       allocate(gfs_phy_tracer%cpi(0:ntrac), stat=status)
           if( status .ne. 0 ) go to 999         
       allocate(gfs_phy_tracer%fscav(ntrac), stat=status)
           if( status .ne. 0 ) go to 999         

!--- fill in met tracers
      gfs_phy_tracer%vname(1) = 'spfh'   
      gfs_phy_tracer%vname(ntoz) = 'o3mr'  
      gfs_phy_tracer%vname(ntcw) = 'clwmr' 
!--- fill in chem tracers
      do i = 1,gfs_phy_tracer%ntrac_chem
       gfs_phy_tracer%vname(i+gfs_phy_tracer%ntrac_met)=reg%vname(i)
      enddo
!--- fill in default values for fscav
      gfs_phy_tracer%fscav(:) = 0.

      endif

! Destroy Chem_Registry
      call Chem_RegistryDestroy ( reg, ierr )

      return

999   print *,'LU_TRC: error in allocate gfs_phy_tracer :',status,me

      end subroutine tracer_config_init

! -------------------------------------------------------------------
! -------------------------------------------------------------------
      function trcindx( specname, tracer )
      implicit none

      character*(*), intent(in)  ::  specname
      type (gfs_phy_tracer_type), intent(in)    ::  tracer

      character*10  ::  name1, name2
      integer       ::  i, trcindx

! -- set default value
      trcindx = -999

! -- convert specname to upper case
      call fixchar(specname, name1, 1)
      do i = 1, tracer%ntrac
        call fixchar(tracer%vname(i), name2, 1)
        if( name1 == name2 ) then
          trcindx = i 
          exit
        endif
      enddo

      return
      end function trcindx

! -------------------------------------------------------------------
      subroutine fixchar(name_in, name_out, option)
      implicit none

      character*(*), intent(in)   ::  name_in
      character*(*), intent(out)  ::  name_out
      integer, intent(in)         ::  option

      character*10                :: temp
      integer                     :: i, ic

      name_out= '          '
      temp = trim(adjustl(name_in))
      do i = 1, len_trim(temp)
        ic = IACHAR(temp(i:i))
        if(option == 1 ) then             !<--- convert to upper case
          if(ic .ge. 97 .and. ic .le. 122) then
            name_out(i:i) = CHAR( IC-32 )
          else
            name_out(i:i) = temp(i:i)
          endif
        endif
        if(option == 2 ) then             !<--- convert to lower case
          if(ic .ge. 65 .and. ic .le. 90) then
            name_out(i:i) = CHAR( IC+32 )
          else
            name_out(i:i) = temp(i:i)
          endif
        endif

      enddo
      name_out=trim(name_out)
      return

      end subroutine fixchar

! ========================================================================= 

      end module gfs_phy_tracer_config
