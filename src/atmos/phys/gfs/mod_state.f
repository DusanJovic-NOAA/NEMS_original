      module mod_state
!
!    New module to supply domain information to the GFS output routines
!    called by wrtout.
!
      use machine,    ONLY: kind_io4
      implicit none
!
!
      real(kind=kind_io4), allocatable :: buff_mult_piece(:),
     1                                    buff_mult_pieces(:)
      real(kind=kind_io4),allocatable,target :: buff_mult_piecef(:,:,:),
     1                                    buff_mult_piecesf(:,:,:,:)
!jw
      real(kind=kind_io4), allocatable,target ::
!jw     1                                    buff_mult_piecea(:,:,:),
     1                                    buff_mult_piecea2d(:,:,:),
     1                                    buff_mult_piecea3d(:,:,:)
!jw      REAL(KIND=KIND_io4) ,pointer ::  buff_mult(:,:,:)
!jw
      integer , allocatable :: ivar_global(:),ivar_global_a(:,:)
     &,                        ivarg_global(:),ivarg_global_a(:,:)
      integer , allocatable :: maskss(:,:,:)
!
      integer ngrid ,ngrida,ngridg
!jw
      integer ngrid2d,ngrid3d

      end module mod_state

