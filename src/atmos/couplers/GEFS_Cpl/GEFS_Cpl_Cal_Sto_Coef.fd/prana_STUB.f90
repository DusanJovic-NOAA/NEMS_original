! MODULE PRANA
! Generate a random proper rotation (orthogonal) matrix. 
!
! The procedure is to fill a square matrix of desired order with independent
! Gaussian random numbers, then to apply Gram-Schmidt orthonormalization.
! The "4.... goto 4" loop allows for reiteration in the very rare case that
! the initial random matrix is (to within round-off) essentially singular.
! The final `coin toss' randomization of orientation of each column along
! the line it lies along overcomes what would otherwise be a bias that
! comes from the way the gram-schmidt code works. The even number of reversals
! ensures the preservation of the determinant.
!=============================================================================
module prana
!=============================================================================
use pran; use peuc
implicit none
private
public:: ranrot
interface ranrot; module procedure ranrot_s,ranrot_d; end interface

contains

!=============================================================================
subroutine ranrot_s(a)
!=============================================================================
real(4),dimension(:,:),intent(OUT)    :: a
!-----------------------------------------------------------------------------
integer                            :: n,i,j,nrank,flip
real(4),dimension(size(a,1),size(a,1)):: as
real(4)                               :: det,xran
end subroutine ranrot_s

!=============================================================================
subroutine ranrot_d(a)
!=============================================================================
real(8),dimension(:,:),intent(OUT)    :: a
!-----------------------------------------------------------------------------
integer                               :: n,i,j,nrank,flip,k
real(8),dimension(size(a,1),size(a,1)):: as
real(8)                               :: det,xran
end subroutine ranrot_d

end module prana
