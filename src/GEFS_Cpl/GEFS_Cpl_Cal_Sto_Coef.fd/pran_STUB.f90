module pran
implicit none
private:: uniform_d,expone_d,gauss_d,gammar,gammar_d,gammad_d,           &
     chisqd,chisqd_d,chisqr_d,gammln_d,gammq_d,gcf_d,gser_d,fgam,fgam_d, &
     randomrot_d,setrot,setrot_d,outerprod3,outerprod3_d

interface chisqq
	module procedure chisqq,chisqq_d, chisqi, chisqi_d
end interface
interface promptseed;module procedure promptseed;        end interface
interface plant;   module procedure plant, plant1;       end interface
interface skewkurt;module procedure skewkurt, skewkurt_d;end interface
interface gammq;   module procedure gammq,  gammq_d;     end interface
interface gammln;  module procedure gammln, gammln_d;    end interface
interface uniform; module procedure uniform,uniform_d;   end interface
interface expone;  module procedure expone, expone_d;    end interface
interface gauss;   module procedure gauss,  gauss_d;     end interface
interface gcf;     module procedure gcf,    gcf_d;       end interface
interface gser;    module procedure gser,   gser_d;      end interface
interface fgam;    module procedure fgam,   fgam_d;      end interface
interface gammad
	module procedure gammar, gammar_d, gammad, gammad_d
end interface
interface chisqr  
   module procedure chisqr, chisqr_d, chisqd, chisqd_d
end interface
interface ranperm; module procedure ranperm;             end interface
interface randomrot; module procedure randomrot,randomrot_d; end interface
contains
subroutine skewkurt(u0,u1,u2,u3,u4,bias,sdev,skew,kurt)
real(4),intent(IN ):: u0,u1,u2,u3,u4
real(4),intent(OUT):: bias,sdev,skew,kurt
real(4)            :: t1,t2,t3,t4,s2
end subroutine skewkurt

subroutine skewkurt_d(u0,u1,u2,u3,u4,bias,sdev,skew,kurt)
real(8),intent(IN ):: u0,u1,u2,u3,u4
real(8),intent(OUT):: bias,sdev,skew,kurt
real(8)            :: t1,t2,t3,t4,s2
end subroutine skewkurt_d

subroutine gammln(g,xx)
real(4),intent(IN) :: xx
real(4),intent(OUT):: g
REAL(8)         :: cof(6),stp,half,one,fpf,x,tmp,ser
integer         :: j
DATA cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0, &
     -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
DATA half,one,fpf/0.5d0,1.0d0,5.5d0/
END subroutine gammln

subroutine gammln_d(g,xx)
real(8),intent(IN ) :: xx
real(8),intent(OUT) :: g
REAL(8)             :: cof(6),stp,half,one,fpf,x,tmp,ser
integer             :: j
END subroutine gammln_d

subroutine gammq(g,a,x)
REAL(4),INTENT(IN) :: a,x
real(4),intent(OUT):: g
REAL(4)            :: gln,gammcf,gamser
END subroutine gammq

subroutine gammq_d(g,a,x)
REAL(8),INTENT(IN ):: a,x
real(8),intent(OUT):: g
REAL(8)            :: gln,gammcf,gamser
END subroutine gammq_d

subroutine promptseed
integer::seed
end subroutine promptseed

subroutine plant
integer,parameter :: nn=2
integer,dimension(nn):: jseed
integer:: iran,jran
end subroutine plant

subroutine plant1(iran)
integer,intent(IN):: iran
integer,parameter :: nn=2
integer,dimension(nn):: jseed
integer:: jran
end subroutine plant1

subroutine uniform(x)
real(4),intent(OUT):: x
real(4)            :: x1,x2
end subroutine uniform

subroutine uniform_d(x)
real(8),intent(OUT):: x
real(8)            :: x1,x2
end subroutine uniform_d

subroutine expone(x)
real(4),intent(OUT):: x
call uniform(x)
x=-log(x)
end subroutine expone

subroutine expone_d(x)
real(8),intent(OUT):: x
end subroutine expone_d

subroutine gauss(x)
real(4),intent(OUT):: x
integer:: iset
real(4):: gset,v1,v2,r,f
end subroutine gauss

subroutine gauss_d(x)
real(8),intent(OUT):: x
integer            :: iset
real(8)            :: gset,v1,v2,r,f
end subroutine gauss_d

subroutine gammad(ia,x)
integer,intent(IN ):: ia
real(4),   intent(OUT):: x
integer:: j
real(4):: am,e,s,v1,v2,y
end subroutine gammad

subroutine gammad_d(ia,x)
integer,intent(IN ):: ia
real(8),intent(OUT):: x
integer            :: j
real(8)            :: am,e,s,v1,v2,y
end subroutine gammad_d

subroutine chisqd(ia,x)
integer,intent(IN ):: ia
real(4),   intent(OUT):: x
integer            :: ih
real(4)               :: y
end subroutine chisqd

subroutine chisqd_d(ia,x)
integer,intent(IN ):: ia
real(8),intent(OUT):: x
integer            :: ih
real(8)            :: y
end subroutine chisqd_d

subroutine chisqr(a,x)
real(4),intent(IN ):: a
real(4),intent(OUT):: x
real(4)            :: ah
end subroutine chisqr

subroutine chisqr_d(a,x)
real(8),intent(IN ):: a
real(8),intent(OUT):: x
real(8)            :: ah
end subroutine chisqr_d

SUBROUTINE gcf(gammcf,a,x,gln)
INTEGER, PARAMETER :: itmax=100
REAL(4),    PARAMETER :: eps=3.e-7
REAL(4),   INTENT(out):: gammcf,gln
REAL(4),   INTENT(in ):: a,x
INTEGER            :: n
REAL(4)               :: gold,a0,a1,b0,b1,fac,an,ana,anf,g 
END SUBROUTINE gcf

SUBROUTINE gcf_d(gammcf,a,x,gln)
INTEGER, PARAMETER :: itmax=100
REAL(8), PARAMETER :: eps=3.d-14
REAL(8),INTENT(out):: gammcf,gln
REAL(8),INTENT(in ):: a,x
INTEGER            :: n
REAL(8)            :: gold,a0,a1,b0,b1,fac,an,ana,anf,g 
END SUBROUTINE gcf_d

SUBROUTINE gser(gamser,a,x,gln)
INTEGER,PARAMETER  :: itmax=100
REAL(4),   PARAMETER  :: eps=3.e-7
REAL(4),   INTENT(IN ):: a,x
REAL(4),   INTENT(OUT):: gamser,gln
INTEGER            :: n
REAL(4)               :: ap,sum,del 
END SUBROUTINE gser

SUBROUTINE gser_d(gamser,a,x,gln)
INTEGER,PARAMETER  :: itmax=100
REAL(8),PARAMETER  :: eps=3.d-14
REAL(8),INTENT(IN ):: a,x
REAL(8),INTENT(OUT):: gamser,gln
INTEGER            :: n
REAL(8)            :: ap,sum,del 
END SUBROUTINE gser_d

subroutine fgam(a,x)
! Generate gamma deviate with fractional degrees of freedom
real(4),intent(IN) :: a
real(4),intent(OUT):: x
real(4)            :: z,one,c,d,am,ai
end subroutine fgam

subroutine fgam_d(a,x)
! Generate gamma deviate with fractional degrees of freedom
real(8),intent(IN) :: a
real(8),intent(OUT):: x
real(8)            :: z,one,c,d,am,ai
end subroutine fgam_d

subroutine gammar(a,x)
real(4),intent(IN ):: a
real(4),intent(OUT):: x
integer         :: ia
real(4)            :: r,z
end subroutine gammar

subroutine gammar_d(a,x)
real(8),intent(IN ):: a
real(8),intent(OUT):: x
integer            :: ia
real(8)            :: r,z
end subroutine gammar_d

subroutine chisqq(g,a,x)
real(4),intent(IN ):: a,x
real(4),intent(OUT):: g
real(4):: b,y

end subroutine chisqq

subroutine chisqq_d(g,a,x)
real(8),intent(IN ):: a,x
real(8),intent(OUT):: g
real(8):: b,y
end subroutine chisqq_d

subroutine chisqi(g,ia,x)
integer,intent(IN):: ia
real(4),intent(IN ):: x
real(4),intent(OUT):: g
real(4):: b,y
end subroutine chisqi

subroutine chisqi_d(g,ia,x)
integer,intent(IN ):: ia
real(8),intent(IN ):: x
real(8),intent(OUT):: g
real(8):: b,y
end subroutine chisqi_d

!=============================================================================
subroutine ranperm(n,p)
!=============================================================================
! Deliver a random permutation of {1:n} in array p
!=============================================================================
implicit none
integer,             intent(IN ):: n
integer,dimension(n),intent(OUT):: p
integer,dimension(n)            :: pool
real(8)                         :: xran
integer                         :: i,j,jc
end subroutine ranperm

!==============================================================================
subroutine randomrot(rot)
!==============================================================================
! Create a random 3*3 rotation matrix
real(4),dimension(3,3),intent(OUT):: rot
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
real(4),dimension(4):: v4
real(4)             :: s,alpha,beta,gamma
integer          :: i
end subroutine randomrot

!==============================================================================
subroutine randomrot_d(rot)
!==============================================================================
! Create a random 3*3 rotation matrix
real(8),dimension(3,3),intent(OUT):: rot
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
real(8),dimension(4):: v4
real(8)             :: s,alpha,beta,gamma
integer             :: i
end subroutine randomrot_d

!==============================================================================
subroutine setrot(alpha,beta,gamma,rot)
!==============================================================================
! Construct a 3*3 orthogonal matrix that represents a rotation by an angle
! gamma counterclockwise about an axis given by the unit vector, "axis3".
! axis3 points in the direction with colatitude alpha and longitude beta.
!==============================================================================
implicit none
real(4),               intent(IN ):: alpha,beta,gamma
real(4),dimension(3,3),intent(OUT):: rot
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
real(4),dimension(3,3)            :: t,q,r,u
real(4),dimension(3)              :: axis1,axis2,axis3,v3
real(4)                           :: ca,sa,cb,sb,cg,sg, &
                                  cca,csa,ssa,ccb,csb,ssb,pi2
end subroutine setrot

!==============================================================================
subroutine setrot_d(alpha,beta,gamma,rot)
!==============================================================================
! Construct a 3*3 orthogonal matrix that represents a rotation by an angle
! gamma counterclockwise about an axis given by the unit vector, "axis3".
! axis3 points in the direction with colatitude alpha and longitude beta.
!==============================================================================
implicit none
real(8),               intent(IN ):: alpha,beta,gamma
real(8),dimension(3,3),intent(OUT):: rot
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
real(8),dimension(3,3)            :: t,q,r,u
real(8),dimension(3)              :: axis1,axis2,axis3,v3
real(8)                           :: ca,sa,cb,sb,cg,sg, &
                                     cca,csa,ssa,ccb,csb,ssb,pi2
end subroutine setrot_d

!=============================================================================
subroutine outerprod3(va,vb,vv)
!=============================================================================
implicit none
real(4),dimension(3),  intent(IN ):: va,vb
real(4),dimension(3,3),intent(OUT):: vv
integer                        :: i
end subroutine outerprod3

!=============================================================================
subroutine outerprod3_d(va,vb,vv)
!=============================================================================
implicit none
real(8),dimension(3),  intent(IN ):: va,vb
real(8),dimension(3,3),intent(OUT):: vv
integer                           :: i
end subroutine outerprod3_d

end module pran
