

#define VERIFY_(A)   if(  A/=0) then; if(present(rc)) rc=A; PRINT *, __LINE__; return; endif
#define ASSERT_(A)   if(.not.A) then; if(present(rc)) rc=1; PRINT *, __LINE__; return; endif
#define RETURN_(A)   if(present(rc)) rc=A; return
#define SUCCESS      0
#define DEALOC_(A)   if(associated(A)) deallocate(A)


!BOP

! !MODULE: MAPL_HorzTransMod
!    A Module to do linear transformations on 2-dimensional arrays


! !INTERFACE:

module MAPL_HorzTransformMod

!  $Id: MAPL_HorzTransform.F90,v 1.7.20.1 2009/01/21 22:25:04 ltakacs Exp $

  implicit none
  private

! PUBLIC TYPES:

  public MAPL_HorzTransform

! PUBLIC MEMBER FUNCTIONS:

  public MAPL_HorzTransformCreate
  public MAPL_HorzTransformDestroy
  public MAPL_HorzTransformRun

  public MAPL_DimTopoCyclic
  public MAPL_DimTopoEdge
  public MAPL_DimTopoCenter


! !DESCRIPTION:

! This package performs a serial linear transformation with a
!   previously computed sparse matrix stored in a MAPL_HorzBinTransform
!   object.  Currently, the Transform objects that can be created
!   are limited to a conservative ``binning'',
!   bilinear intepolation, or biquadratic interpolation
!   of uniform-resolution, 2-dimensional
!   (LogicalRectangular) grids that span the same domain. 
!
! Either dimension of the array can be cyclic or bounded. If bounded
!   the grids can have either an edge or a center at the domain boundary,
!   the only difference being that an "edge" boundary does not allow
!   any extrapolation on input values below the first or above the last
!   point. Edge is useful only for meridional interpolation between the
!   two poles.
!   
! To use the pachage, a transform is first created using {\bf MAPL_HorzTransformCreate} and
!   then used repeatedly by passing it to the interpolation routine {\bf MAPL_HorzTransformRun}
!   together with input and output arrays. When a transfor is no longer needed, it can be
!   destroyed with {\bf MAPL_HorzTransformDestroy}. 
!
! In addition to these three methods and the type of the class object {\bf MAPL_HorzTransform},
!   the only puclic quantities are the allowed values of an enumeration that   
!   specifies the 'topology' of a dimension (cyclic, edge, or centered).
!   
! The package is not parallelized, so all the data must be on a single processor. It  
!   does not rely on an other modules or libraries.
!   
!EOP

! Topology enumeration:
!   Edge and Center refer to the location of the 
!   first and last values in a dimension. The distance between points
!   is always the same, but with edge the first and last points represent 
!   only half a grid box, as with the meridional dimension in the grid
!   used by FV.

  integer, parameter :: MAPL_DimTopoCyclic = 0
  integer, parameter :: MAPL_DimTopoEdge   = -1
  integer, parameter :: MAPL_DimTopoCenter = 1

! The code is mostly general for N dimensions, but in fact works only
! for transformations of 2-dimensional arrays.

  integer, parameter :: NUMDIMS = 2

  type Weights
     integer       :: n0,n1
     real, pointer :: f(:)
  end type Weights

  type Mapping
     type(Weights), pointer :: WeightList(:)
  end type Mapping

  type MAPL_HorzTransform
     private
     integer                :: Order
     integer                :: N_in  (NUMDIMS), N_out  (NUMDIMS)
     integer                :: topoIN(NUMDIMS), topoOUT(NUMDIMS)
     real                   :: XmaxIn(NUMDIMS), XmaxOut(NUMDIMS)
     real                   :: XminIn(NUMDIMS), XminOut(NUMDIMS)
     type(Mapping)          :: DimMapping(NUMDIMS)
  end type MAPL_HorzTransform

  interface MAPL_HorzTransformCreate
     module procedure MAPL_HorzTransformCreateBySize
     module procedure MAPL_HorzTransformCreateGEOS
  end interface

  interface MAPL_HorzTransformRun
     module procedure MAPL_HorzTransformRun2
     module procedure MAPL_HorzTransformRun3
  end interface

  real, parameter :: wc = 0.7

contains


  subroutine MAPL_HorzTransformCreateGEOS  (Trans,           &
                                       im_in,  jm_in,   &
                                       im_out, jm_out,  &
                                       XYOFFSET, Order, &
                                                      rc)

    type (MAPL_HorzTransform), intent(OUT) :: Trans
    integer,              intent(IN ) :: im_in, jm_in, im_out, jm_out
    integer, optional,    intent(IN ) :: XYOFFSET
    integer, optional,    intent(IN ) :: Order
    integer, optional,    intent(OUT) :: rc

    integer :: xyoffset_
    integer :: order_

! Specialized interface for use in GEOS5. Input is always on an FV grid;
! if DePe is absent or false, output is also on an FV grid; if DePe is
! present and true, it is on a global, dateline-edge-pole-edge grid. The
! FV grids is dateline-center with .

    order_ = 0
    if(present(Order)) order_ =Order

    Xyoffset_ = 0
    if(present(Xyoffset)) Xyoffset_ = Xyoffset

    select case(Xyoffset_)
    case(0)
       call MAPL_HorzTransformCreateBySize(Trans,                      &
                     (/ im_in              , jm_in            /), &
                     (/ im_out             , jm_out           /), &
                     (/ MAPL_DimTopoCyclic , MAPL_DimTopoEdge /), &
                     (/ MAPL_DimTopoCyclic , MAPL_DimTopoEdge /), &
                     (/ -180.              , -90.             /), &
                     (/ -180.              , -90.             /), &
                     (/  180.-(360./IM_in ),  90.             /), &
                     (/  180.-(360./IM_out),  90.             /), &
                     order_,                                      &
                                                               rc )

    case(1)
       call MAPL_HorzTransformCreateBySize(Trans,                      &
                     (/im_in              , jm_in             /), &
                     (/im_out             , jm_out            /), &
                     (/MAPL_DimTopoCyclic , MAPL_DimTopoEdge  /), &
                     (/MAPL_DimTopoCyclic , MAPL_DimTopoEdge  /), &
                     (/-180.              , -90.              /), &
                     (/-180.+(180./im_out), -90.              /), &
                     (/ 180.-(360./IM_in ),  90.              /), &
                     (/ 180.-(180./IM_out),  90.              /), &
                     order_,                                      &
                                                               rc )
    case(2)
       call MAPL_HorzTransformCreateBySize(Trans,                      &
                     (/im_in              , jm_in             /), &
                     (/im_out             , jm_out            /), &
                     (/MAPL_DimTopoCyclic , MAPL_DimTopoEdge  /), &
                     (/MAPL_DimTopoCyclic , MAPL_DimTopoCenter/), &
                     (/-180.              , -90.              /), &
                     (/-180.              , -90.+(90./jm_out) /), &
                     (/ 180.-(360./IM_in ),  90.              /), &
                     (/ 180.-(360./IM_out),  90.-(90./jm_out) /), &
                     order_,                                      &
                                                               rc )
    case(3)
       call MAPL_HorzTransformCreateBySize(Trans,                      &
                     (/im_in              , jm_in             /), &
                     (/im_out             , jm_out            /), &
                     (/MAPL_DimTopoCyclic , MAPL_DimTopoEdge  /), &
                     (/MAPL_DimTopoCyclic , MAPL_DimTopoCenter/), &
                     (/-180.              , -90.              /), &
                     (/-180.+(180./im_out), -90.+(90./jm_out) /), &
                     (/ 180.-(360./IM_in ),  90.              /), &
                     (/ 180.-(180./IM_out),  90.-(90./jm_out) /), &
                     order_,                                      &
                                                               rc )
    case default
       ASSERT_(.false.)
    end select

  end subroutine MAPL_HorzTransformCreateGEOS



  subroutine MAPL_HorzTransformCreateBySize(Trans,           &
                                       N_in ,  N_out,   &
                                       topoIN, topoOUT, &
                                       FirstIN,FirstOUT,&
                                       LastIN, LastOUT, &
                                       order,           &
                                                    rc)
    type (MAPL_HorzTransform), intent(OUT) :: Trans
    integer,              intent(IN ) :: N_in    (NUMDIMS)
    integer,              intent(IN ) :: N_out   (NUMDIMS)
    integer,              intent(IN ) :: topoIN  (NUMDIMS) 
    integer,              intent(IN ) :: topoOUT (NUMDIMS) 
    real,                 intent(IN ) :: FirstIN (NUMDIMS)
    real,                 intent(IN ) :: LastIN  (NUMDIMS)
    real,                 intent(IN ) :: FirstOUT(NUMDIMS)
    real,                 intent(IN ) :: LastOUT (NUMDIMS)
    integer,              intent(IN ) :: order
    integer, optional,    intent(OUT) :: rc

    integer :: status
    integer :: i, j
    real    :: rngIn, rngOut

    real, allocatable :: xin(:), xout(:)

    if(present(rc)) rc = 0

    Trans%N_in    = N_in
    Trans%N_out   = N_out
    Trans%topoIN  = topoIN
    Trans%topoOUT = topoOUT
    Trans%XminIn  = FirstIN
    Trans%XmaxIn  = LastIN  
    Trans%XminOut = FirstOUT
    Trans%XmaxOut = LastOUT
    Trans%Order   = order

    ASSERT_(all(Trans%N_in >1))
    ASSERT_(all(Trans%N_out>1))

    do i=1,NUMDIMS

       if(topoIN(i)==MAPL_DimTopoEdge) then
          ASSERT_(Trans%XminIn(i)<=Trans%XminOut(i))
          ASSERT_(Trans%XmaxIn(i)>=Trans%XmaxOut(i))
       end if

       if(topoIN(i)==MAPL_DimTopoCyclic) then
          ASSERT_(topoOUT(i)==MAPL_DimTopoCyclic)
       end if

       if(topoOUT(i)==MAPL_DimTopoCyclic) then
          ASSERT_(topoIN(i)==MAPL_DimTopoCyclic)
       end if

       allocate(Trans%DimMapping(i)%WeightList(n_out(i)), stat=STATUS)
       VERIFY_(STATUS)

       if(topoIN(i)==MAPL_DimTopoCyclic) then
          allocate(xin (2*Trans%N_in (i)+1),stat=STATUS)
          VERIFY_(STATUS)
          allocate(xout(  Trans%N_out(i)+1),stat=STATUS)
          VERIFY_(STATUS)
       elseif(order==0) then
          allocate(xin (  Trans%N_in (i)+1),stat=STATUS)
          VERIFY_(STATUS)
          allocate(xout(  Trans%N_out(i)+1),stat=STATUS)
          VERIFY_(STATUS)
       else
          allocate(xin (  Trans%N_in (i)  ),stat=STATUS)
          VERIFY_(STATUS)
          allocate(xout(  Trans%N_out(i)  ),stat=STATUS)
          VERIFY_(STATUS)
       endif

       call GetX(xin , Trans%N_in(i) , Trans%XminIn (i), Trans%XmaxIn (i),  &
            Trans%TopoIn (i), order, rc=status)
       VERIFY_(STATUS)

       call GetX(xout, Trans%N_out(i), Trans%XminOut(i), Trans%XmaxOut(i),  &
            Trans%TopoOut(i), order, rc=status)
       VERIFY_(STATUS)

       if(topoIN(i)==MAPL_DimTopoCyclic) then
          rngIn  = ((Trans%XmaxIn (i)-Trans%XminIn (i))*Trans%N_in (i))/(Trans%N_in (i)-1)
          rngOut = ((Trans%XmaxOut(i)-Trans%XminOut(i))*Trans%N_out(i))/(Trans%N_out(i)-1)

          ASSERT_(abs( (rngIn-rngOut)/rngIn ) < 1.e-5)

          if(xout(1)<xin(1)) then
             xout  = xout + int((xin(1)-xout(1))/rngIn+1)*rngIn
          else
             xout  = xout + int((xin(1)-xout(1))/rngIn)*rngIn
          end if
       end if

       ASSERT_(xin(size(xin)) >= xout(size(xout)))
       ASSERT_(xin(        1) <= xout(         1))

       select case (order)
       case(0)
          call ComputeDimBinWeights(Trans%DimMapping(i)%WeightList,Xin,Xout, &
           HasPoles=(topoIN(i)==MAPL_DimTopoEdge.and.topoOUT(i)==MAPL_DimTopoEdge),&
           rc=status)
          VERIFY_(STATUS)
       case(1)
          call ComputeDimLinWeights(Trans%DimMapping(i)%WeightList,Xin,Xout,rc=status)
          VERIFY_(STATUS)
       case default
          ASSERT_(.false.)
       end select

       ASSERT_(all(Trans%DimMapping(i)%WeightList(:)%n0<=2*Trans%N_in(i)))

       deallocate(Xin )
       deallocate(Xout)

    end do

    RETURN_(SUCCESS)

  end subroutine MAPL_HorzTransformCreateBySize

  subroutine GetX(X,N,xmin,xmax,topo,order,rc)

    real,                 intent(OUT) :: x(:)
    integer,              intent(IN ) :: N
    integer,              intent(IN ) :: topo
    integer,              intent(IN ) :: order
    real,                 intent(IN ) :: xmin
    real,                 intent(IN ) :: xmax
    integer, optional,    intent(OUT) :: rc

    integer :: status, j, jm
    real    :: dx

    jm = size(X)
    dx = (Xmax-xmin) / (N-1)

    if(order==1) then
       x(1)       = xmin
    else
       x(1)       = xmin-0.5*dx
    end if

    do j=2,JM
       x(j) = x(1) + (j - 1)*dx
    end do

    if(topo==MAPL_DimTopoEdge  ) then
       x( 1) = xmin
       x(JM) = xmax
    end if

    RETURN_(SUCCESS)
  end subroutine GetX


  subroutine ComputeDimLinWeights(Weight,Xin,Xout,rc)
    
    type(Weights),     intent(INOUT) :: Weight(:)
    real,              intent(IN   ) :: Xin(:), Xout(:)
    integer, optional, intent(OUT  ) :: rc


    ! Compute weights for binned interpolation along a dimension.
    ! Xout are the N_in + 1 input bin edges.
    ! Xin  are the N_out + 1 output bin edges
    ! Weigths are the mapping


    integer :: j_out, j0, j1
    integer :: N_in
    integer :: status
    real, pointer :: b(:)

    N_in  = size(Xin )

    do j_out=1,size(Weight)
       j0 = 1
       do           
          if(Xout(j_out  )<=Xin(j0+1)) exit
          j0=j0+1
          ASSERT_(j0 < N_in )
       end do
       j1 = j0 + 1

       allocate(b(j0:j1), stat=STATUS)
       VERIFY_(STATUS)

       b(j0  ) = (Xin(j1)-Xout(j_out))/(Xin(j1)-Xin(j0))
       b(j0+1) = 1.0 - b(j0)

       Weight(j_out)%n0 =  j0
       Weight(j_out)%n1 =  j1
       Weight(j_out)%f  => b

    end do

  end subroutine ComputeDimLinWeights

  subroutine ComputeDimBinWeights(Weight,Xin,Xout,HasPoles,rc)
    
    type(Weights),     intent(INOUT) :: Weight(:)
    real,              intent(IN   ) :: Xin(:), Xout(:)
    logical,           intent(IN   ) :: HasPoles
    integer, optional, intent(OUT  ) :: rc


    ! Compute weights for binned interpolation along a dimension.
    ! Xout are the N_in + 1 input bin edges.
    ! Xin  are the N_out + 1 output bin edges
    ! Weigths are the mapping


    integer :: j_out, j0, j1, j
    integer :: N_in, N_out
    integer :: status
    real    :: dx, ff
    real, pointer :: b(:)

    N_in  = size(Xin )-1
    N_out = size(Weight)

    do j_out=1,N_out
       j0 = 1
       do           
          if(Xout(j_out  )>=Xin(j0) .and. Xout(j_out  )<=Xin(j0+1)) exit
          j0=j0+1
          ASSERT_(j0 <= N_in )
       end do

       j1 = j0
       do
          if(Xout(j_out+1)>=Xin(j1) .and. Xout(j_out+1)<=Xin(j1+1)) exit
          j1=j1+1
          ASSERT_(j1 <= N_in )
       end do

       allocate(b(j0:j1), stat=STATUS)
       VERIFY_(STATUS)

       if(j0==j1) then
          b(j0) = 1.
       else
          dx    = Xin(j0+1)-Xout(j_out)
          ff    = dx
          b(j0) = dx
          do j=j0+1,j1-1
             dx   = Xin(j+1) - Xin(j)
             ff   = ff + dx
             b(j) = dx
          end do
          dx    = Xout(j_out+1)-Xin(j1)
          ff    = ff + dx
          b(j1) = dx
          b     = b/ff
       end if

       Weight(j_out)%n0 =  j0
       Weight(j_out)%n1 =  j1
       Weight(j_out)%f  => b

    end do

    if(HasPoles) then
       Weight(    1)%n0 =  1
       Weight(    1)%n1 =  1
       Weight(    1)%f  =  1.

       Weight(N_out)%n0 =  N_in
       Weight(N_out)%n1 =  N_in
       Weight(N_out)%f  =  1.
    endif

  end subroutine ComputeDimBinWeights

  subroutine DestroyMapping(MAP, rc)
    type (Mapping),    intent(INOUT) :: Map
    integer, optional, intent(  OUT) :: rc

    integer :: i
    
    do i=1,size(MAP%WeightList)
       DEALOC_(MAP%WeightList(i)%f)
    end do

    DEALOC_(Map%Weightlist)

    RETURN_(SUCCESS)
  end subroutine DestroyMapping

  subroutine MAPL_HorzTransformDestroy(Trans, rc)
    type (MAPL_HorzTransform), intent(INOUT) :: Trans
    integer, optional,    intent(  OUT) :: rc

    integer :: status
    integer :: i, j

    do I=1,NUMDIMS
       call DestroyMapping(Trans%DimMapping(i), RC=STATUS)
       VERIFY_(STATUS)
    enddo

    RETURN_(SUCCESS)
  end subroutine MAPL_HorzTransformDestroy

  subroutine MAPL_HorzTransformRun2(Trans, qin, qout, undef, gridtype, rc)

!    Trans ....  Precomputed transform
!      qin ....  Input Variable
!      qout....  Output Variable
!    undef ....  UNDEF Value

    type (MAPL_HorzTransform), intent(IN   ) :: Trans
    real,                    intent(IN   ) :: qin (:,:)
    real,                    intent(OUT  ) :: qout(:,:)
    real,    optional,       intent(IN   ) :: undef
    character(len=*), optional, intent(IN) :: gridtype
    integer, optional,       intent(  OUT) :: rc

    real            :: undef_
    real            :: q, w, f
    integer         :: STATUS
    integer         :: i,j
    integer         :: i0,i1,j0,j1
    integer         :: ii,jj,jx,ix
    real, pointer   :: fx(:), fy(:)
    
    logical :: doCube2Latlon
    integer :: npx, npy
    integer :: nlon, nlat

    if(present(rc)) rc = 0
    
    if(present(undef)) then
       undef_ = undef
    else
       undef_ = huge(undef_)
    end if

    ASSERT_(all(shape(qin )==Trans%N_in ))
    ASSERT_(all(shape(qout)==Trans%N_out))

    doCube2Latlon=.false.
    if (present(gridtype)) then
       if (gridtype=='Cubed-Sphere') doCube2Latlon=.true.
    endif

    if (doCube2Latlon) then 

       npx  = Trans%N_in(1)
       npy  = Trans%N_in(2)
       nlon = Trans%N_out(1)
       nlat = Trans%N_out(2)
       call cube2latlon(npx, npy, nlon, nlat, qin, qout)

    else

    do j=1,Trans%N_out(2)
       j0 = Trans%DimMapping(2)%WeightList(j)%n0
       j1 = Trans%DimMapping(2)%WeightList(j)%n1
       fy =>Trans%DimMapping(2)%WeightList(j)%f

       do i=1,Trans%N_out(1)
          i0 = Trans%DimMapping(1)%WeightList(i)%n0
          i1 = Trans%DimMapping(1)%WeightList(i)%n1
          fx =>Trans%DimMapping(1)%WeightList(i)%f

          q = 0.0
          w = 0.0

          do jj=j0,j1
             if(jj>Trans%N_in(2)) then
                jx = jj - Trans%N_in(2)
             else
                jx = jj
             end if

             do ii=i0,i1
                if(ii>Trans%N_in(1)) then
                   ix = ii - Trans%N_in(1)
                else
                   ix = ii
                end if

                if(qin(ix,jx) /= undef_) then
                   f = fx(ii)*fy(jj)
                   q = q + f*qin(ix,jx)
                   w = w + f           
                end if
             end do
          end do

          if ( w >= wc ) then
             qout(i,j) = q / w
          else
             qout(i,j) = undef_
          end if

       end do
    end do

    endif

    RETURN_(SUCCESS)
  end subroutine MAPL_HorzTransformRun2

  subroutine MAPL_HorzTransformRun3(Trans, qin, qout, undef, gridtype, rc)

!    Trans ....  Precomputed transform
!      qin ....  Input Variable
!      qout....  Output Variable
!    undef ....  UNDEF Value

    type (MAPL_HorzTransform), intent(IN   ) :: Trans
    real,                    intent(IN   ) :: qin (:,:,:)
    real,                    intent(OUT  ) :: qout(:,:,:)
    real,    optional,       intent(IN   ) :: undef
    character(len=*), optional, intent(IN) :: gridtype
    integer, optional,       intent(  OUT) :: rc

    real          :: undef_
    real          :: q(size(qin,3))
    real          :: w(size(qin,3))
    real          :: f
    integer       :: STATUS
    integer       :: i,j,k
    integer       :: i0,i1,j0,j1
    integer       :: ii,jj,jx,ix
    integer       :: sh(3)

    logical :: doCube2Latlon
    integer :: npx, npy, npz
    integer :: nlon, nlat

    real, pointer :: fx(:), fy(:)

    if(present(rc)) rc = 0
    
    if(present(undef)) then
       undef_ = undef
    else
       undef_ = huge(undef_)
    end if

    sh = shape(qin )
    ASSERT_(all(sh(1:2)==Trans%N_in ))
    sh = shape(qout)
    ASSERT_(all(sh(1:2)==Trans%N_out))

    ASSERT_(size(qin,3)==size(qout,3))

    doCube2Latlon=.false.
    if (present(gridtype)) then
       if (gridtype=='Cubed-Sphere') doCube2Latlon=.true.
    endif

    if (doCube2Latlon) then

       npx  = Trans%N_in(1)
       npy  = Trans%N_in(2)
       npz  = Trans%N_in(3)
       nlon = Trans%N_out(1)
       nlat = Trans%N_out(2)
       do k=1,npz
          call cube2latlon(npx, npy, nlon, nlat, qin(:,:,k), qout(:,:,k))
       enddo

    else

    do j=1,Trans%N_out(2)
       j0 = Trans%DimMapping(2)%WeightList(j)%n0
       j1 = Trans%DimMapping(2)%WeightList(j)%n1
       fy =>Trans%DimMapping(2)%WeightList(j)%f

       do i=1,Trans%N_out(1)
          i0 = Trans%DimMapping(1)%WeightList(i)%n0
          i1 = Trans%DimMapping(1)%WeightList(i)%n1
          fx =>Trans%DimMapping(1)%WeightList(i)%f

          q = 0.0
          w = 0.0

          do jj=j0,j1
             if(jj>Trans%N_in(2)) then
                jx = jj - Trans%N_in(2)
             else
                jx = jj
             end if

             do ii=i0,i1
                if(ii>Trans%N_in(1)) then
                   ix = ii - Trans%N_in(1)
                else
                   ix = ii
                end if

                f = fx(ii)*fy(jj)

                where(qin(ix,jx,:) /= undef_)
                   q = q + f*qin(ix,jx,:)
                   w = w + f           
                end where
             end do
          end do

          where( w >= wc )
             qout(i,j,:) = q / w
          elsewhere
             qout(i,j,:) = undef_
          end where

       end do
    end do

    endif

    RETURN_(SUCCESS)
  end subroutine MAPL_HorzTransformRun3

end module MAPL_HorzTransformMod
