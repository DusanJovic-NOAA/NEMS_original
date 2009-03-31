!-----------------------------------------------------------------------
                        module module_exchange
!-----------------------------------------------------------------------
!
!***  MODULE_EXCHANGE contains the halo exchange routines.  There is a
!***  unique routine for every combination of 2-D and 3-D real array
!***  exchanges being done at one time.  Each subroutine name begins with
!***  "exch" which is then followed by a string of integers.  Each "2" 
!***  in the string indicates exchange being done for a 2-D array.  
!***  Similarly each "3" or "4" in the string indicates exchange being done
!***  for a 3-D or 4-D  array.
!***  Currently there are routines for these combinations:
!***
!***  2, 22, 23, 223, 3, 33, 333, 3333, 4
!
!***  A generic interface exists so that all of the routines
!***  may be called with the name "halo_exch".  If new routines
!***  are added because new combinations are needed then also
!***  add the routine's name to the interface block.
!
!
!***  Buffer arrays are used during the exchange process.  Set the size
!***  below in the parameter ibufexch.  If an error occurs where the 
!***  MPI library indicates that the receive buffer is too small then
!***  increase the size of ibufexch.
!
!***  The 4-element IHANDLE array is used for the nonblocking requests
!***  for all the ISENDS/IRECVS and their MPI_WAITS.  Here is the key
!***  to their use:
!***
!***  IRECV/store from north --> IHANDLE(1)
!***  IRECV/store from south --> IHANDLE(2)
!***  ISEND to north --> IHANDLE(3)
!***  ISEND to south --> IHANDLE(4)
!***
!***  IRECV/store from west --> IHANDLE(1)
!***  IRECV/store from east --> IHANDLE(2)
!***  ISEND to east --> IHANDLE(3)
!***  ISEND to west --> IHANDLE(4)
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      interface halo_exch
        module procedure exch2
        module procedure exch22
        module procedure exch23
        module procedure exch223
        module procedure exch3
        module procedure exch33
        module procedure exch333
        module procedure exch3333
        module procedure exch4
      end interface
!
!-----------------------------------------------------------------------
!
      contains
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch2(arr1,ll1,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for a single 2-D array.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 2-D array (=1)
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
!
 arr1                   ! array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!
      end subroutine exch2
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch3(arr1,ll1,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for a single 3-D array.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
!
 arr1                   ! array whose haloes are exchanged
!
!
!-----------------------------------------------------------------------
!
      end subroutine exch3
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch23(arr1,ll1,arr2,ll2,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for 2-D and 3-D arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 2-D array (=1)
,ll2 &                  ! vertical dimension of 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1                   ! 2-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-D array whose haloes are exchanged
!
!
!-----------------------------------------------------------------------
!
      end subroutine exch23
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch22(arr1,ll1,arr2,ll2,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for two 2-D arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 2-D array (=1)
,ll2 &                  ! vertical dimension of 2nd 2-D array (=1)
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1 &                 ! 2-D array whose haloes are exchanged
,arr2                   ! 2-D array whose haloes are exchanged
!
!
!-----------------------------------------------------------------------
!
      end subroutine exch22
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch223(arr1,ll1,arr2,ll2,arr3,ll3,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for two 2-D arrays and one 3-D array.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 2-D array (=1)
,ll2 &                  ! vertical dimension of 2nd 2-D array (=1)
,ll3 &                  ! vertical dimension of 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1 &                 ! 2-D array whose haloes are exchanged
,arr2                   ! 2-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-D array whose haloes are exchanged
!
!
!-----------------------------------------------------------------------
!
      end subroutine exch223
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch22333(arr1,ll1,arr2,ll2,arr3,ll3,arr4,ll4,arr5,ll5 &
                          ,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for two 2-D arrays and three 3-D arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 2-D array (=1)
,ll2 &                  ! vertical dimension of 2nd 2-D array (=1)
,ll3 &                  ! vertical dimension of 1st 3-D array
,ll4 &                  ! vertical dimension of 2nd 3-D array
,ll5 &                  ! vertical dimension of 3rd 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1 &                 ! 2-D array whose haloes are exchanged
,arr2                   ! 2-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll4),intent(inout) :: &
 arr4                   ! 3-D array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll5),intent(inout) :: &
 arr5                   ! 3-D array whose haloes are exchanged
!
!
!-----------------------------------------------------------------------
!
      end subroutine exch22333
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch33(arr1,ll1,arr2,ll2,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes two 3-D real arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 3-D array
,ll2 &                  ! vertical dimension of 2nd 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
 arr1                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-D array whose haloes are exchanged
!
!
!-----------------------------------------------------------------------
!
      end subroutine exch33
!
!--------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch333(arr1,ll1,arr2,ll2,arr3,ll3,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes three 3-D real arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 3-D array
,ll2 &                  ! vertical dimension of 2nd 3-D array
,ll3 &                  ! vertical dimension of 3rd 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
 arr1                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-D array whose haloes are exchanged
!
!
!-----------------------------------------------------------------------
!
      end subroutine exch333
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch3333(arr1,ll1,arr2,ll2,arr3,ll3,arr4,ll4,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes four 3-D real arrays.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 3-D array
,ll2 &                  ! vertical dimension of 2nd 3-D array
,ll3 &                  ! vertical dimension of 3rd 3-D array
,ll4 &                  ! vertical dimension of 4th 3-D array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
 arr1                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-D array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll4),intent(inout) :: &
 arr4                   ! 3-D array whose haloes are exchanged
!
!
!-----------------------------------------------------------------------
!
      end subroutine exch3333
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch4(arr1,ll1,nl1,nstart,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  Exchange haloes for a single 4-D array.
!
!-----------------------------------------------------------------------
!***  Argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of array
,nl1 &                  ! 4th dimension of array
,nstart &               ! index of the 4th dimension to start exchange
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1,nl1),intent(inout) :: &
!
 arr1                   ! array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!
      end subroutine exch4
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------
!
      end module module_exchange
!
!-----------------------------------------------------------------------
