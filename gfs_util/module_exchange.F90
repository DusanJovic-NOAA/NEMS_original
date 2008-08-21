!-----------------------------------------------------------------------
                        module module_exchange_gfs
!-----------------------------------------------------------------------
!
!***  module_exchange contains the halo exchange routines.  there is a
!***  unique routine for every combination of 2-d and 3-d real array
!***  exchanges being done at one time.  each subroutine name begins with
!***  "exch" which is then followed by a string of integers.  each "2" 
!***  in the string indicates exchange being done for a 2-d array.  
!***  similarly each "3" in the string indicates exchange being done for
!***  a 3-d array.  currently there are routines for these combinations:
!***
!***  2, 23, 223, 3, 33, 333, 3333
!
!***  a generic interface exists so that all of the routines
!***  may be called with the name "halo_exch".  if new routines
!***  are added because new combinations are needed then also
!***  add the routine's name to the interface block.
!
!
!***  buffer arrays are used during the exchange process.  set the size
!***  below in the parameter ibufexch.  if an error occurs where the 
!***  mpi library indicates that the receive buffer is too small then
!***  increase the size of ibufexch.
!
!***  the 4-element ihandle array is used for the nonblocking requests
!***  for all the isends/irecvs and their mpi_waits.  here is the key
!***  to their use:
!***
!***  irecv/store from north --> ihandle(1)
!***  irecv/store from south --> ihandle(2)
!***  isend to north --> ihandle(3)
!***  isend to south --> ihandle(4)
!***
!***  irecv/store from west --> ihandle(1)
!***  irecv/store from east --> ihandle(2)
!***  isend to east --> ihandle(3)
!***  isend to west --> ihandle(4)
!
!-----------------------------------------------------------------------
!
use module_include
use module_dm_parallel_gfs,only : its,ite,jts,jte &
                             ,ims,ime,jms,jme &
                             ,ids,ide,jds,jde &
                             ,mype_share,my_neb,mpi_comm_comp
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer(kind=kint),parameter :: ibufexch=2500000
      integer(kind=kint),private :: mype
!
      real(kind=kfpt),dimension(ibufexch) :: buf0,buf1,buf2,buf3
!
!-----------------------------------------------------------------------
!
      interface halo_exch
        module procedure exch2
        module procedure exch23
        module procedure exch223
        module procedure exch3
        module procedure exch33
        module procedure exch333
        module procedure exch3333
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
!***  exchange haloes for a single 2-d array.
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 2-d array (=1)
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
!
 arr1                   ! array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!
!***  the array of neighbors called my_neb is filled in subroutine
!***  decomp in module_dm_parallel.  recall that the directional
!***  designations for my_neb are:
!
!***      north: 1
!***       east: 2
!***      south: 3
!***       west: 4
!***  northeast: 5
!***  southeast: 6
!***  southwest: 7
!***  northwest: 8
!
!***  if my_neb(n) holds the task id of each neighbor.  if there is
!***  no neighbor due to the presence of a global boundary then the
!***  value of my_neb(n) in that direction is -1.
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  north/south
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!
      if(my_neb(1)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
        call mpi_isend(buf2,ic,mpi_real,my_neb(1),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
        call mpi_isend(buf3,ic,mpi_real,my_neb(3),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  east/west
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
        call mpi_isend(buf2,ic,mpi_real,my_neb(2),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
       call mpi_isend(buf3,ic,mpi_real,my_neb(4),mype &
                     ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
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
!***  exchange haloes for a single 3-d array.
!
!-----------------------------------------------------------------------
!***  argument variables
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
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  north/south
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf2,ic,mpi_real,my_neb(1),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf3,ic,mpi_real,my_neb(3),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  east/west
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf2,ic,mpi_real,my_neb(2),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
       call mpi_isend(buf3,ic,mpi_real,my_neb(4),mype &
                     ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
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
!***  exchange haloes for 2-d and 3-d arrays.
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 2-d array (=1)
,ll2 &                  ! vertical dimension of 3-d array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1                   ! 2-d array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-d array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  north/south
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf2,ic,mpi_real,my_neb(1),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf3,ic,mpi_real,my_neb(3),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  east/west
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf2,ic,mpi_real,my_neb(2),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
       call mpi_isend(buf3,ic,mpi_real,my_neb(4),mype &
                     ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch23
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch223(arr1,ll1,arr2,ll2,arr3,ll3,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  exchange haloes for two 2-d arrays and one 3-d array.
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 2-d array (=1)
,ll2 &                  ! vertical dimension of 2nd 2-d array (=1)
,ll3 &                  ! vertical dimension of 3-d array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout) :: &
 arr1 &                 ! 2-d array whose haloes are exchanged
,arr2                   ! 2-d array whose haloes are exchanged
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-d array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  north/south
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf2,ic,mpi_real,my_neb(1),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,jts+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf3,ic,mpi_real,my_neb(3),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  east/west
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf2,ic,mpi_real,my_neb(2),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
       call mpi_isend(buf3,ic,mpi_real,my_neb(4),mype &
                     ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf0(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf0(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j)=buf1(ic)
        enddo
        enddo
!
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j)=buf1(ic)
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch223
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine exch33(arr1,ll1,arr2,ll2,ihalo,jhalo)
!
!-----------------------------------------------------------------------
!
!***  exchange haloes two 3-d real arrays.
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 3-d array
,ll2 &                  ! vertical dimension of 2nd 3-d array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
 arr1                   ! 3-d array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-d array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  north/south
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j,k)
        enddo
        enddo
        enddo
!
        call mpi_isend(buf2,ic,mpi_real,my_neb(1),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j,k)
        enddo
        enddo
        enddo
!
        call mpi_isend(buf3,ic,mpi_real,my_neb(3),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  east/west
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        call mpi_isend(buf2,ic,mpi_real,my_neb(2),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
       call mpi_isend(buf3,ic,mpi_real,my_neb(4),mype &
                     ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
      endif
!   
!-----------------------------------------------------------------------
!***  store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
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
!***  exchange haloes three 3-d real arrays.
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 3-d array
,ll2 &                  ! vertical dimension of 2nd 3-d array
,ll3 &                  ! vertical dimension of 3rd 3-d array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
 arr1                   ! 3-d array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-d array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-d array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  north/south
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j,k)
        enddo
        enddo
        enddo
!
        call mpi_isend(buf2,ic,mpi_real,my_neb(1),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,jts+j,k)
        enddo
        enddo
        enddo
!
        call mpi_isend(buf3,ic,mpi_real,my_neb(3),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
      endif
!
!-----------------------------------------------------------------------
!***  store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  east/west
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        call mpi_isend(buf2,ic,mpi_real,my_neb(2),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        call mpi_isend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
      endif
!   
!-----------------------------------------------------------------------
!***  store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
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
!***  exchange haloes four 3-d real arrays.
!
!-----------------------------------------------------------------------
!***  argument variables
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
!
 ll1 &                  ! vertical dimension of 1st 3-d array
,ll2 &                  ! vertical dimension of 2nd 3-d array
,ll3 &                  ! vertical dimension of 3rd 3-d array
,ll4 &                  ! vertical dimension of 4th 3-d array
,ihalo &                ! number of halo rows in i direction
,jhalo                  ! number of halo rows in j direction
!
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll1),intent(inout) :: &
 arr1                   ! 3-d array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll2),intent(inout) :: &
 arr2                   ! 3-d array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll3),intent(inout) :: &
 arr3                   ! 3-d array whose haloes are exchanged
      real(kind=kfpt),dimension(ims:ime,jms:jme,ll4),intent(inout) :: &
 arr4                   ! 3-d array whose haloes are exchanged
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i,ibeg,ic,iend,ierr,irecv,isend,j,jbeg,jend,k
!
      integer(kind=kint),dimension(mpi_status_size) :: istat
!
      integer(kind=kint),dimension(4) :: ihandle
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  north/south
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(1),my_neb(1) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(3),my_neb(3) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to north
!-----------------------------------------------------------------------
!     
      ibeg=max(its-ihalo,ids)
      iend=min(ite+ihalo,ide)
!     
      if(my_neb(1)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,jte-j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr4(i,jte-j,k)
        enddo
        enddo
        enddo
!
        call mpi_isend(buf2,ic,mpi_real,my_neb(1),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to south
!-----------------------------------------------------------------------
!    
      if(my_neb(3)>=0)then
        ic=0
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,jts+j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr4(i,jts+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf3,ic,mpi_real,my_neb(3),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store results from south
!-----------------------------------------------------------------------
!
      if(my_neb(3)>=0)then
        ic=0 
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr4(i,jts-j-1,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  store from north
!-----------------------------------------------------------------------
!
      if(my_neb(1)>=0)then
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr1(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr2(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr3(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=0,jhalo-1
        do i=ibeg,iend
          ic=ic+1
          arr4(i,jte+j+1,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(1)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
      if(my_neb(3)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***
!***  east/west
!***
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  receive from west
!-----------------------------------------------------------------------
!    
      if(my_neb(4)>=0)then
        call mpi_irecv(buf0,ibufexch,mpi_real,my_neb(4),my_neb(4) &
                      ,mpi_comm_comp,ihandle(1),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  receive from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        call mpi_irecv(buf1,ibufexch,mpi_real,my_neb(2),my_neb(2) &
                      ,mpi_comm_comp,ihandle(2),irecv)
      endif
!
!-----------------------------------------------------------------------
!***  send to east
!-----------------------------------------------------------------------
!      
      jbeg=max(jts-jhalo,jds)
      jend=min(jte+jhalo,jde)
!      
      if(my_neb(2)>=0)then
        ibeg=ite-ihalo+1
        iend=ite
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf2(ic)=arr4(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf2,ic,mpi_real,my_neb(2),mype &
                      ,mpi_comm_comp,ihandle(3),isend)
      endif
!
!-----------------------------------------------------------------------
!***  send to west
!-----------------------------------------------------------------------
!       
      if(my_neb(4)>=0)then
        ibeg=its
        iend=its+ihalo-1
        ic=0
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr1(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr2(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr3(i,j,k)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          buf3(ic)=arr4(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(buf3,ic,mpi_real,my_neb(4),mype &
                      ,mpi_comm_comp,ihandle(4),isend)
      endif
!
!-----------------------------------------------------------------------
!***  store from west
!-----------------------------------------------------------------------
!
      if(my_neb(4)>=0)then
        ibeg=its-ihalo
        iend=its-1
        ic=0
        call mpi_wait(ihandle(1),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr4(i,j,k)=buf0(ic)
        enddo
        enddo
        enddo
      endif
!   
!-----------------------------------------------------------------------
!***  store from east
!-----------------------------------------------------------------------
!    
      if(my_neb(2)>=0)then
        ibeg=ite+1
        iend=ite+ihalo
        ic=0
        call mpi_wait(ihandle(2),istat,ierr)
!
        do k=1,ll1
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr1(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll2
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr2(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll3
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr3(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
!
        do k=1,ll4
        do j=jbeg,jend
        do i=ibeg,iend
          ic=ic+1
          arr4(i,j,k)=buf1(ic)
        enddo
        enddo
        enddo
      endif
!
      if(my_neb(4)>=0)then
        call mpi_wait(ihandle(4),istat,ierr)
      endif
!
      if(my_neb(2)>=0)then
        call mpi_wait(ihandle(3),istat,ierr)
      endif
!
!-----------------------------------------------------------------------
!
      end subroutine exch3333
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!--------------------------------------------------------------------
!
      end module module_exchange_gfs
!
!-----------------------------------------------------------------------
