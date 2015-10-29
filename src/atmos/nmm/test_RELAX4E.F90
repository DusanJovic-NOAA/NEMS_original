program test_relax4e
  use module_relax4e

  integer, parameter :: ni=39,nj=21
  integer :: i,j, mask(ni,nj),ibar,jbar
  real :: data(ni,nj),old(ni,nj),maxd,mind
  integer :: idata(ni,nj)

  ibar=(ni+1)/2
  jbar=(nj+1)/2

  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,nj
     do i=1,ni
        data(i,j)=0
        mask(i,j)=1
     enddo
  enddo
  !$OMP END PARALLEL DO

  do j=1,nj
     data(ibar,j)=real(nj/2-abs(j-jbar))/nj*2
     mask(ibar,j)=0
  enddo

  do i=1,ni
     data(i,jbar)=real(ni/2-abs(i-ibar))/ni*2
     mask(i,jbar)=0
  enddo

  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=1,nj
     do i=1,ni
        old(i,j)=data(i,j)
     enddo
  enddo

  call printit('Before relax:')

  call relax4e(data,mask,0.3678,ni*nj/4,&
       1,ni,1,nj, 1,ni,1,nj, 1,ni,1,nj)

  call printit('After relax:')

  call relax4e(data,mask,0.3678,10000,&
       1,ni,1,nj, 1,ni,1,nj, 1,ni,1,nj)

  call printit('Very relaxed:')

19 format('ERROR: boundary value changed: old(',I0,',',I0,')=',&
        F0.3,' new = ',F0.3)

  do j=1,nj
     if(data(ibar,j)/=old(ibar,j)) then
        write(0,19) ibar,j,old(ibar,j),data(ibar,j)
        stop 9
     endif
  enddo

  do i=1,ni
     if(data(i,jbar)/=old(i,jbar)) then
        write(0,19) i,jbar,old(i,jbar),data(i,jbar)
        stop 9
     endif
  enddo

contains

  subroutine printit(message)
    character*(*), intent(in) :: message

    maxd=data(1,1)
    mind=maxd
    !$OMP PARALLEL DO PRIVATE(i,j) REDUCTION(max:maxd) REDUCTION(min:mind)
    do j=1,nj
       do i=1,ni
          maxd=max(maxd,data(i,j))
          mind=min(mind,data(i,j))
       enddo
    enddo
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(i,j)
    do j=1,nj
       do i=1,ni
          idata(i,j)=max(0,min(100,nint((data(i,j)-mind)/(maxd-mind)*100)))
       enddo
    enddo
    !$OMP END PARALLEL DO

    print '(A)', message
39 format (39(I3,' '))
    do j=1,nj
       print 39,idata(:,j)
    enddo
  end subroutine printit
end program test_relax4e
