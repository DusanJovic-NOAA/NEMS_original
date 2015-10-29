program test
  use mpi
  use module_TEST_GRID
  use module_REDUCTION
  implicit none

  type(solver_internal_state) :: sis
  integer, parameter :: ni=200,nj=50,nk=150
  integer, parameter :: imin=7,jmin=17,kmin=100
  integer, parameter :: imax=130,jmax=11,kmax=5
  integer, parameter :: ibar=ni/2,jbar=nj/2

  integer :: i1,i2,j1,j2,i,j,k,ifnd,kms,kme
  real :: rfnd,rout
  integer :: iloc,jloc,kloc,rankloc,ierr,iout
  double precision :: dout

  character(len=80) :: success,failure

  call MPI_Init(ierr)
  call sis_init(sis,ni,nj,nk,MPI_COMM_WORLD)
  i1=sis%its
  i2=sis%ite
  j1=sis%jts
  j2=sis%jte

  sis%r3d=9e19
  sis%r2d=9e19
  sis%i2d=-20000
  sis%i3d=-20000
  do k=1,nk
     !$OMP PARALLEL DO PRIVATE(i,j)
     do j=j1,j2
        do i=i1,i2
           sis%r3d(i,j,k) = real(k-1)*ni*nj + real(j-1)*ni + real(i-1) + 1
           sis%i3d(i,j,k) = i-1+ni*(j-1+nj*(k-1)) +1
        enddo
     enddo
     !$OMP END PARALLEL DO
  enddo
  !$OMP PARALLEL DO PRIVATE(i,j)
  do j=j1,j2
     do i=i1,i2
        sis%r2d(i,j) = real(j-1)*ni + real(i-1) + 1
        sis%i2d(i,j) = i-1+ni*(j-1) +1
     enddo
  enddo
  !$OMP END PARALLEL DO

  if( imax>=i1 .and. imax<=i2 .and. jmax>=j1 .and. jmax<=j2) then
14   format('Setting i=',I0,' j=',I0,' k=',I0,' to 777777777')
     print 14,imax,jmax,kmax
     sis%r3d(imax,jmax,kmax)=777777777
     sis%i3d(imax,jmax,kmax)=777777777
     sis%r2d(imax,jmax)=777777777
     sis%i2d(imax,jmax)=777777777
  endif

  if( imin>=i1 .and. imin<=i2 .and. jmin>=j1 .and. jmin<=j2) then
15   format('Setting i=',I0,' j=',I0,' k=',I0,' to -888888888')
     print 15,imin,jmin,kmin
     sis%r3d(imin,jmin,kmin)=-888888888
     sis%i3d(imin,jmin,kmin)=-888888888
     sis%r2d(imin,jmin)=-888888888
     sis%i2d(imin,jmin)=-888888888
  endif

  ! FIND 3D TESTS

  call find(sis,rfnd,iloc,jloc,kloc,rankloc,sis%r3d,&
            FIND_MAX,sis%ims,sis%ime,sis%jms,sis%jme,sis%kms,sis%kme)
38 format('R3D max found at ',I0,',',I0,',',I0,' rank ',I0,' value ',F0.7)
  !print 38,iloc,jloc,kloc,rankloc,rfnd
  call fail_if(sis,'Wrong R3D max location.','Correct R3D max location.',&
       iloc/=imax .or. jloc/=jmax .or. kloc/=kmax)

  call find(sis,ifnd,iloc,jloc,kloc,rankloc,sis%i3d,&
            FIND_MAX,sis%ims,sis%ime,sis%jms,sis%jme,sis%kms,sis%kme)
33 format('I3D max found at ',I0,',',I0,',',I0,' rank ',I0,' value ',I0)
  !print 33,iloc,jloc,kloc,rankloc,ifnd
  call fail_if(sis,'Wrong I3D max location.','Correct I3D max location.',&
       iloc/=imax .or. jloc/=jmax .or. kloc/=kmax)

  call find(sis,rfnd,iloc,jloc,kloc,rankloc,sis%r3d,&
            FIND_MIN,sis%ims,sis%ime,sis%jms,sis%jme,sis%kms,sis%kme)
58 format('R3D min found at ',I0,',',I0,',',I0,' rank ',I0,' value ',F0.7)
  !print 58,iloc,jloc,kloc,rankloc,rfnd
  call fail_if(sis,'Wrong R3D min location','Correct R3D min location.',&
       iloc/=imin .or. jloc/=jmin .or. kloc/=kmin)

  call find(sis,ifnd,iloc,jloc,kloc,rankloc,sis%i3d,&
            FIND_MIN,sis%ims,sis%ime,sis%jms,sis%jme,sis%kms,sis%kme)
53 format('I3D min found at ',I0,',',I0,',',I0,' rank ',I0,' value ',I0)
  !print 53,iloc,jloc,kloc,rankloc,ifnd
  call fail_if(sis,'Wrong I3D min location','Correct I3D min location.',&
       iloc/=imin .or. jloc/=jmin .or. kloc/=kmin)

  ! FIND 2D TESTS

  call find(sis,rfnd,iloc,jloc,rankloc,sis%r2d,&
            FIND_MAX,sis%ims,sis%ime,sis%jms,sis%jme)
238 format('R2D max found at ',I0,',',I0,' rank ',I0,' value ',F0.7)
  !print 238,iloc,jloc,rankloc,rfnd
  call fail_if(sis,'Wrong R2D max location.','Correct R2D max location.',&
       iloc/=imax .or. jloc/=jmax)

  call find(sis,ifnd,iloc,jloc,rankloc,sis%i2d,&
            FIND_MAX,sis%ims,sis%ime,sis%jms,sis%jme)
233 format('I2D max found at ',I0,',',I0,' rank ',I0,' value ',I0)
  !print 233,iloc,jloc,rankloc,ifnd
  call fail_if(sis,'Wrong I2D max location.','Correct I2D max location.',&
       iloc/=imax .or. jloc/=jmax)

  call find(sis,rfnd,iloc,jloc,rankloc,sis%r2d,&
            FIND_MIN,sis%ims,sis%ime,sis%jms,sis%jme)
258 format('R2D min found at ',I0,',',I0,' rank ',I0,' value ',F0.7)
  !print 258,iloc,jloc,rankloc,rfnd
  call fail_if(sis,'Wrong R2D min location','Correct R2D min location.',&
       iloc/=imin .or. jloc/=jmin)

  call find(sis,ifnd,iloc,jloc,rankloc,sis%i2d,&
            FIND_MIN,sis%ims,sis%ime,sis%jms,sis%jme)
253 format('I2D min found at ',I0,',',I0,' rank ',I0,' value ',I0)
  !print 253,iloc,jloc,rankloc,ifnd
  call fail_if(sis,'Wrong I2D min location','Correct I2D min location.',&
       iloc/=imin .or. jloc/=jmin)

  ! Reduce tests 
#define REDUCE_3D_TEST(OUT,VAR,METHOD,EXPECT,WHAT) \
  call reduce(sis,OUT,sis%VAR,METHOD,sis%ims,sis%ime,sis%jms,sis%jme,sis%kms,sis%kme) ; \
  write(success,'(A," test succeeded.")') WHAT ; \
  write(failure,'(A," test failed.")') WHAT ; \
  call fail_if(sis,failure,success,.not. EXPECT)

#define REDUCE_2D_TEST(OUT,VAR,METHOD,EXPECT,WHAT) \
  call reduce(sis,OUT,sis%VAR,METHOD,sis%ims,sis%ime,sis%jms,sis%jme) ; \
  write(success,'(A," test succeeded.")') WHAT ; \
  write(failure,'(A," test failed.")') WHAT ; \
  call fail_if(sis,failure,success,.not. EXPECT)

  REDUCE_3D_TEST(rout,r3d,REDUCE_MAX,abs(rout-777777777.)<10000,'R3D max reduction')
  REDUCE_3D_TEST(rout,r3d,REDUCE_MIN,abs(rout+888888888.)<10000,'R3D min reduction')
  REDUCE_3D_TEST(dout,r3d,REDUCE_MAX,abs(dout-777777777.)<10000,'R3D dbl max reduction')
  REDUCE_3D_TEST(dout,r3d,REDUCE_MIN,abs(dout+888888888.)<10000,'R3D dbl min reduction')

  REDUCE_2D_TEST(rout,r2d,REDUCE_MAX,abs(rout-777777777.)<10000,'R2D max reduction')
  REDUCE_2D_TEST(rout,r2d,REDUCE_MIN,abs(rout+888888888.)<10000,'R2D min reduction')
  REDUCE_2D_TEST(dout,r2d,REDUCE_MAX,abs(dout-777777777.)<10000,'R2D dbl max reduction')
  REDUCE_2D_TEST(dout,r2d,REDUCE_MIN,abs(dout+888888888.)<10000,'R2D dbl min reduction')

  sis%r3d=1
  if( imin>=i1 .and. imin<=i2 .and. jmin>=j1 .and. jmin<=j2) then
     sis%r3d(imin,jmin,kmin)=2
  endif

  REDUCE_3D_TEST(rout,r3d,REDUCE_ADD,rout/=iloc*jloc*kloc+1,'R3D add')
  REDUCE_3D_TEST(rout,r3d,REDUCE_ADD,rout/=iloc*jloc*kloc+1,'R3D add')
  REDUCE_3D_TEST(dout,r3d,REDUCE_ADD_DOUBLE,dout/=iloc*jloc*kloc+1,'R3D add_double')
  REDUCE_3D_TEST(dout,r3d,REDUCE_ADD_DOUBLE,dout/=iloc*jloc*kloc+1,'R3D add_double')

  REDUCE_2D_TEST(rout,r2d,REDUCE_ADD,rout/=iloc*jloc+1,'R2D add')
  REDUCE_2D_TEST(rout,r2d,REDUCE_ADD,rout/=iloc*jloc+1,'R2D add')
  REDUCE_2D_TEST(dout,r2d,REDUCE_ADD_DOUBLE,dout/=iloc*jloc+1,'R2D add_double')
  REDUCE_2D_TEST(dout,r2d,REDUCE_ADD_DOUBLE,dout/=iloc*jloc+1,'R2D add_double')

  call MPI_Finalize(ierr)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine fail_if(sis,fail_message,success_message,condition)
    type(solver_internal_state), intent(in) :: sis
    character*(*), intent(in) :: fail_message,success_message
    logical, intent(in) :: condition
    integer ierr
    if(condition) then
       write(0,*) trim(fail_message)
       print *,trim(fail_message)
       FLUSH 5
       FLUSH 6
       call MPI_Abort(MPI_COMM_WORLD,2,ierr)
       stop 2
    else if(sis%MYPE==0) then
       print *,trim(success_message)
    end if
  end subroutine fail_if
end program test
