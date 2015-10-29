module module_TEST_GRID
  implicit none
  private
  public :: sis_init

  type, public :: solver_internal_state
     integer :: ims,ime,jms,jme,kms,kme,&
          ids,ide,jds,jde,kds,kde,&
          its,ite,jts,jte,kts,kte,&
          MYPE,MPI_COMM_COMP,NUM_PES
     integer, pointer :: i2d(:,:), i3d(:,:,:)
     real, pointer :: r2d(:,:), r3d(:,:,:)
   contains
     procedure init => sis_init
     procedure debug => sis_debug
     procedure message => sis_message
     procedure error => sis_error
     procedure fail_if => sis_fail_if
  end type solver_internal_state

contains

  subroutine sis_debug(sis,message)
    class(solver_internal_state), intent(in) :: sis
    character*(*), intent(in) :: message
101 format('NMM',I5,': ',A)
    print 101,sis%MYPE,message
  end subroutine sis_debug

  subroutine sis_message(sis,message)
    class(solver_internal_state), intent(in) :: sis
    character*(*), intent(in) :: message
101 format('NMM',I5,': ',A)
    if(sis%MYPE==0) &
       print 101,sis%MYPE,message
  end subroutine sis_message

  subroutine sis_error(sis,message)
    class(solver_internal_state), intent(in) :: sis
    character*(*), intent(in) :: message
101 format('NMM',I5,': ',A)
    print 101,sis%MYPE,message
    write(0,101) sis%MYPE,message
  end subroutine sis_error

  subroutine sis_fail_if(sis,failure_message,success_message,condition)
    use mpi
    class(solver_internal_state), intent(in) :: sis
    character*(*), intent(in) :: failure_message,success_message
    logical, intent(in) :: condition
    integer ierr
    if(condition) then
       call sis%error(failure_message)
       flush 6
       flush 0
       call MPI_Abort(MPI_COMM_WORLD,2,ierr)
       stop 2
    elseif(sis%MYPE==0) then
       call sis%message(success_message)
    endif
  end subroutine sis_fail_if

  subroutine sis_init(sis,ni,nj,nk,comm)
    use mpi
    class(solver_internal_state), intent(inout) :: sis
    integer, intent(in) :: ni,nj,nk,comm

    integer :: ierr,iranks,jranks,imax,jmax,i,j,irank,jrank,i1,i2,j1,j2
    real :: best_badness,badness

    sis%MPI_COMM_COMP=comm
    sis%MYPE=0
    sis%NUM_PES=0
    call MPI_Comm_rank(comm,sis%MYPE,ierr)
    call MPI_Comm_size(comm,sis%NUM_PES,ierr)

    iranks=sis%NUM_PES
    jranks=1
    best_badness = abs( (ni/iranks) - (nj/jranks) )
    do i=sis%NUM_PES,1,-1
       jmax=nint(real(sis%NUM_PES)/i)
       do j=1,jmax
          if(j*i==sis%NUM_PES) then
48           format(I0,'=',I0,'*',I0,': ',A)
             badness=abs(ni/i - nj/j)
             if(badness<best_badness) then
                !print 48,sis%NUM_PES,i,j,'best so far'
                iranks=i
                jranks=j
                best_badness=badness
             else
                !print 48,sis%NUM_PES,i,j,'not as good'
             endif
          endif
       enddo
    enddo

    if(ni<iranks*2 .or. nj<jranks*2 .or. nk<5) then
51     format('Bad ranks: grid=',I0,'x',I0,'x',I0,' ranks=',I0,'x',I0)
       write(0,51) ni,nj,nk,iranks,jranks
       call MPI_Abort(MPI_COMM_WORLD,2,ierr)
    end if

    jrank=sis%MYPE/iranks
    irank=mod(sis%MYPE,iranks)

    i1=nint(real(ni) * irank/iranks) + 1
    i2=nint(real(ni) * (irank+1)/iranks)

    j1=nint(real(nj) * jrank/jranks) + 1
    j2=nint(real(nj) * (jrank+1)/jranks)

11  format('Rank ',I02,' is at irank=',I02,'/',I0,' jrank=',I02,'/',I0,&
         ' i=',I0,':',I0,' j=',I0,':',I0)
    print 11,sis%MYPE,irank,iranks,jrank,jranks,i1,i2,j1,j2

  ! MEMORY DIMS  ; MPI DIMS   ; DOMAIN DIMS
    sis%ims=i1-3 ; sis%its=i1 ; sis%ids=1
    sis%ime=i2+3 ; sis%ite=i2 ; sis%ide=ni
    sis%jms=j1-3 ; sis%jts=j1 ; sis%jds=1
    sis%jme=j2+3 ; sis%jte=j2 ; sis%jde=nj
    sis%kms=-2   ; sis%kts=1  ; sis%kds=1
    sis%kme=nk+3 ; sis%kte=nk ; sis%kde=nk

43  format('Rank ',I0,': allocate sis%i3d(',I0,':',I0,', ',I0,':',I0,', ',I0,':',I0,')')
    print 43,sis%mype,sis%ims,sis%ime,sis%jms,sis%jme,sis%kms,sis%kme

    flush 5
    flush 6
    flush 0

    allocate(sis%i3d(sis%ims:sis%ime,sis%jms:sis%jme,sis%kms:sis%kme))
    allocate(sis%r3d(sis%ims:sis%ime,sis%jms:sis%jme,sis%kms:sis%kme))

    allocate(sis%i2d(sis%ims:sis%ime,sis%jms:sis%jme))
    allocate(sis%r2d(sis%ims:sis%ime,sis%jms:sis%jme))

  end subroutine sis_init

end module module_TEST_GRID
