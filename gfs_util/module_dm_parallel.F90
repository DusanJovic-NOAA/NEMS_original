!-----------------------------------------------------------------------
                        module module_dm_parallel_gfs
!-----------------------------------------------------------------------
!
!***  this module contains all codes directly related to distributed
!***  memory issues except for halo exchange although note that the
!***  halo widths must be set here.
!
!-----------------------------------------------------------------------
!
      use module_include
      use esmf_mod
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!---domain decomposition info-------------------------------------------
!-----------------------------------------------------------------------
      integer :: isee=055,jsee=138,lsee=10,pesee=06
!
integer(kind=kint),parameter :: &
!
 ihalo=3 &                   ! halo width in i direction
,jhalo=3 &                   ! halo width in j direction
!
!***  hardwire to 100 the maximum number of server groups allowed.
!***  this is far greater than should ever be needed.
!
,max_groups=100 &            ! max number of quilt server groups
!
!***  for now, set the number of threads here to 1.
!***  clearly, this must be reconciled before actually using threading.
!
,num_tiles=1                 ! number of threads
!
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 ide &                       ! ending data index, x direction
,ids &                       ! starting data index, x direction
,ims &                       ! the starting memory i for each task
,ime &                       ! the ending memory i for each task
,its &                       ! the starting integration i for each task
,ite &                       ! the ending integration i for each task
,iquilt_group &              ! number of i/o server groups
,jde &                       ! ending data index, y direction
,jds &                       ! starting data index, y direction
,jms &                       ! the starting memory j for each task
,jme &                       ! the ending memory j for each task
,jts &                       ! the starting integration j for each task
,jte &                       ! the ending integration j for each task
,mpi_comm_comp &             ! local mpi communicator
,mpi_comm_inter &            ! intercommunicator for the i/o servers
,mype_share &                ! my task id to be seen by any uses
,npes &                      ! total number of forecast tasks
,num_pts_max                 ! max points in any task's subdomain
!
integer(kind=kint) :: &
 its_b1 &                    ! its and 1 point from global boundary
,its_b2 &                    ! its and 2 points from global boundary
,its_h1 &                    ! its and 1 point into halo
,its_h2 &                    ! its and 2 points into halo
,its_b1_h1 &                 ! its and _b1 and _h1
,its_b1_h2 &                 ! its and _b1 and _h2
,ite_b1 &                    ! ite and 1 point from global boundary
,ite_b2 &                    ! ite and 2 points from global boundary
,ite_h1 &                    ! ite and 1 point into halo
,ite_h2 &                    ! ite and 2 points into halo
,ite_b1_h1 &                 ! ite and _b1 and _h1
,ite_b1_h2 &                 ! ite and _b1 and _h2
,jts_b1 &                    ! jts and 1 point from global boundary
,jts_b2 &                    ! jts and 2 points from global boundary
,jts_h1 &                    ! jts and 1 point into halo
,jts_h2 &                    ! jts and 2 points into halo
,jts_b1_h1 &                 ! jts and _b1 and _h1
,jts_b1_h2 &                 ! jts and _b1 and _h2
,jte_b1 &                    ! jte and 1 point from global boundary
,jte_b2 &                    ! jte and 2 points from global boundary
,jte_h1 &                    ! jte and 1 point into halo
,jte_h2 &                    ! jte and 2 points into halo
,jte_b1_h1 &                 ! jte and _b1 and _h1
,jte_b1_h2                   ! jte and _b1 and _h2
!
integer(kind=kint) :: &
 ide_m1 & 
,ide_m2 & 
,ids_p1 & 
,jde_m1 & 
,jde_m2 & 
,jds_p1 
integer(kind=kint),dimension(8) :: &
 my_neb                      ! my task's eight neighbors 
!
integer(kind=kint),dimension(max_groups) :: &
 mpi_comm_inter_array &      ! intercommunicators for the integration tasks
,num_serv_per_grp            ! number of servers in each group
!
integer(kind=kint),allocatable,dimension(:) :: &
 local_iend &
,local_istart &
,local_jend &
,local_jstart
!
!-----------------------------------------------------------------------
!
       contains
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine setup_servers(mype,inpes,jnpes,npes,mpi_intra)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  setup_servers splits the communicator between integration 
!***  and output tasks (servers).
!
!-----------------------------------------------------------------------
!
!   input argument list:
!         mype - my task id
!        inpes - number of mpi tasks in the x direction
!        jnpes - number of mpi tasks in the y direction
!         npes - total number of mpi tasks provided to the job
!                npes at least equal the product of inpes*jnpes otherwise the
!                integration cannot proceed. the difference between the product
!                npes_mod and npes is the number of mpi tasks that are available
!                for i/o serving. this can be zero, in which case output will
!                write a direct access file that can be separately quilted. 
!                in order to skip the separate quilting step, make sure that
!                the number of mpi tasks that the code is initiated with is at
!                least one greater than npes_mod.
!                later in the routine, npes is reset to the number of fcst tasks.
!    mpi_intra - the global communicator.
!
!   input files:
!         none but the code does attempt to read the environment variable "server_groups".
!         this is the number of independent groups of server tasks. the default is one
!         and should be ok for most applications of the code. if one set of i/o
!         servers can not complete before the next ouput time then additional i/o server
!         groups would be useful.
!
!-----------------------------------------------------------------------
!***  argument variables.
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: &
 mype &                     ! each task's id
,inpes &                    ! number of compute tasks in x direction
,jnpes &                    ! number of compute tasks in y direction
,mpi_intra                  ! global communicator
!
      integer(kind=kint),intent(inout) :: &
 npes                       ! total number of tasks provided then
                            ! converted to the number of fcst tasks
!-----------------------------------------------------------------------
!***  local variables.
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: comdup,i,icc,icolor,iendq,iendxx,ierr &
                           ,igroup,igroup_x,iqserver,irlr,iss,issl &
                           ,istaq,istaxx,iworld,iworld_minus,ixx,jj,kk & 
                           ,npes_mod,one
!
      integer(kind=kint),allocatable,dimension(:) :: irank
!
      logical :: yes
!
      character(4) :: get
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
!***  let npes_mod be the product of inpes and jnpes (namelist variables).
!***  this is the number of mpi tasks the executable has been built for.
!***  npes, returned from mpi_comm_size, must be at least this size
!***  otherwise the integration cannot proceed. the difference between
!***  npes_mod and npes is the number of mpi tasks that are available
!***  for i/o serving. this can be zero, in which case output will
!***  write a direct access file that can be separately quilted.
!***  in order to skip the separate quilting step, make sure that
!***  the number of mpi tasks that the code is initiated with is at
!***  least one greater than npes_mod.
!
!-----------------------------------------------------------------------
!
!!!   mype=mype_share
      npes_mod=inpes*jnpes
      mpi_comm_comp=mpi_intra
!
!-----------------------------------------------------------------------
!
!***  at this point npes is the total number of mpi tasks. we will
!***  reset this at the end of the subroutine to the number of mpi
!***  tasks that are working on the model integration.
!
!***  first, however, we need to make sure that a sufficient number
!***  of mpi tasks have been initiated. if not, we will stop.
!
!-----------------------------------------------------------------------
      if(npes<npes_mod)then
         write(0,*)' ***********************************************'
         write(0,*)' ***********************************************'
         write(0,*)' *************major problem*********************'
         write(0,*)' *************major problem*********************'
         write(0,*)' *************major problem*********************'
         write(0,*)' *************major problem*********************'
         write(0,*)' '
         write(0,*)' there are insufficient mpi tasks to continue'
         write(0,*)' you must specify at least ',npes_mod,' tasks'
         write(0,*)' stopping now'
         write(0,*)' '
         write(0,*)' *************major problem*********************'
         write(0,*)' *************major problem*********************'
         write(0,*)' *************major problem*********************'
         write(0,*)' *************major problem*********************'
         write(0,*)' ***********************************************'
         write(0,*)' ***********************************************'
!!!      call mpi_abort(mpi_comm_world,1,ierr)
         call mpi_abort(mpi_intra,1,ierr)
      endif
!-----------------------------------------------------------------------
!
!***  ok, we have a sufficient number of mpi tasks to continue.
!
!***  how many groups of servers?  the default is 1 group.
!***  the environment variable, server_groups, can be used to
!***  specify more server groups if desired.
!
!-----------------------------------------------------------------------
      get='1'
      call getenv('server_groups',get)
      read(get,fmt='(i4)')iquilt_group
!     iquilt_group=max(iquilt_group,1)
      one=1
      iquilt_group=max(iquilt_group,one)
!-----------------------------------------------------------------------
!
!***  error check for number of groups - maximum is 100 - this is alot!
!
!-----------------------------------------------------------------------
      if(iquilt_group>100)then
        write(0,*)' ***** iquilt_group is greater than 100'
        write(0,*)' ***** do you really want this ?'
        write(0,*)' ***** if so then increase size in mpp.h'
        write(0,*)' ***** also, change if check in setup_servers'
        write(0,*)' ***** resetting the number of server groups to 100'
        write(0,*)' ***** we are continuing ....   '
        iquilt_group=max_groups
      endif
!
      if(mype==0)then
        write(0,*)' number of server groups: ',iquilt_group
      endif
!-----------------------------------------------------------------------
!
!***  compute the number of quilt servers per group.
!***  all mpi tasks beyond npes_mod will be quilt servers.
!***  if the number of servers is not equally divisible by
!***  the number of groups of servers then some groups may have
!***  more servers then others.  this is fine.
!***  note that we require at least one server per group.
!***  we may need to reduce the number of server groups if
!***  it exceeds the number of servers.
!
!-----------------------------------------------------------------------
      iqserver=npes-npes_mod
!
      if(iqserver==0)then
        if(mype==0)then
          write(0,*)' *** you specified 0 i/o servers '
          write(0,*)' output will write a direct access file'
        endif
         iquilt_group = 0
      endif
!
      if(iquilt_group>iqserver)then
        iquilt_group=iqserver
        write(0,*)' ***** not enough servers'
        write(0,*)' ***** we need to reduce the number of server groups'
        write(0,*)' ***** number of server groups is ',iquilt_group
      endif
!
      do i=0,iquilt_group-1
        call para_range(one,iqserver,iquilt_group,i,istaq,iendq)
        num_serv_per_grp(i+1)=iendq-istaq+1
        if(mype==0)then
          write(0,*)' number of servers for group ',i+1,' is ',num_serv_per_grp(i+1)
        endif
      enddo
!
!-----------------------------------------------------------------------
!***  set up the "color" for mpi_comm_split.
!***  those tasks which will do model integration will be color 0.
!***  the quilt server tasks will have the color of the group number to
!***  which they will belong.
!-----------------------------------------------------------------------
!
      if(mype<npes_mod)then
        icolor=0
      else 
        istaxx=npes_mod
        do i=1,iquilt_group
          iendxx=istaxx+num_serv_per_grp(i)-1
          if(mype>=istaxx.and.mype<=iendxx)then
            icolor=i
          endif
          istaxx=iendxx+1
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  split the communicator - the new intracommunicator for all tasks
!***  is mpi_comm_comp. mpi_comm_world (mpi_intra) is still available but it 
!***  refers to all the mpi tasks (model integration and i/o serving).
!-----------------------------------------------------------------------
!        
!!!   call mpi_comm_dup(mpi_comm_world,comdup,ierr)
      call mpi_comm_dup(mpi_intra,comdup,ierr)
      call mpi_comm_split(comdup,icolor,mype,mpi_comm_comp,ierr)
!
!-----------------------------------------------------------------------
!***  at this point we have a new communicator, mpi_comm_comp,
!***  that can be used by the forecast tasks and the quilt server tasks
!***  for their internal communications. on to the intercommunicators ...
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  now we must create the intercommunicators for use between the mpi
!***  tasks doing the model integration and the mpi tasks for each 
!***  server group.  the first step is to exclude the tasks that do not
!***  belong.  we will do this for each server group by excluding the
!***  tasks from all of the other server groups.
!-----------------------------------------------------------------------
!
      allocate(irank(iqserver))
!
      ixx=npes_mod
!
!-----------------------------------------------------------------------
      inter_comm : do i=1,iquilt_group
!-----------------------------------------------------------------------
        yes=.true.
!
        if(mype<npes_mod)then
          irlr=ixx
        else
          irlr=0
        endif
!
        icc=0
        iss=npes_mod
!
!***  this is the first possible task id that could be excluded.
!
        do jj=1,iquilt_group
          if(jj/=i)then
            issl=iss
            do kk=1,num_serv_per_grp(jj)
              icc=icc+1
              irank(icc)=issl
              if(mype==issl)yes=.false.
              issl=issl+1
            enddo
          endif
          iss=iss+num_serv_per_grp(jj)
        enddo
!
!-----------------------------------------------------------------------
!***  at this point we have an array, irank, with task ids to exclude.
!***  there are icc of them.
!***  create a new group with the tasks from the other server groups
!***  excluded and then create a new communicator (iworld_minus) that
!***  contains only the mpi tasks doing the model integration and the
!***  tasks that belong to the server group we are considering.
!-----------------------------------------------------------------------
!
!!!   iworld=mpi_comm_world
      iworld=mpi_intra
      call mpi_comm_group(iworld,igroup,ierr)
      call mpi_group_excl(igroup,icc,irank,igroup_x,ierr)
      call mpi_comm_create(iworld,igroup_x,iworld_minus,ierr)
      call mpi_group_free(igroup,ierr)
      call mpi_group_free(igroup_x,ierr)
!
!-----------------------------------------------------------------------
!***  at this point we have a communicator that excludes the tasks we dont want.
!***  create an intercommunicator for use between the mpi tasks doing the model
!***  integration and the i/o server group we are considering. this process is
!***  a collective routine so it can only be done by the tasks that have not 
!***  been excluded. save this new communicator in mpi_comm_inter for use by
!***  the tasks that belong to the server group that we are considering. the
!***  tasks that are performing the model integration will reference
!***  mpi_comm_inter_array() since we will need to select which server
!***  group we wish to communicate with.
!-----------------------------------------------------------------------
!
      if(yes)then
        call mpi_intercomm_create(mpi_comm_comp,0,iworld_minus,irlr,0 &
                                 ,mpi_comm_inter_array(i),ierr)
         mpi_comm_inter=mpi_comm_inter_array(i)
      endif
!
!!!   call mpi_barrier(mpi_comm_world,ierr)
      call mpi_barrier(mpi_intra,ierr)
!
!-----------------------------------------------------------------------
      enddo  inter_comm
!-----------------------------------------------------------------------
!***
!***  set npes to the number of tasks working on the model integration.
!***
      npes=npes-iqserver
!
      if(mype==0)then
        write(0,*)' number of integration tasks: ',npes
        write(0,*)' total number of quilt servers: ',iqserver
        write(0,*)' exit setup_servers npes=',npes
      endif
!***
      deallocate (irank)
!-----------------------------------------------------------------------
!
      end subroutine setup_servers
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine para_range &
!
      (n1,n2,nprocs,irank,ista,iend)
!
!-----------------------------------------------------------------------
!
!   input argument list:
!     n1 - first interate value
!     n2 - last interate value
!     nprocs - number of mpi tasks
!     irank - my task id
!
!   output argument list:
!     ista - first loop value
!     iend - last loop value
!
!   output files:  none
!
!   subprograms called:  none
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) ::                                  &
        irank                                                           &
       ,n1                                                              &
       ,n2                                                              &
       ,nprocs
!
      integer(kind=kint),intent(out) ::                                 &
        iend                                                            &
       ,ista
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
      integer(kind=kint) ::                                             &
        iwork1                                                          &
       ,iwork2
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      iwork1=(n2-n1+1)/nprocs
      iwork2=mod (n2-n1+1,nprocs)
      ista=irank*iwork1+n1+min(irank,iwork2)
      iend=ista+iwork1-1
      if(iwork2>irank)iend=iend+1
!-----------------------------------------------------------------------
!
      end subroutine para_range
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine decomp &
!
      (inpes,jnpes,im,jm,lm,global,ijcount,grid_comp)
!
!-----------------------------------------------------------------------
!
!***  decomp specifies the domain decomposition.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  argument variables.
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: im &
                                      ,jm &
                                      ,lm &
                                      ,inpes &
                                      ,jnpes 
!
      logical,intent(in) :: global
!
      integer(kind=kint),dimension(2),intent(in) :: ijcount
!
      type(esmf_gridcomp),intent(in) :: grid_comp
!
!-----------------------------------------------------------------------
!***  local variables.
!-----------------------------------------------------------------------
!
      integer(kind=kint) :: i &
                           ,i_add &
                           ,icol &
                           ,iend &
                           ,ierr &
                           ,iguess &
                           ,ipe &
                           ,irecv &
                           ,iremain &
                           ,irtn &
                           ,isend &
                           ,istart &
                           ,istat &
                           ,j &
                           ,j_add  &
                           ,jend &
                           ,jguess &
                           ,jremain &
                           ,jrow &
                           ,jstart &
                           ,k2 &
                           ,l_remain &
                           ,lyr_frac &
                           ,my_e &
                           ,my_n &
                           ,my_ne &
                           ,my_nw &
                           ,my_s &
                           ,my_se &
                           ,my_sw &
                           ,my_w &
                           ,myi &
                           ,myj &
                           ,mype &
                           ,n &
                           ,npe &
                           ,num_pts 
!
      integer(kind=kint),dimension(4) :: limits
!
      integer(kind=kint),dimension(mpi_status_size) :: jstat
!
      integer(kind=kint),allocatable,dimension(:,:) :: ijcount_all      &
                                                      ,itemp
!
      type(esmf_vm) :: vm
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!***  retrieve mype from the esmf virtual machine
!
      call esmf_gridcompget(grid_comp,vm=vm,rc=irtn)
      call esmf_vmget(vm,localpet=mype,rc=irtn)
!
      mype_share=mype  ! the gather/scatter routines and other parallel
                       ! routines will get mype from here and thus
                       ! remain esmf-neutral themselves
!
      npes=inpes*jnpes
!
!-----------------------------------------------------------------------
!***
!***  compute the index limits within each mpi task.
!***  we will divide the number of points in the i and j directions
!***  evenly and then give as many of the initial tasks in each
!***  direction single additional points until any remainders are
!***  used up.
!***  the task ids will start with 0 in the lower left corner and
!***  increase in the i direction and then in the j direction.
!***
!-----------------------------------------------------------------------
!***  the full dimensions of the integration domain.
!-----------------------------------------------------------------------
!
      if(global)then
        ids=1
        ide=im+2
        jds=1
        jde=jm+2
      else
        ids=1
        ide=im
        jds=1
        jde=jm
      endif
!
!-----------------------------------------------------------------------
!***  find the remainders of points in each direction that will be
!***  incrementally added to each of the final tasks in each direction.
!-----------------------------------------------------------------------
!
      iguess=(ide-ids+1)/inpes
      iremain=(ide-ids+1)-iguess*inpes
      jguess=(jde-jds+1)/jnpes
      jremain=(jde-jds+1)-jguess*jnpes
!
!-----------------------------------------------------------------------
!***  let every task know where all other tasks start and end
!***  on the full grid.
!***  each task will save its own start/end values.
!-----------------------------------------------------------------------
!
      allocate(local_istart(0:npes-1),stat=istat)
      allocate(local_iend(0:npes-1),stat=istat)
      allocate(local_jstart(0:npes-1),stat=istat)
      allocate(local_jend(0:npes-1),stat=istat)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(mype==0)then
        allocate(ijcount_all(2,0:npes-1),stat=istat)
        ijcount_all(1,0)=ijcount(1)
        ijcount_all(2,0)=ijcount(2)
        local_istart(0)=ids
        local_iend(0)=ids+ijcount(1)-1
        local_jstart(0)=jds
        local_jend(0)=jds+ijcount(2)-1
!
        do npe=1,npes-1
          call mpi_recv(ijcount_all(1,npe),2,mpi_integer,npe,npe        &
                       ,mpi_comm_comp,jstat,irecv)
!
!-----------------------------------------------------------------------
!***  find the start and end i values on this task's 
!***  primary integration region (inside the haloes).
!-----------------------------------------------------------------------
!
          icol=mod(npe,inpes)+1
          if(icol==1)then
            local_istart(npe)=ids
          else
            local_istart(npe)=local_istart(npe-1)+ijcount_all(1,npe-1)
          endif
          local_iend(npe)=local_istart(npe)+ijcount_all(1,npe)-1
!
!-----------------------------------------------------------------------
!***  find the start and end j values on this task's 
!***  primary integration region (inside the haloes).
!-----------------------------------------------------------------------
!
          jrow=npe/inpes+1
          if(jrow==1)then
            local_jstart(npe)=jds
          else
            local_jstart(npe)=local_jstart(npe-inpes)+ijcount_all(2,npe-inpes)
          endif
          local_jend(npe)=local_jstart(npe)+ijcount_all(2,npe)-1
        enddo
!
      else
        call mpi_send(ijcount,2,mpi_integer,0,mype,mpi_comm_comp,isend)
      endif
!
      call mpi_bcast(local_istart,npes,mpi_integer,0,mpi_comm_comp,ierr)
      call mpi_bcast(local_iend,npes,mpi_integer,0,mpi_comm_comp,ierr)
      call mpi_bcast(local_jstart,npes,mpi_integer,0,mpi_comm_comp,ierr)
      call mpi_bcast(local_jend,npes,mpi_integer,0,mpi_comm_comp,ierr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!
      local_ij: do npe=0,npes-1
!
        if(mype==npe)then
          its=local_istart(npe)
          ite=local_iend(npe)
          jts=local_jstart(npe)
          jte=local_jend(npe)
        endif
!
      enddo local_ij
!
!-----------------------------------------------------------------------
!***  the memory or storage dimensions include the halo.
!-----------------------------------------------------------------------
!
!     ims=max(its-ihalo,ids)
!     ime=min(ite+ihalo,ide)
!     jms=max(jts-jhalo,jds)
!     jme=min(jte+jhalo,jde)
      ims=its-ihalo
      ime=ite+ihalo
      jms=jts-jhalo
      jme=jte+jhalo
!
!-----------------------------------------------------------------------
!***  additional loop limits regarding global boundary and haloes.
!***  if "_bn" is appended to a start/end limit then that means
!***  that the mpi tasks along the global boundary must stay n
!***  points away from that boundary.
!***  if "_hn" is appended to a start/end limit then that means
!***  that each task will compute n points into the halo unless
!***  being blocked by the global boundary.
!-----------------------------------------------------------------------
!
      ids_p1=max(its,ids+1)
      ide_m1=min(ite,ide-1)
      ide_m2=min(ite,ide-2)
      jds_p1=max(jts,jds+1)
      jde_m1=min(jte,jde-1)
      jde_m2=min(jte,jde-2)
!
      its_b1=max(its,ids+1)
      ite_b1=min(ite,ide-1)
      its_b2=max(its,ids+2)
      ite_b2=min(ite,ide-2)
      jts_b1=max(jts,jds+1)
      jte_b1=min(jte,jde-1)
      jts_b2=max(jts,jds+2)
      jte_b2=min(jte,jde-2)
!
      its_h1=max(its-1,ids)
      ite_h1=min(ite+1,ide)
      its_h2=max(its-2,ids)
      ite_h2=min(ite+2,ide)
      jts_h1=max(jts-1,jds)
      jte_h1=min(jte+1,jde)
      jts_h2=max(jts-2,jds)
      jte_h2=min(jte+2,jde)
!
      its_b1_h1=max(its-1,ids+1)
      ite_b1_h1=min(ite+1,ide-1)
      ite_b1_h2=min(ite+2,ide-1)
      jts_b1_h1=max(jts-1,jds+1)
      jte_b1_h1=min(jte+1,jde-1)
      jte_b1_h2=min(jte+2,jde-1)
!
      if(mype==0)then
        write(0,*)' ids=',ids,' ide=',ide,' jds=',jds,' jde=',jde
      endif
!
      do npe=0,npes-1
        if(mype==npe)then
          write(0,*)' pe=',mype
          write(0,*)' its=',its,' ite=',ite,' jts=',jts,' jte=',jte
          write(0,*)' ims=',ims,' ime=',ime,' jms=',jms,' jme=',jme
          write(0,*)' ids=',ids,' ide=',ide,' jds=',jds,' jde=',jde
        endif
        call mpi_barrier(mpi_comm_comp,irtn)
      enddo
!-----------------------------------------------------------------------
!***  find the maximum horizontal size of each task's subdomain
!***  since task 0 will need that in subroutine dstrb.
!-----------------------------------------------------------------------
!
      num_pts_max=0
!
      if(mype==0)then
        do ipe=1,npes-1
          call mpi_recv(limits,4,mpi_integer,ipe,ipe,mpi_comm_comp      &
     &,                 jstat,irecv)
!
          istart=limits(1)
          iend=limits(2)
          jstart=limits(3)
          jend=limits(4)
!
          num_pts=(iend-istart+1)*(jend-jstart+1)
          if(num_pts>num_pts_max)then
            num_pts_max=num_pts
          endif
        enddo
!
      else
!
        limits(1)=its
        limits(2)=ite
        limits(3)=jts
        limits(4)=jte
!
        call mpi_send(limits,4,mpi_integer,0,mype,mpi_comm_comp,isend)
      endif
!
!-----------------------------------------------------------------------
!***  let each task determine who its eight neighbors are because we
!***  will need to know that for the halo exchanges.  the direction
!***  to each neighbor will be designated by the following integers:
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
!***  if a task has no neighbor in a particular direction because of
!***  the presence of the global domain boundary then that element
!***  of my_neb is set to -1.
!-----------------------------------------------------------------------
!
      allocate(itemp(inpes,jnpes),stat=istat)
      ipe=0
!
      do j=1,jnpes
      do i=1,inpes
        itemp(i,j)=ipe
        if(ipe==mype)then
          myi=i
          myj=j
        endif
        ipe=ipe+1
      enddo
      enddo
!
      my_n=-1
      if(myj+1<=jnpes)my_n=itemp(myi,myj+1)
!
      my_e=-1
      if(myi+1<=inpes)my_e=itemp(myi+1,myj)
!
      my_s=-1
      if(myj-1>=1)my_s=itemp(myi,myj-1)
!
      my_w=-1
      if(myi-1>=1)my_w=itemp(myi-1,myj)
!
      my_ne=-1
      if((myi+1<=inpes).and.(myj+1<=jnpes)) &
         my_ne=itemp(myi+1,myj+1)
!
      my_se=-1
      if((myi+1<=inpes).and.(myj-1>=1)) &
         my_se=itemp(myi+1,myj-1)
!
      my_sw=-1
      if((myi-1>=1).and.(myj-1>=1)) &
         my_sw=itemp(myi-1,myj-1)
!
      my_nw=-1
      if((myi-1>=1).and.(myj+1<=jnpes)) &
         my_nw=itemp(myi-1,myj+1)
!
      my_neb(1)=my_n
      my_neb(2)=my_e
      my_neb(3)=my_s
      my_neb(4)=my_w
      my_neb(5)=my_ne
      my_neb(6)=my_se
      my_neb(7)=my_sw
      my_neb(8)=my_nw
!
      deallocate(itemp)
!
!-----------------------------------------------------------------------
!
!     write(0,*)' exit decomp ite=',ite
      end subroutine decomp
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine dstrb &
      (arrayg,arrayl,lgs,lge,lls,lle,l1)
!-----------------------------------------------------------------------
!     dstrb distributes the elements of real global array arrayg to the
!     real local array arrayl. 
!-----------------------------------------------------------------------
!     input argument list:
!       arrayg - global array
!       lgs    - starting vertical index of global array
!       lge    - ending vertical index of global array
!       lls    - starting vertical index of local array
!       lle    - ending vertical index of local array
!       l1     - vertical level of arrayl being filled in this call
!                (used only when lge=1 and lle>1, i.e. when the global
!                 array is actually just one level of a multi_level
!                 array)
!
!     output argument list:
!       arrayl - local array
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***
!***  argument variables
!***
      integer(kind=kint),intent(in) :: l1,lge,lgs,lle,lls
!
      real(kind=kfpt),dimension(ids:ide,jds:jde,lgs:lge),intent(in) :: &
                                                                  arrayg
      real(kind=kfpt),dimension(ims:ime,jms:jme,lls:lle),intent(out) :: &
                                                                  arrayl
!
!-----------------------------------------------------------------------
!***
!***  local variables
!***
!
      integer(kind=kfpt) :: i,iend,ipe,irecv,irtn,isend,istart,j,jend &
                           ,jstart,knt,l,mype,numvals
      integer,dimension(4) :: limits
      integer,dimension(mpi_status_size) :: jstat
!
      real(kind=kfpt),allocatable,dimension(:) :: arrayx
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  task 0 fills its own local domain then parcels out all the other 
!***  pieces to the other tasks.
!-----------------------------------------------------------------------
!
      mype=mype_share
!
      tasks : if(mype==0)then
!
        if(lge==lgs)then
          do j=jts,jte
          do i=its,ite
            arrayl(i,j,l1)=arrayg(i,j,lgs)
          enddo
          enddo
!
        else
!
          do l=lgs,lge
            do j=jts,jte
            do i=its,ite
              arrayl(i,j,l)=arrayg(i,j,l)
            enddo
            enddo
          enddo
        endif
!
!***  task 0 creates an array to hold all the values from the other
!***  tasks' pieces of the global array and then sends those pieces 
!***  out to the appropriate tasks.
!
        numvals=num_pts_max*(lge-lgs+1)
        allocate(arrayx(numvals),stat=i)
!
        do ipe=1,npes-1
!
          call mpi_recv(limits,4,mpi_integer,ipe,ipe,mpi_comm_comp &
                       ,jstat,irecv)
!
          istart=limits(1)
          iend=limits(2)
          jstart=limits(3)
          jend=limits(4)
          knt=0
!
          do l=lgs,lge
            do j=jstart,jend
            do i=istart,iend
              knt=knt+1
              arrayx(knt)=arrayg(i,j,l)
            enddo
            enddo
          enddo
!
          call mpi_send(arrayx,knt,mpi_real,ipe,ipe,mpi_comm_comp,isend)
!
        enddo
!
        deallocate(arrayx)
!
!-----------------------------------------------------------------------
!***  all other tasks tell task 0 what their horizontal limits are and
!***  receive their piece of the global array from task 0.
!-----------------------------------------------------------------------
!
      else
!
        limits(1)=its
        limits(2)=ite
        limits(3)=jts
        limits(4)=jte
!
        call mpi_send(limits,4,mpi_integer,0,mype,mpi_comm_comp,isend)
!
        knt=(ite-its+1)*(jte-jts+1)*(lge-lgs+1)
        allocate(arrayx(knt),stat=i)
!
        call mpi_recv(arrayx,knt,mpi_real,0,mype,mpi_comm_comp &
                     ,jstat,irecv)
!
        knt=0
        if(lge==lgs)then
          do j=jts,jte
          do i=its,ite
            knt=knt+1
            arrayl(i,j,l1)=arrayx(knt)
          enddo
          enddo
        else
          do l=lgs,lge
            do j=jts,jte
            do i=its,ite
              knt=knt+1
              arrayl(i,j,l)=arrayx(knt)
            enddo
            enddo
          enddo
        endif
!
        deallocate(arrayx)
!
!-----------------------------------------------------------------------
!
      endif tasks
!
!-----------------------------------------------------------------------
      call mpi_barrier(mpi_comm_comp,irtn)
!-----------------------------------------------------------------------
!
      end subroutine dstrb
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine idstrb &
      (iarrayg,iarrayl)
!-----------------------------------------------------------------------
!     idstrb distributes the elements of integer global array iarrayg
!     to the integer local array iarrayl. 
!-----------------------------------------------------------------------
!     input argument list:
!       iarrayg - global array
!
!     output argument list:
!       iarrayl - local array
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***
!***  argument variables
!***
      integer(kind=kint),dimension(ids:ide,jds:jde),intent(in) :: &
                                                                 iarrayg
      integer(kind=kint),dimension(ims:ime,jms:jme),intent(out) :: &
                                                                 iarrayl
!-----------------------------------------------------------------------
!***
!***  local variables
!***
!
      integer(kind=kfpt) :: i,iend,ipe,irecv,irtn,isend,istart,j,jend &
                           ,jstart,knt,l,mype,numvals
      integer,dimension(4) :: limits
      integer,dimension(mpi_status_size) :: jstat
!
      integer(kind=kint),allocatable,dimension(:) :: iarrayx
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!***  initialize the output array.
!
      do j=jms,jme
      do i=ims,ime
        iarrayl(i,j)=0.
      enddo
      enddo
!
!-----------------------------------------------------------------------
!***  task 0 fills its own local domain then parcels out all the other 
!***  pieces to the other tasks.
!-----------------------------------------------------------------------
!
      tasks : if(mype==0)then
!
        do j=jts,jte
        do i=its,ite
          iarrayl(i,j)=iarrayg(i,j)
        enddo
        enddo
!
!***  task 0 creates an array to hold all the values from the other
!***  tasks' pieces of the global array and then sends those pieces 
!***  out to the appropriate tasks.
!
        numvals=num_pts_max
        allocate(iarrayx(numvals),stat=i)
!
        do ipe=1,npes-1
!
          call mpi_recv(limits,4,mpi_integer,ipe,ipe,mpi_comm_comp &
                       ,jstat,irecv)
!
          istart=limits(1)
          iend=limits(2)
          jstart=limits(3)
          jend=limits(4)
          knt=0
!
          do j=jstart,jend
          do i=istart,iend
            knt=knt+1
            iarrayx(knt)=iarrayg(i,j)
          enddo
          enddo
!
          call mpi_send(iarrayx,knt,mpi_integer,ipe,ipe,mpi_comm_comp &
                       ,isend)
!
        enddo
!
        deallocate(iarrayx)
!
!-----------------------------------------------------------------------
!***  all other tasks tell task 0 what their horizontal limits are and
!***  receive their piece of the global array from task 0.
!-----------------------------------------------------------------------
!
      else
!
        limits(1)=its
        limits(2)=ite
        limits(3)=jts
        limits(4)=jte
!
        call mpi_send(limits,4,mpi_integer,0,mype,mpi_comm_comp,isend)
!
        knt=(ite-its+1)*(jte-jts+1)
        allocate(iarrayx(knt),stat=i)
!
        call mpi_recv(iarrayx,knt,mpi_integer,0,mype,mpi_comm_comp &
                     ,jstat,irecv)
!
        knt=0
        do j=jts,jte
        do i=its,ite
          knt=knt+1
          iarrayl(i,j)=iarrayx(knt)
        enddo
        enddo
!
        deallocate(iarrayx)
!
!-----------------------------------------------------------------------
!
      endif tasks
!
!-----------------------------------------------------------------------
      call mpi_barrier(mpi_comm_comp,irtn)
!-----------------------------------------------------------------------
!
      end subroutine idstrb
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine dstrb_soil &
      (arrayg,arrayl,lgs,lge,lls,lle)
!-----------------------------------------------------------------------
!     dstrb distributes the elements of real global array arrayg to the
!     real local array arrayl. 
!-----------------------------------------------------------------------
!     input argument list:
!       arrayg - global soil array
!       lgs    - starting vertical index of global array
!       lge    - ending vertical index of global array
!       lls    - starting vertical index of local array
!       lle    - ending vertical index of local array
!
!     output argument list:
!       arrayl - local soil array
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***
!***  argument variables
!***
      integer(kind=kint),intent(in) :: lge,lgs,lle,lls
!
      real(kind=kfpt),dimension(lgs:lge,ids:ide,jds:jde),intent(in) :: &
                                                                  arrayg
      real(kind=kfpt),dimension(lls:lle,ims:ime,jms:jme),intent(out) :: &
                                                                  arrayl
!-----------------------------------------------------------------------
!***
!***  local variables
!***
!
      integer(kind=kfpt) :: i,iend,ipe,irecv,irtn,isend,istart,j,jend &
                           ,jstart,knt,l,mype,numvals
      integer,dimension(4) :: limits
      integer,dimension(mpi_status_size) :: jstat
!
      real(kind=kfpt),allocatable,dimension(:) :: arrayx
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!***  initialize the output array.
!
      do j=jms,jme
      do i=ims,ime
      do l=lls,lle
        arrayl(l,i,j)=0.
      enddo
      enddo
      enddo
!
!-----------------------------------------------------------------------
!***  task 0 fills its own local domain then parcels out all the other 
!***  pieces to the other tasks.
!-----------------------------------------------------------------------
!
      tasks : if(mype==0)then
!
        do j=jts,jte
        do i=its,ite
          do l=lgs,lge
            arrayl(l,i,j)=arrayg(l,i,j)
          enddo
        enddo
        enddo
!
!***  task 0 creates an array to hold all the values from the other
!***  tasks' pieces of the global array and then sends those pieces 
!***  out to the appropriate tasks.
!
        numvals=num_pts_max*(lge-lgs+1)
        allocate(arrayx(numvals),stat=i)
!
        do ipe=1,npes-1
!
          call mpi_recv(limits,4,mpi_integer,ipe,ipe,mpi_comm_comp &
                       ,jstat,irecv)
!
          istart=limits(1)
          iend=limits(2)
          jstart=limits(3)
          jend=limits(4)
          knt=0
!
          do j=jstart,jend
          do i=istart,iend
            do l=lgs,lge
              knt=knt+1
              arrayx(knt)=arrayg(l,i,j)
            enddo
          enddo
          enddo
!
          call mpi_send(arrayx,knt,mpi_real,ipe,ipe,mpi_comm_comp,isend)
!
        enddo
!
        deallocate(arrayx)
!
!-----------------------------------------------------------------------
!***  all other tasks tell task 0 what their horizontal limits are and
!***  receive their piece of the global array from task 0.
!-----------------------------------------------------------------------
!
      else
!
        limits(1)=its
        limits(2)=ite
        limits(3)=jts
        limits(4)=jte
!
        call mpi_send(limits,4,mpi_integer,0,mype,mpi_comm_comp,isend)
!
        knt=(ite-its+1)*(jte-jts+1)*(lge-lgs+1)
        allocate(arrayx(knt),stat=i)
!
        call mpi_recv(arrayx,knt,mpi_real,0,mype,mpi_comm_comp &
                     ,jstat,irecv)
!
        knt=0
        do j=jts,jte
        do i=its,ite
          do l=lgs,lge
            knt=knt+1
            arrayl(l,i,j)=arrayx(knt)
          enddo
        enddo
        enddo
!
        deallocate(arrayx)
!
!-----------------------------------------------------------------------
!
      endif tasks
!
!-----------------------------------------------------------------------
      call mpi_barrier(mpi_comm_comp,irtn)
!-----------------------------------------------------------------------
!
      end subroutine dstrb_soil
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine gather_layers &
      (field,lm,npes,msize_dummy_fft &
      ,lm_fft,k1_fft,k2_fft &
      ,local_istart,local_iend &
      ,local_jstart,local_jend &
      ,jstart_fft,jend_fft &
      ,ipe_start,ipe_end &
      ,my_domain_has_fft_lats &
      ,array_lyrs)
!-----------------------------------------------------------------------
!***  gather_layers distributes all the elements of field between layers
!***  k1 and k2 inclusive to the appropriate task for subsequent 
!***  application of ffts.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 ipe_end &
,ipe_start &
,jstart_fft &
,jend_fft &
,lm &
,lm_fft &
,msize_dummy_fft &
,npes 
!
integer(kind=kint),dimension(0:npes-1),intent(in) :: &
 k1_fft &
,k2_fft &
,local_iend &
,local_istart &
,local_jend &
,local_jstart 
!
real(kind=kfpt),dimension(ims:ime,jms:jme,lm),intent(in) :: &
 field
!
real(kind=kfpt),dimension(ids:ide,jstart_fft:jend_fft,1:lm_fft),intent(out) :: &
 array_lyrs
!
logical(kind=klog),dimension(0:npes-1),intent(in) :: &
 my_domain_has_fft_lats
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,i_extent &
,iend &
,irecv &
,isend &
,istart &
,j &
,j_extent &
,jend_fft_local &
,jend_fill_mine &
,jstart_fft_local &
,jstart_fill_mine &
,k &
,k1 &
,k2 &
,k_extent &
,mpe &
,mype &
,n &
,nn &
,npe &
,numvals
!
integer(kind=kint),dimension(mpi_status_size) :: &
 jstat  
!
real(kind=kfpt),dimension(msize_dummy_fft) :: &
 dummy
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!***  each task is responsible for multiple full layers (the entire
!***  horizontal expanse of points in a given model layer in a given
!***  hemisphere that are on latitude circles where ffts are applied),
!***  a single layer, or even a partial layer if there are more mpi
!***  compute tasks than there are layers (by the definition above,
!***  there are 2*lm layers).
!***  each task must gather full latitude circles for its designated
!***  layers and rows and send to the other tasks its pieces of the
!***  latitude circles they need for the model layers they will be
!***  handling.
!
!***  k1_fft and k2_fft provide the vertical limits for each task's
!***  group of model layers.
!
!***  we gather into array array_lyrs which will hold full latitude
!***  circles of data in each task's own set of assigned model layers.
!-----------------------------------------------------------------------
!
!     write(0,*)' enter gather_layers'
!       do n=0,npes-1
!         write(0,*)' gather n=',n,' k1_fft=',k1_fft(n),' k2_fft=',k2_fft(n)
!       enddo
      jstart_fill_mine=max(jts,jstart_fft)
      jend_fill_mine=min(jte,jend_fft)
      j_extent=jend_fft-jstart_fft+1
!
!-----------------------------------------------------------------------
!***
!***  loop through the tasks that will be handling layers with ffts
!***  for the appropriate hemisphere.
!***
!-----------------------------------------------------------------------
      receivers: do mpe=ipe_start,ipe_end
!     write(0,*)' gather_layers point 1 mpe=',mpe
!-----------------------------------------------------------------------
!
!***  we need the model layers handled by receiver mpe.
!
        k1=k1_fft(mpe)
        k2=k2_fft(mpe)
        k_extent=k2-k1+1
!     write(0,*)' k1(mpe)=',k1,' k2(mpe)=',k2,' k_extent=',k_extent
!
!-----------------------------------------------------------------------
        if(mype==mpe)then  ! this task is a receiver
!-----------------------------------------------------------------------
!
!***  first, the receiver fills the working array with its own data
!***  if its subdomain contains fft latitude circles.
!
!     write(0,*)' gather_layers point 2 mype=',mpe
          if(my_domain_has_fft_lats(mype))then
!     write(0,*)' gather_layers point 3'
            n=0
            do k=k1,k2
              n=n+1
!     write(0,*)' gather_layers point 3.1 k=',k
!     write(0,*)' jstart_fft=',jstart_fft,' jend_fft=',jend_fft,' lm_fft=',lm_fft,' lm=',lm
!     write(0,*)' jstart_fill_mine=',jstart_fill_mine,' jend_fill_mine=',jend_fill_mine
              do j=jstart_fill_mine,jend_fill_mine
              do i=its,ite
                array_lyrs(i,j,n)=field(i,j,k)
              enddo
              enddo
            enddo
!     write(0,*)' gather_layers point 4'
          endif
!     write(0,*)' gather_layers point 5'
!
!-----------------------------------------------------------------------
          recv_from_npe: do npe=ipe_start,ipe_end
!     write(0,*)' gather_layers point 6 npe=',npe
!-----------------------------------------------------------------------
!
!***  need the number of points being sent by sender npe.
!***  the receiver uses its k extent, the sender's full i extent,
!***  and the sender's j extent over which its subdomain has
!***  fft latitude circles.
!
            if(mype/=npe.and.my_domain_has_fft_lats(npe))then
!     write(0,*)' gather_layers point 7'
              istart=local_istart(npe)
              iend=local_iend(npe)
              jstart_fft_local=max(local_jstart(npe),jstart_fft)
              jend_fft_local=min(local_jend(npe),jend_fft)
!
!     write(0,*)' npe=',npe,' local_istart=',local_istart(npe),' local_iend=',local_iend(npe)
              i_extent=iend-istart+1
              j_extent=jend_fft_local-jstart_fft_local+1
              numvals=j_extent*k_extent*i_extent
!     write(0,*)' msize_dummy_fft=',msize_dummy_fft
!     write(0,*)' i_extent=',i_extent,' j_extent=',j_extent,' k_extent=',k_extent
!     write(0,*)' gather mpe=',mpe,' to recv numvals=',numvals,' from npe=',npe
!
              call mpi_recv(dummy,numvals,mpi_real,npe,npe &
                           ,mpi_comm_comp,jstat,irecv)
!     write(0,*)' gather_layers point 8'
!
              n=0
              nn=0
              do k=k1,k2
                n=n+1
                do j=jstart_fft_local,jend_fft_local
                do i=istart,iend
                  nn=nn+1
                  array_lyrs(i,j,n)=dummy(nn)
                enddo
                enddo
              enddo
!     write(0,*)' gather_layers point 9'
!
            endif
!     write(0,*)' gather_layers point 10'
!
          enddo recv_from_npe
!     write(0,*)' gather_layers point 11'
!
!-----------------------------------------------------------------------
        else  ! these tasks are senders
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
          senders: do npe=ipe_start,ipe_end
!     write(0,*)' gather_layers point 12 npe=',npe
!-----------------------------------------------------------------------
            if(mype==npe.and.my_domain_has_fft_lats(mype))then
!     write(0,*)' gather_layers point 13'
!
!***  need the number of points sent by sender npe.
!***  the sender uses the receiver's k extent, its own i extent,
!***  and the j extent for which its subdomain has fft latitudes.
!
              jstart_fft_local=max(jts,jstart_fft)
              jend_fft_local=min(jte,jend_fft)
!
!             i_extent=ite-its+1
!     write(0,*)' i_extent=',i_extent,' its=',its,' ite=',ite
!             j_extent=jend_fft_local-jstart_fft_local+1
!             numvals=j_extent*k_extent*i_extent
!     write(0,*)' gather_layers point 14'
!
              n=0
              nn=0
              do k=k1,k2
                n=n+1
                do j=jstart_fft_local,jend_fft_local
                do i=its,ite
                  nn=nn+1
                  dummy(nn)=field(i,j,k)
                enddo
                enddo
              enddo
!     write(0,*)' msize_dummy_fft=',msize_dummy_fft
!     write(0,*)' i_extent=',i_extent,' j_extent=',j_extent,' k_extent=',k_extent
!     write(0,*)' numvals for send=',numvals,' nn=',nn
!     write(0,*)' gather mype=',mype,' to send numvals=',nn,' to mpe=',mpe
!     write(0,*)' gather_layers point 15'
!
!!!!!         call mpi_send(dummy,numvals,mpi_real,mpe,mype &
              call mpi_send(dummy,nn,mpi_real,mpe,mype &
                           ,mpi_comm_comp,isend)
!     write(0,*)' gather_layers point 16'
!
            endif
!     write(0,*)' gather_layers point 17'
!
          enddo senders
!     write(0,*)' gather_layers point 18'
!
!-----------------------------------------------------------------------
        endif
!     write(0,*)' gather_layers point 19'
!-----------------------------------------------------------------------
!
      enddo receivers
!     write(0,*)' exit gather_layers'
!
!-----------------------------------------------------------------------
!
      end subroutine gather_layers
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine scatter_layers &
      (array_lyrs,lm,npes,msize_dummy_fft &
      ,lm_fft,k1_fft,k2_fft &
      ,local_istart,local_iend &
      ,local_jstart,local_jend &
      ,jstart_fft,jend_fft &
      ,ipe_start,ipe_end &
      ,my_domain_has_fft_lats &
      ,field)
!-----------------------------------------------------------------------
!***  scatter_layers distributes the elements of array_lyrs between 
!***  layers k1 and k2 inclusive to the appropriate tasks that actually
!***  own the fft latitude rows.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 ipe_end &
,ipe_start &
,jstart_fft &
,jend_fft &
,lm &
,lm_fft &
,msize_dummy_fft &
,npes 
!
integer(kind=kint),dimension(0:npes-1),intent(in) :: &
 k1_fft &
,k2_fft &
,local_iend &
,local_istart &
,local_jend &
,local_jstart
!
real(kind=kfpt),dimension(ids:ide,jstart_fft:jend_fft,1:lm_fft),intent(in) :: &
 array_lyrs
!
real(kind=kfpt),dimension(ims:ime,jms:jme,lm),intent(out) :: &
 field
!
logical(kind=klog),dimension(0:npes-1),intent(in) :: &
 my_domain_has_fft_lats
!
!-----------------------------------------------------------------------
!***  local variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,i_extent &
,iend &
,irecv &
,isend &
,istart &
,j &
,j_extent &
,jend_fft_local &
,jend_fill_mine &
,jstart_fft_local &
,jstart_fill_mine &
,k &
,k_extent &
,k1 &
,k2 &
,mpe &
,mype &
,n &
,nn &
,npe &
,numvals
!
integer(kind=kint),dimension(mpi_status_size) :: &
 jstat 
!
real(kind=kfpt),dimension(msize_dummy_fft) :: &
 dummy
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!***  each task holds full latitude circles of data within its own
!***  subset of model layers that it was assigned.
!
!***  k1_fft and k2_fft provide the vertical limits for each task's
!***  group of model layers.
!
!-----------------------------------------------------------------------
!
!     write(0,*)' enter scatter msize_dummy_fft=',msize_dummy_fft
      jstart_fill_mine=max(jts,jstart_fft)
      jend_fill_mine=min(jte,jend_fft)
      j_extent=jend_fft-jstart_fft+1
!
!-----------------------------------------------------------------------
!***  loop through the tasks in the hemisphere which in turn will
!***  send to those tasks whose subdomains contain the appropriate
!***  fft latitude rows.
!-----------------------------------------------------------------------
      receivers: do mpe=ipe_start,ipe_end
!-----------------------------------------------------------------------
!     write(0,*)' i am mype=',mype,' mpe=',mpe,' ipe_start=',ipe_start,' ipe_end=',ipe_end
!
!-----------------------------------------------------------------------
        if(mype==mpe)then  ! this task might receive
!-----------------------------------------------------------------------
!
!     write(0,*)' mpe=',mpe,' mype=',mype,' my_domain_has_fft_lats(mype)=',my_domain_has_fft_lats(mype)
          if(my_domain_has_fft_lats(mype))then  ! receive if subdomain has
                                                ! fft latitudes
!     write(0,*)' scatter point 3 my_domain_has_fft_lats=',my_domain_has_fft_lats(mype)
!
!-----------------------------------------------------------------------
!
!***  first, the receiver fills the prognostic array with its own data.
!
            k1=k1_fft(mpe)
            k2=k2_fft(mpe)
!     write(0,*)' scatter point 4 k1_fft=',k1_fft(mpe),' k2_fft=',k2_fft(mpe)
!     write(0,*)' jstart_fill_mine=',jstart_fill_mine,' jend_fill_mine=',jend_fill_mine
            n=0
!
            do k=k1,k2
              n=n+1
              do j=jstart_fill_mine,jend_fill_mine
!     write(0,*)' scatter point 5 k=',k,' j=',j
              do i=its,ite
                field(i,j,k)=array_lyrs(i,j,n)
              enddo
              enddo
            enddo
!
!-----------------------------------------------------------------------
            recv_from_npe: do npe=ipe_start,ipe_end
!-----------------------------------------------------------------------
!
              if(mype/=npe)then
!
!***  need the number of points being sent by each of the other tasks
!***  for their subset of model layers.
!***  the receiver uses the sender's k extent and its own i extent
!***  and the j extent for which its subdomain has fft latitudes.
!
                k1=k1_fft(npe)
                k2=k2_fft(npe)
                k_extent=k2-k1+1
                i_extent=ite-its+1
                jstart_fft_local=max(jts,jstart_fft)
                jend_fft_local=min(jte,jend_fft)
                j_extent=jend_fft_local-jstart_fft_local+1
!
                numvals=j_extent*k_extent*i_extent
!     write(0,*)' npe=',npe,' jstart_fft_local=',jstart_fft_local,' jend_fft_local=',jend_fft_local
!     write(0,*)' jts=',jts,' jstart_fft=',jstart_fft,' jte=',jte,' jend_fft=',jend_fft
!     write(0,*)' i_extent=',i_extent,' j_extent=',j_extent,' k_extent=',k_extent
!     write(0,*)' its=',its,' ite=',ite
!     write(0,*)' mpe=',mpe,' recving numvals=',numvals,' from npe=',npe
!
                call mpi_recv(dummy,numvals,mpi_real,npe,npe &
                             ,mpi_comm_comp,jstat,irecv)
!     write(0,*)' mpe=',mpe,' recvd from npe=',npe
!     write(0,*)' jstart_fft=',jstart_fft,' jend_fft=',jend_fft,' k1=',k1,' k2=',k2
!
                n=0
                nn=0
                do k=k1,k2
                  n=n+1
                  do j=jstart_fft_local,jend_fft_local
                  do i=its,ite
                    nn=nn+1
                    field(i,j,k)=dummy(nn)
                  enddo
                  enddo
                enddo
!     write(0,*)' mpe=',mpe,' filled the piece of field from npe=',npe
!
              endif
!
            enddo recv_from_npe
!
          endif
!
!-----------------------------------------------------------------------
        else  ! this task is a sender
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
          senders: do npe=ipe_start,ipe_end
!-----------------------------------------------------------------------
            if(mype==npe.and.my_domain_has_fft_lats(mpe))then
!
!***  need the number of points sent by sender npe.
!***  the sender uses its own k extent and the receiver's i extent
!***  and the receiver's j extent that spans fft laltitude rows.
!
              k1=k1_fft(npe)
              k2=k2_fft(npe)
              k_extent=k2-k1+1
!
              istart=local_istart(mpe)
              iend=local_iend(mpe)
              i_extent=iend-istart+1
              jstart_fft_local=max(local_jstart(mpe),jstart_fft)
              jend_fft_local=min(local_jend(mpe),jend_fft)
              j_extent=jend_fft_local-jstart_fft_local+1
!             numvals=j_extent*k_extent*i_extent
!     write(0,*)' npe=',npe,' mpe=',mpe,' i_extent=',i_extent,' j_extent=',j_extent,' k_extent=',k_extent
!     write(0,*)' jstart_fft_local=',jstart_fft_local,' jend_fft_local=',jend_fft_local
!     write(0,*)' local_jstart(mpe)=',local_jstart(mpe),' jstart_fft=',jstart_fft
!     write(0,*)' local_jend(mpe)=',local_jend(mpe),' jend_fft=',jend_fft
!
              n=0
              nn=0
              do k=k1,k2
                n=n+1
                do j=jstart_fft_local,jend_fft_local
                do i=istart,iend
                  nn=nn+1
                  dummy(nn)=array_lyrs(i,j,n)
                enddo
                enddo
              enddo
!     write(0,*)' npe=',npe,' sending numvals=',nn,' to mpe =',mpe
!
!             call mpi_send(dummy,numvals,mpi_real,mpe,mype &
              call mpi_send(dummy,nn,mpi_real,mpe,mype &
                           ,mpi_comm_comp,isend)
!     write(0,*)' npe=',npe,' sent to mpe =',mpe
!
            endif
!
          enddo senders
!
!-----------------------------------------------------------------------
        endif
!-----------------------------------------------------------------------
!
      enddo receivers
!
!-----------------------------------------------------------------------
!
      end subroutine scatter_layers
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                       end module module_dm_parallel_gfs
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
