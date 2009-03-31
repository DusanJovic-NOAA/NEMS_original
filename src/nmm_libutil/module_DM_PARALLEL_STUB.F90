!-----------------------------------------------------------------------
                        module module_dm_parallel
!-----------------------------------------------------------------------
!
!***  This module contains all codes directly related to distributed
!***  memory issues except for halo exchange although note that the
!***  halo widths must be set here.
!
!-----------------------------------------------------------------------
!
!
!      use module_include
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!---domain decomposition info-------------------------------------------
!-----------------------------------------------------------------------
!
integer(kind=kint),parameter :: &
 ihalo=3 &                   ! halo width in I direction
,jhalo=3 &                   ! halo width in J direction
!
!***  Hardwire to 100 the maximum number of server groups allowed.
!***  This is far greater than should ever be needed.
!
,max_groups=100 &            ! max number of quilt server groups
!
!***  For now, set the number of threads here to 1.
!***  Clearly, this must be reconciled before actually using threading.
!
,num_tiles=1                 ! number of threads
!
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 ide &                       ! ending data index, x direction
,ids &                       ! starting data index, x direction
,ims &                       ! the starting memory I for each task
,ime &                       ! the ending memory I for each task
,its &                       ! the starting integration I for each task
,ite &                       ! the ending integration I for each task
,jde &                       ! ending data index, y direction
,jds &                       ! starting data index, y direction
,jms &                       ! the starting memory J for each task
,jme &                       ! the ending memory J for each task
,jts &                       ! the starting integration J for each task
,jte &                       ! the ending integration J for each task
,lm  &                       ! the number of atmospheric model layers
,mpi_comm_comp &             ! local mpi communicator
,mpi_comm_inter &            ! intercommunicator for the quilt/write tasks
,mype_share &                ! my task ID to be seen by any USEs
,npes &                      ! total number of forecast tasks after SETUP_SERVERS
!
,num_pts_max                 ! max points in any task's subdomain
!
integer(kind=kint) :: &
 its_b1 &                    ! its AND 1 point from global boundary
,its_b2 &                    ! its AND 2 points from global boundary
,its_h1 &                    ! its AND 1 point into halo
,its_h2 &                    ! its AND 2 points into halo
,its_b1_h1 &                 ! its AND _b1 AND _h1
,its_b1_h2 &                 ! its AND _b1 AND _h2
,ite_b1 &                    ! ite AND 1 point from global boundary
,ite_b2 &                    ! ite AND 2 points from global boundary
,ite_h1 &                    ! ite AND 1 point into halo
,ite_h2 &                    ! ite AND 2 points into halo
,ite_b1_h1 &                 ! ite AND _b1 AND _h1
,ite_b1_h2 &                 ! ite AND _b1 AND _h2
,jts_b1 &                    ! jts AND 1 point from global boundary
,jts_b2 &                    ! jts AND 2 points from global boundary
,jts_h1 &                    ! jts AND 1 point into halo
,jts_h2 &                    ! jts AND 2 points into halo
,jts_b1_h1 &                 ! jts AND _b1 AND _h1
,jts_b1_h2 &                 ! jts AND _b1 AND _h2
,jte_b1 &                    ! jte AND 1 point from global boundary
,jte_b2 &                    ! jte AND 2 points from global boundary
,jte_h1 &                    ! jte AND 1 point into halo
,jte_h2 &                    ! jte AND 2 points into halo
,jte_b1_h1 &                 ! jte AND _b1 AND _h1
,jte_b1_h2                   ! jte AND _b1 AND _h2
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
,num_serv_per_grp            ! number of tasks in each group
!
integer(kind=kint),allocatable,dimension(:) :: &
 local_iend &
,local_istart &
,local_jend &
,local_jstart
!
integer :: mpi_intra
!-----------------------------------------------------------------------
!
       contains
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      subroutine setup_servers(mype,inpes,jnpes,npes                    &
                              ,ngroups_write,write_tasks_per_group      &
                              ,mpi_intra)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!***  SETUP_SERVERS splits the communicator between integration 
!***  and output tasks.
!
!-----------------------------------------------------------------------
!
!   input argument list:
!    mype          - My task ID
!    inpes         - Number of mpi tasks in the X direction
!    jnpes         - Number of mpi tasks in the Y direction
!    npes          - Total number of mpi tasks provided to the job.  As input 
!                    to SETUP_SERVERS it includes the forecast tasks plus all
!                    write tasks in all groups of write tasks.  npes must at least
!                    equal the product of inpes*jnpes otherwise the integration
!                    cannot proceed. The difference between the product npes_fcst
!                    and npes is the number of mpi tasks that are available
!                    for i/o serving. This can be zero, in which case output will
!                    write a direct access file that can be separately quilted. 
!                    In order to skip the separate quilting step, make sure that
!                    the number of mpi tasks that the code is initiated with is at
!                    least one greater than npes_fcst.
!                    Later in the routine, npes is reset to the number of fcst tasks.
!    ngroups_write - Number of groups of write tasks.
!    mpi_intra     - The global communicator.
!    write_tasks_per_group - # of write tasks per write group
!
!-----------------------------------------------------------------------
!***  Argument variables.
!-----------------------------------------------------------------------
!
      implicit none
      integer(kind=kint),intent(in) :: &
 mype &                     ! each task's ID
,inpes &                    ! number of compute tasks in X direction
,jnpes &                    ! number of compute tasks in Y direction
,ngroups_write &            ! number of groups of write tasks
,write_tasks_per_group &    ! number of groups of write tasks per group
,mpi_intra                  ! global communicator
!
      integer(kind=kint),intent(inout) :: &
 npes                       ! total number of tasks provided
                            ! then converted to the number of fcst tasks
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
!!!   (n1,n2,nprocs,irank,ista,iend)
      (jnpes,ntasks_per_group,n_write_task,jrow_first,jrow_last)
!
!-----------------------------------------------------------------------
!***  Each write task will receive history data from a subset of the
!***  forecast tasks.  Entire rows (not partial rows) of forecast
!***  tasks will send to a write task.  Determine the range of rows
!***  of forecast tasks that will be sending data to each write
!***  task.  Row 1 is the row of forecast tasks along the domain's
!***  southern boundary.
!-----------------------------------------------------------------------
!
!   input argument list:
!     jnpes            - # of forecast tasks in the J direction
!     ntasks_per_group - # of write tasks per write group
!     n_write_task     - the "index" of the write task being
!                        considered given that every write group
!                        contains tasks 1-->ntasks_per_group
!
!   output argument list:
!     jrow_first - the first row of forecast tasks to send history
!                  data to this write task
!     jrow_last  - the last row of forecast tasks to send history
!                  data to this write task
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) ::                                  &
        jnpes                                                           &
       ,n_write_task                                                    &
       ,ntasks_per_group
!
      integer(kind=kint),intent(out) ::                                 &
        jrow_first                                                      &
       ,jrow_last
!
!
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
      (mype,inpes,jnpes,npes_fcst,im,jm,lmx,global,ijcount)
!
!-----------------------------------------------------------------------
!
!***  DECOMP specifies the domain decomposition.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  Argument variables.
!-----------------------------------------------------------------------
!
      integer(kind=kint),intent(in) :: im &
                                      ,jm &
                                      ,lmx &
                                      ,inpes &
                                      ,jnpes &
                                      ,mype &
                                      ,npes_fcst
!
      logical,intent(in) :: global
!
      integer(kind=kint),dimension(2),intent(in) :: ijcount
!
!
!-----------------------------------------------------------------------
!
      end subroutine decomp
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine dstrb &
      (arrayg,arrayl,lgs,lge,lls,lle,l1)
!-----------------------------------------------------------------------
!     DSTRB distributes the elements of real global array arrayg to the
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
!
      end subroutine dstrb
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine idstrb &
      (iarrayg,iarrayl)
!-----------------------------------------------------------------------
!     IDSTRB distributes the elements of integer global array iarrayg
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
      integer(kind=kint),dimension(ids:ide,jds:jde),intent(in) :: &
                                                                 iarrayg
      integer(kind=kint),dimension(ims:ime,jms:jme),intent(out) :: &
                                                                 iarrayl
      
!
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
!     DSTRB distributes the elements of real global array arrayg to the
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
!***  GATHER_LAYERS distributes all the elements of field between layers
!***  k1 and k2 inclusive to the appropriate task for subsequent 
!***  application of FFTs.
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
!***  SCATTER_LAYERS distributes the elements of array_lyrs between 
!***  layers k1 and k2 inclusive to the appropriate tasks that actually
!***  own the FFT latitude rows.
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
!
!-----------------------------------------------------------------------
!
      end subroutine scatter_layers
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                       end module module_dm_parallel
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
