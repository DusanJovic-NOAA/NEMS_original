!-----------------------------------------------------------------------
!
      MODULE MODULE_WRITE_ROUTINES_GFS
!
!-----------------------------------------------------------------------
!***  THIS MODULE CONTAINS ROUTINES NEEDED BY THE RUN STEP OF THE
!***  WRITE GRIDDED COMPONENT IN WHICH HISTORY OUTPUT DATA FROM
!***  THE FORECAST TASKS ARE ASSEMBLED AND WRITTEN TO HISTORY FILES
!***  BY THE WRITE TASKS.
!-----------------------------------------------------------------------
!***
!***  HISTORY   
!***
!       07 May 2009:  J. Wang - adopt write _routine from NMMB_io
!                                modified for GFS
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
      USE MODULE_WRITE_INTERNAL_STATE_GFS
!
      USE MODULE_GFS_MPI_DEF, ONLY :  MPI_COMM_INTER_ARRAY                 &
                                     ,N_GROUP
!
      USE MODULE_INCLUDE_GFS
!
      USE MODULE_ERR_MSG,ONLY: ERR_MSG,MESSAGE_CHECK
!
      USE NEMSIO_MODULE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------

      PRIVATE
!
      PUBLIC :: FIRST_PASS_GFS                                          &
               ,WRITE_ASYNC_GFS                                         &
               ,WRITE_INIT_GFS                                          &
               ,WRITE_NEMSIO_OPEN                                       
!
!-----------------------------------------------------------------------
!
      INTEGER :: MYPE                                                     !<-- My MPI task ID
!
!-----------------------------------------------------------------------
!
      CONTAINS
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE FIRST_PASS_GFS(IMP_STATE_WRITE                         &
                               ,WRT_INT_STATE                           &
                               ,NTASKS                                  &
                               ,MYPE                                    &
                               ,NBDL                                    &
                               ,NCURRENT_GROUP                          &
                                )
! 
!-----------------------------------------------------------------------
!***  EACH TIME A NEW GROUP OF WRITE TASKS IS INVOKED FOR THE FIRST
!***  TIME THIS ROUTINE WILL PERFORM CERTAIN COMPUTATIONS AND TASKS
!***  THAT ONLY NEED TO BE DONE ONCE FOR EACH WRITE GROUP.
!***  THE ROUTINE WILL UNLOAD SOME VALUES FROM THE WRITE COMPONENT'S
!***  IMPORT STATE.  BECAUSE THESE QUANTITIES DO NOT CHANGE WITH
!***  FORECAST TIME, THE ROUTINE IS NOT NEEDED FOR SUBSEQUENT USE
!***  OF EACH WRITE GROUP.  THE DATA BEING UNLOADED CONSISTS OF
!***  EVERYTHING EXCEPT THE 2D/3D GRIDDED FORECAST ARRAYS.
!***  ALSO BASIC INFORMATION IS PROVIDED TO FORECAST AND WRITE TASKS
!***  THAT WILL BE NECESSARY IN THE QUILTING OF LOCAL 2-D GRIDDED
!***  HISTORY DATA INTO FULL DOMAIN 2-D FIELDS.
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: NTASKS
      INTEGER,INTENT(IN) :: MYPE
      INTEGER,INTENT(IN) :: NCURRENT_GROUP
      INTEGER,INTENT(IN) :: NBDL
!
      TYPE(ESMF_State)          ,INTENT(INOUT) :: IMP_STATE_WRITE
      TYPE(WRITE_INTERNAL_STATE_GFS),INTENT(INOUT) :: WRT_INT_STATE
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_FieldBundle)       :: FILE_BUNDLE
      INTEGER                      :: I,IERR,IM,J,JM,L,K                &
                                     ,N,NN,NUM_ATTRIB,NWTPG             &
                                     ,RC,RC_WRT
!
      INTEGER,DIMENSION(:),POINTER :: INPES,JNPES,ITMP
      REAL,DIMENSION(:),POINTER    :: RTMP
      TYPE(ESMF_Logical),DIMENSION(:),POINTER :: LTMP
!
      INTEGER                      :: JROW_FIRST,JROW_LAST,JROWS
!
      INTEGER                      :: LAST_FCST_TASK                    &
                                     ,LEAD_WRITE_TASK                   &
                                     ,LAST_WRITE_TASK
!
      INTEGER,SAVE                 :: NCHAR_I1D                         &
                                     ,NCHAR_R1D                         &
                                     ,NCHAR_LOG
!
      INTEGER                      :: JEND_WRITE,JSTA_WRITE
!
      INTEGER                      :: KOUNT_I1D                         &
                                     ,KOUNT_I2D                         &
                                     ,KOUNT_R1D                         &
                                     ,KOUNT_R2D                         &
                                     ,KOUNT_LOG
!
      INTEGER                      :: LENGTH                            &
                                     ,LENGTH_SUM_I1D                    &
                                     ,LENGTH_SUM_R1D                    &
                                     ,LENGTH_SUM_LOG
!
      INTEGER                      :: NPOSN_START,NPOSN_END
!
      INTEGER                      :: NUM_FIELD_NAMES                   &
                                     ,NUM_PES_FCST                      &
                                     ,MYPE_LOCAL   
!
      INTEGER,DIMENSION(MPI_STATUS_SIZE) :: JSTAT
!
      INTEGER,DIMENSION(:),POINTER :: LOCAL_ISTART                      &
                                     ,LOCAL_IEND                        &
                                     ,LOCAL_JSTART                      &
                                     ,LOCAL_JEND
!jw:gfs
      INTEGER                      :: nbelt,nremain,ISTART,LAT
      INTEGER,DIMENSION(:),allocatable :: fcst_lat_to_write_task
      CHARACTER(NAME_MAXSTR),DIMENSION(:),allocatable :: field_name
      CHARACTER(NAME_MAXSTR*MAX_DATA_R2D)             :: NAMETMP
!
      INTEGER,DIMENSION(:),POINTER :: NCHAR_I2D                         &
                                     ,NCHAR_R2D
!
      INTEGER,DIMENSION(:),POINTER :: WORK_ARRAY_I1D
!
      REAL(4),DIMENSION(:),POINTER :: WORK_ARRAY_R1D
      REAL(8),DIMENSION(:),POINTER :: WORK_ARRAY_R1D8
!
      CHARACTER(ESMF_MAXSTR)       :: ATTRIB_NAME
!
      TYPE(ESMF_TypeKind)          :: DATATYPE
!
      TYPE(ESMF_Field)             :: FIELD_WORK1
!
      TYPE(ESMF_Logical)              :: WORK_LOGICAL
      TYPE(ESMF_Logical),DIMENSION(1) :: NO_FIELDS
!
      TYPE(ESMF_VM)                :: VM
!
!-----------------------------------------------------------------------
!
      REAL(KIND=kind_evod)              :: btim,btim0
!
!-----------------------------------------------------------------------
!***** part 1: necessary info for fcst tasks to send data to write tasks
!***********************************************************************
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  FIRST WE NEED THE NUMBER OF WRITE TASKS IN EACH GROUP.
!-----------------------------------------------------------------------
!
      NWTPG=wrt_int_state%WRITE_TASKS_PER_GROUP
      NUM_PES_FCST=wrt_int_state%INPES*wrt_int_state%JNPES               !<-- Number of fcst tasks
!
      LAST_FCST_TASK =NTASKS-NWTPG-1
      LEAD_WRITE_TASK=LAST_FCST_TASK+1
      LAST_WRITE_TASK=NTASKS-1

      write(0,*)'in first pass, LAST_FCST_TASK=',LAST_FCST_TASK,  &
        'LEAD_WRITE_TASK=',LEAD_WRITE_TASK,'LAST_WRITE_TASK=',  &
        LAST_WRITE_TASK,'NWTPG=',NWTPG,'NUM_PES_FCST=',NUM_PES_FCST
!
!-----------------------------------------------------------------------
!***  THE INTEGER QUANTITIES 'INPES' AND 'JNPES' MUST BE POINTERS
!***  OF LENGTH 1 SINCE THEY WILL BE PASSED THROUGH ESMF Send/Recv.
!***  LIKEWISE WITH THE TOTAL LENGTH OF ALL NAMES OF 2D DATA 
!***  AND THE HALO DEPTHS.
!-----------------------------------------------------------------------
!
      ALLOCATE(INPES(1))
      ALLOCATE(JNPES(1))
      ALLOCATE(ITMP(10000))
      ALLOCATE(RTMP(50000))
      ALLOCATE(LTMP(5000))
!jw      ALLOCATE(IHALO(1))
!jw      ALLOCATE(JHALO(1))
      ALLOCATE(NCHAR_I2D(1))
      ALLOCATE(NCHAR_R2D(1))
!
!-----------------------------------------------------------------------
!***  EXTRACT THE FULL DOMAIN LIMITS FROM THE COMPONENT'S IMPORT 
!***  STATE.  THESE ARE NEEDED FOR ALLOCATING THE WORKING ARRAYS 
!***  THAT WILL MOVE THE 2-D AND 3-D FIELDS FROM THE IMPORT TO THE
!***  EXPORT STATE.
!-----------------------------------------------------------------------
!
        write(0,*)'in first_pass 0, INPES=',wrt_int_state%INPES,  &
         'jnpes=',wrt_int_state%JNPES,'im=',wrt_int_state%im(1), &
         'jm=',wrt_int_state%jm(1),'lm=',wrt_int_state%lm(1)
!
!-----------------------------------------------------------------------
!
      domain_limits: IF(MYPE<=LAST_FCST_TASK)THEN                          !<-- This selects only forecast tasks to do extractions
                                                                           !    since only they know what is in the import state
!-----------------------------------------------------------------------
!*** get bundle
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract History Bundle from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_StateGet(state      =IMP_STATE_WRITE                  &  !<-- The write component's import state
                          ,itemName   =wrt_int_state%filename_base(NBDL)&  !<-- The name of the history data Bundle
                          ,fieldbundle=FILE_BUNDLE                      &  !<-- The history data Bundle inside the import state
                          ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
       write(0,*)'in first_pass, filebundle name=',wrt_int_state%filename_base(NBDL),'RC=',RC
!
!--------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Global Parameters from History Bundle" 
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(bundle    =FILE_BUNDLE                &  !<-- The Bundle of history data
                              ,name      ='lonf'                        &  !<-- Name of the Attribute to extract
                              ,count     =1                             &  !<-- Length of Attribute
                              ,valueList =wrt_int_state%IM              &  !<-- Extract this Attribute from History Bundle
                              ,rc        =RC)
!
!jw        CALL ESMF_AttributeGet(bundle    =FILE_BUNDLE                &  !<-- The Bundle of history data
!jw                              ,name      ='latg'                        &  !<-- Name of the Attribute to extract
!jw                              ,count     =1                             &  !<-- Length of Attribute
!jw                              ,valueList =wrt_int_state%JM              &  !<-- Extract this Attribute from History Bundle
!jw                              ,rc        =RC)
!
        CALL ESMF_AttributeGet(bundle    =FILE_BUNDLE                &  !<-- The Bundle of history data
                              ,name      ='levs'                        &  !<-- Name of the Attribute to extract
                              ,count     =1                             &  !<-- Length of Attribute
                              ,valueList =wrt_int_state%LM              &  !<-- Extract this Attribute from History Bundle
                              ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  NOW EXTRACT LOCAL SUBDOMAIN LIMITS.
!***  THESE WILL BE USED TO ALLOCATE THE WORKING ARRAY TO HOLD FIELDS
!***  ON EACH SUBDOMAIN PRIOR TO QUILTING THEM TOGETHER.
!***  WE FIRST NEED THE NUMBER OF FORECAST TASKS SINCE THAT
!***  DETERMINES THE SIZE OF THE ARRAYS HOLDING THE LOCAL
!***  SUBDOMAIN LIMITS.
!
!***  ALSO EXTRACT THE HALO DEPTHS SINCE THEY ARE NEEDED FOR
!***  EXCLUDING HALO POINTS FROM THE FINAL HISTORY DATA.
!
!***  THESE VALUES ARE NOT TO BE WRITTEN TO THE HISTORY FILES SO
!***  THEY WERE NOT INSERTED INTO THE HISTORY DATA Bundle INSIDE
!***  THE WRITE COMPONENT'S IMPORT STATE.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Local Quilting Info from Write Import State"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        write(0,*)'in first_pass, INPES=',wrt_int_state%INPES,  &
         'jnpes=',wrt_int_state%JNPES,'im=',wrt_int_state%im(1), &
         'jm=',wrt_int_state%jm(1),'lm=',wrt_int_state%lm(1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!

        write(0,*)'in first_pass, ipt_lats_node_a=',wrt_int_state%ipt_lats_node_a,  &
         'lats_node_a=',wrt_int_state%lats_node_a,'global_lats_a=', &
           wrt_int_state%global_lats_a(1:wrt_int_state%jm(1)),'jm=',wrt_int_state%jm(1)
!
!-----------------------------------------------------------------------
!#####jw #############     for 1 tasks
!-----------------------------------------------------------------------
!
        if1tasks:  if(NTASKS>1 ) then
!
!-----------------------------------------------------------------------
!*** set up an array of the write_tasks for each fcst task to send data so
!*** the data on write PEs will be continuously distributed in the belted 
!*** subdomain.
!-----------------------------------------------------------------------
        if(.not.allocated(wrt_int_state%nlat_write_task)) then
          allocate(wrt_int_state%nlat_write_task(NWTPG))
          allocate(wrt_int_state%nstart_write_task(NWTPG))
          allocate(wrt_int_state%fcst_lat_for_write_task(wrt_int_state%lats_node_a))
          allocate(wrt_int_state%nwrttask_on_fcst(wrt_int_state%lats_node_a))
        endif

        allocate(fcst_lat_to_write_task(wrt_int_state%lats_node_a))
        NBELT=wrt_int_state%jm(1)/NWTPG
        Nremain=mod(wrt_int_state%jm(1),NWTPG)
!
        do i=1,wrt_int_state%lats_node_a
          lat=wrt_int_state%global_lats_a(wrt_int_state%ipt_lats_node_a-1+i)
          if(lat<=nremain*(nbelt+1) ) then
            fcst_lat_to_write_task(i)=(lat-1)/(nbelt+1)+1
          else
            fcst_lat_to_write_task(i)=nremain+ (lat-nremain*(nbelt+1)-1)/Nbelt+1
          endif
        enddo
        write(0,*)'nbelt=',nbelt,'nremain=',nremain,'lats_nodes=', &
         wrt_int_state%lats_node_a,'ipt_lats_node_a=',wrt_int_state%ipt_lats_node_a, &
         'global_lats_a=', &
         wrt_int_state%global_lats_a(wrt_int_state%ipt_lats_node_a: &
         wrt_int_state%ipt_lats_node_a-1+wrt_int_state%lats_node_a),'fcst_lat_to_write_task=',  &
         fcst_lat_to_write_task
!
!jw*** on each forecast task, specify the starting point and the number of lats in global_lat_a to
!jw*** to each write task, which will decide the position of lat on write tasks
!jw*** 
        n=0
        wrt_int_state%nstart_write_task(1)=1
        wrt_int_state%fcst_lat_for_write_task(:)=-1.0E6
        do i=1,NWTPG
          wrt_int_state%nlat_write_task(i)=0
          do j=1,wrt_int_state%lats_node_a
            if( fcst_lat_to_write_task(j)==i) then
              n=n+1
              wrt_int_state%nlat_write_task(i)=wrt_int_state%nlat_write_task(i)+1
              lat=j-1+wrt_int_state%ipt_lats_node_a
              wrt_int_state%fcst_lat_for_write_task(n)=wrt_int_state%global_lats_a(lat)
              wrt_int_state%nwrttask_on_fcst(n)=j
            endif
          enddo 
          if(i>1) then
            wrt_int_state%nstart_write_task(i)=wrt_int_state%nstart_write_task(i-1)+ &
                        wrt_int_state%nlat_write_task(i-1)
          endif
       enddo
       write(0,*)'nstart_write_task=',wrt_int_state%nstart_write_task(1:NWTPG),'nlat_write_task=', &
          wrt_int_state%nlat_write_task(1:nwtpg),'fcst_lat_for_write_task=', &
          wrt_int_state%fcst_lat_for_write_task(1:wrt_int_state%lats_node_a), &
          'nwrttask_on_fcst=',wrt_int_state%nwrttask_on_fcst(1:wrt_int_state%lats_node_a)
!
        deallocate(fcst_lat_to_write_task)
!
!-----------------------------------------------------------------------
!#####jw #############     for 1 tasks
!-----------------------------------------------------------------------
!
        endif if1tasks
!
!-----------------------------------------------------------------------
!
      ENDIF domain_limits
!
!
!-----------------------------------------------------------------------
!#####jw #############     for 1 tasks
!-----------------------------------------------------------------------
!
        if1task1:  if(NTASKS>1 ) then
!

!----------------------------------------------------------------------
!***  FORECAST TASK 0 SENDS THE DOMAIN SIZE INFORMATION
!***  TO THE FIRST WRITE TASK IN EACH WRITE GROUP BECAUSE 
!***  THE WRITE TASKS NEED TO KNOW THIS TO ASSEMBLE THE
!***  FINAL GRIDDED DATA.
!***  FIRST WE NEED THE VM.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get the Current VM"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_VMGetCurrent(vm=VM                                      &  !<-- The ESMF virtual machine for this group of tasks
                            ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!                            -- IM --
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Forecast task 0 sends
        write(0,*)'in first_pass, IM=',wrt_int_state%IM(1),NCURRENT_GROUP
        DO N=0,NWTPG-1
          CALL MPI_SEND(wrt_int_state%IM                                &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Failed to send IM from fcst task0 to write tasks'
!
        ENDDO
      ENDIF 

      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- Write tasks in this group receive
        CALL MPI_RECV(wrt_int_state%IM                                  &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Recv from fcst 0
                     ,0                                                 &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
        write(0,*)'in first_pass,recv IM=',wrt_int_state%IM(1),NCURRENT_GROUP
!
        IF(IERR/=0)WRITE(0,*)' Write tasks failed to receive IM from fcst task0'
!
      ENDIF
!
!-----------------------------------------------------------------------
!                          -- JM --
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Forecast task 0 sends
        DO N=0,NWTPG-1                                            
          CALL MPI_SEND(wrt_int_state%JM                                &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Failed to send JM from fcst task0 to write tasks'
!
        ENDDO
      ENDIF 
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- Write tasks in this group receive
        CALL MPI_RECV(wrt_int_state%JM                                  &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Recv from fcst 0
                     ,0                                                 &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
!
        IF(IERR/=0)WRITE(0,*)' Write tasks failed to receive JM from fcst task0'
!
      ENDIF
!
!-----------------------------------------------------------------------
!                          -- LM --
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Forecast task 0 sends
        DO N=0,NWTPG-1                                            
          CALL MPI_SEND(wrt_int_state%LM                                &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,0                                               &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          IF(IERR/=0)WRITE(0,*)' Failed to send LM from fcst task0 to write tasks'
!
        ENDDO
      ENDIF 
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- Write tasks in this group receive
        CALL MPI_RECV(wrt_int_state%LM                                  &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,0                                                 &  !<-- Recv from fcst 0
                     ,0                                                 &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
!
        IF(IERR/=0)WRITE(0,*)' Write tasks failed to receive LM from fcst task0'
!
      ENDIF
!
      IM=wrt_int_state%IM(1)
      JM=wrt_int_state%JM(1)
      write(0,*)'in first pass, IM=',IM,'JM=',JM
!
!-----------------------------------------------------------------------
!*** set up the dimension for write pes
!-----------------------------------------------------------------------
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                    !<-- Write tasks in this group receive

       if(.not. allocated(wrt_int_state%nlat_write_task) )then
       allocate(wrt_int_state%jstart_write(NWTPG))
       allocate(wrt_int_state%jend_write(NWTPG))
       allocate(wrt_int_state%nlat_write_task(NUM_PES_FCST))
       allocate(wrt_int_state%nstart_write_task(NUM_PES_FCST))
       endif
!
       nbelt=jm/NWTPG
       nremain=mod(jm,NWTPG)

       DO I=1,NWTPG
         if(mod(i-1,NWTPG)<nremain) then
             wrt_int_state%JSTART_WRITE(i)=(i-1)*(nbelt+1) +1
             wrt_int_state%JEND_WRITE(i)=wrt_int_state%JSTART_WRITE(i)+nbelt
          else
             wrt_int_state%JSTART_WRITE(i)=nremain*(nbelt+1)+(i-1-nremain)*nbelt+1
             wrt_int_state%JEND_WRITE(i)=wrt_int_state%JSTART_WRITE(i)+nbelt-1
          endif
        ENDDO
!
       JSTA_WRITE=wrt_int_state%JSTART_WRITE(mype-lead_write_task+1)
       JEND_WRITE=wrt_int_state%JEND_WRITE(mype-lead_write_task+1)
       if(.not. allocated(wrt_int_state%fcst_lat_on_write_task) )then
         allocate(wrt_int_state%fcst_lat_on_write_task(JSTA_WRITE:JEND_WRITE))
       endif
       write(0,*)'in first pass,jstart_write=',wrt_int_state%JSTART_WRITE, &
         'jend_write=',wrt_int_state%JEND_WRITE
!
      ENDIF
!
!-----------------------------------------------------------------------
!*** get the subdomain information for write tasks
!-----------------------------------------------------------------------
!
     IF(MYPE<NUM_PES_FCST) THEN
!
        DO N=0,NWTPG-1
          CALL MPI_SEND(wrt_int_state%nlat_write_task(N+1)               &  !<-- Send this data
                       ,1                                               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,mype+1001                                       &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
!
          if (wrt_int_state%nlat_write_task(N+1)>0) then
            ISTART=wrt_int_state%nstart_write_task(N+1)
            CALL MPI_SEND(wrt_int_state%fcst_lat_for_write_task(ISTART)             &  !<-- Send this data
                       ,wrt_int_state%nlat_write_task(N+1)               &  !<-- Number of words sent
                       ,MPI_INTEGER                                     &  !<-- Datatype
                       ,N                                               &  !<-- Send to each of the write tasks (local IDs)
                       ,mype+1001                                       &  !<-- An MPI tag
                       ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)            &  !<-- MPI communicator
                       ,IERR)
            write(0,*)'nlat on fcst=',wrt_int_state%nlat_write_task(N+1), &
             'nstart=',wrt_int_state%nstart_write_task(N+1),&
             'fcst_lat=',wrt_int_state%fcst_lat_for_write_task( &
             wrt_int_state%nstart_write_task(N+1):wrt_int_state%nstart_write_task(N+1)+ &
             wrt_int_state%nlat_write_task(N+1)-1)
!
            IF(IERR/=0)WRITE(0,*)' Failed to send IM from fcst task0 to write tasks'
          endif
!
        ENDDO
     ENDIF
!
     IF(MYPE>=LEAD_WRITE_TASK) THEN
       DO I=1,NUM_PES_FCST
           CALL MPI_RECV(wrt_int_state%nlat_write_task(I)                 &  !<-- Recv this data
                     ,1                                                 &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,I-1                                               &  !<-- Recv from fcst 0
                     ,I+1000                                            &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR)
           write(0,*)'wrt task,nlat=',wrt_int_state%nlat_write_task(I)
!
           if(I==1) THEN
             wrt_int_state%nstart_write_task(I)=JSTA_WRITE
           else
             wrt_int_state%nstart_write_task(I)=                        &
                wrt_int_state%nstart_write_task(I-1)+                   &
                wrt_int_state%nlat_write_task(I-1)
           endif
           ISTART=wrt_int_state%nstart_write_task(I)
           if(wrt_int_state%nlat_write_task(I)>0) then
            CALL MPI_RECV(wrt_int_state%fcst_lat_on_write_task(ISTART)  &  !<-- Recv this data
                     ,wrt_int_state%nlat_write_task(I)                  &  !<-- Words received
                     ,MPI_INTEGER                                       &  !<-- Datatype
                     ,I-1                                               &  !<-- Recv from fcst 0
                     ,I+1000                                            &  !<-- An MPI tag
                     ,MPI_COMM_INTER_ARRAY(NCURRENT_GROUP)              &  !<-- MPI communicator
                     ,JSTAT                                             &  !<-- MPI status object
                     ,IERR) 
!
            IF(IERR/=0)WRITE(0,*)' Failed to send IM from fcst task0 to write tasks'

           write(0,*)'wrt task,nstart=',istart,'fcst_lat=', &
             wrt_int_state%fcst_lat_on_write_task(ISTART:  &
             ISTART+wrt_int_state%nlat_write_task(I)-1)
          endif
!
       ENDDO
       write(0,*)'nlat_write_task from fcst=',wrt_int_state%nlat_write_task, &
        'nstart_wrt_task=',wrt_int_state%nstart_write_task, &
        'fcst_lat_on_wrttask=',wrt_int_state%fcst_lat_on_write_task
!       
      ENDIF
!
!-----------------------------------------------------------------------
!#####jw #############     for 1 tasks
!-----------------------------------------------------------------------
!
       else
         IM=wrt_int_state%IM(1)
         JM=wrt_int_state%JM(1)
         write(0,*)'in first pass, IM=',IM,'JM=',JM
!
        endif if1task1
!
!-----------------------------------------------------------------------
!***************** part 2: bundles *************************************
!***********************************************************************
!*** each bundle is corresponding to one file, and 
!*** num_file is the total number of files that need to be opened
!-----------------------------------------------------------------------
!*** get bundle
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  THE NUMBER OF Attributes (FOR SCALARS AND 1D ARRAYS) AND
!***  Fields (FOR GRIDDED 2D ARRAYS) IN THE WRITE COMPONENT'S
!***  IMPORT STATE ARE NOT KNOWN A PRIORI.  IN ORDER TO TRANSFER
!***  THEM TO THE WRITE TASKS, EXTRACT THE NUMBER OF EACH OF
!***  THEM ALONG WITH THEIR NAMES.  THE SCALARS CAN BE LUMPED IN
!***  WITH THE 1D ARRAYS AT THIS POINT.
!
!***  EVEN THOUGH THESE COUNTS ARE JUST SCALAR INTEGERS THEIR
!***  POINTERS WERE ALLOCATED IN WRT_INIT TO LENGTH 1 SINCE THEY
!***  WILL BE USED IN ESMF_Send/Recv WHICH REQUIRE THEM TO BE 
!***  CONTIGUOUS DATA ARRAYS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ALL INTEGER QUANTITIES (AS 1D ARRAYS) AND 1D AND 2D REAL
!***  QUANTITIES WILL BE STRUNG TOGETHER IN SINGLE ARRAYS OF 
!***  EACH PARTICULAR TYPE.  ARRAYS THAT WILL HOLD THE LENGTH OF  
!***  EACH OF THE QUANTITIES IN THESE 'STRINGS' WERE ALLOCATED
!***  IN WRT_INIT.
!-----------------------------------------------------------------------
!
      KOUNT_I1D=0
      KOUNT_R1D=0
      KOUNT_LOG=0
!
      LENGTH_SUM_I1D=0
      LENGTH_SUM_R1D=0
      LENGTH_SUM_LOG=0
!
      NCHAR_I1D=0
      NCHAR_R1D=0
      NCHAR_LOG=0
!
!-----------------------------------------------------------------------
!
      fcst_tasks: IF(MYPE<=LAST_FCST_TASK)THEN                             !<-- Only forecast tasks will extract output information
                                                                           !    from the import state because only they participated
                                                                           !    in filling the import state in the Dynamics/Physics
                                                                           !    components.
!
!-----------------------------------------------------------------------
!***  FIRST FIND THE NUMBER OF Attributes IN THE HISTORY DATA Bundle
!***  IN THE IMPORT STATE AND THEN FIND THEIR NAMES, LENGTHS, AND
!***  DATATYPES.
!***  EXTRACT THE INTEGER AND REAL DATA AND PACK IT INTO INTEGER
!***  AND REAL BUFFERS.  LATER THE BUFFERS WILL BE SENT FROM THE
!***  FORECAST TASKS (THE ONLY ONES WHO CAN SEE THE ORIGINAL
!***  DATA) TO THE WRITE TASKS.
!
!***  THE FACT THAT THE Attribute HISTORY DATA IS BEING COLLECTED
!***  HERE IN A BLOCK THAT EXECUTES ONLY ONCE PER WRITE GROUP
!***  IMPLIES THE ASSUMPTION THAT ONLY THE 2D/3D DATA
!***  ASSOCIATED WITH THE FORECAST GRID CAN CHANGE WITH TIME.
!***  IF ANY SCALAR/1D Attribute DATA CHANGE WITH TIME THEN
!***  THIS MUST BE MOVED OUT OF THIS 'FIRST' BLOCK.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Attribute Count from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_AttributeGet(bundle =FILE_BUNDLE                      &  !<-- The write component's history data Bundle
                              ,count  =NUM_ATTRIB                       &  !<-- # of Attributes in the history data Bundle
                              ,rc     =RC)       
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
        attribute_loop: DO N=1,NUM_ATTRIB                                  !<-- Loop through all the Attributes
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Attribute Names, Datatypes, Lengths"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_AttributeGet(bundle         =FILE_BUNDLE            &  !<-- The write component's history data Bundle
                                ,attributeIndex =N                      &  !<-- Index of each Attribute
                                ,name           =ATTRIB_NAME            &  !<-- Each Attribute's name
                                ,typekind       =DATATYPE               &  !<-- Each Attribute's ESMF Datatype
                                ,count          =LENGTH                 &  !<-- Each Attribute's length
                                ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            write(0,*)'get var ',trim(ATTRIB_NAME),'from bundle,length=', &
            length,'DATATYPE=',DATATYPE,'NBDL=',NBDL
!
!-----------------------------------------------------------------------
!                 -- SCALAR AND 1D INTEGER HISTORY DATA --
!-----------------------------------------------------------------------
!
          IF(DATATYPE==ESMF_TYPEKIND_I4)THEN                               !<-- Extract integer data with rank <2
!
            ALLOCATE(WORK_ARRAY_I1D(LENGTH),stat=RC)                       !<-- This length is from the preceding call 
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Scalar/1-D Integer Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(bundle    =FILE_BUNDLE               &  !<-- The write component's history data Bundle
                                  ,name      =ATTRIB_NAME               &  !<-- Name of the Attribute to extract
                                  ,count     =LENGTH                    &  !<-- Length of Attribute
                                  ,valueList =WORK_ARRAY_I1D            &  !<-- Place the Attribute here
                                  ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_I1D=KOUNT_I1D+1                                          !<-- Count # of integer Attributes
!
            NPOSN_END=KOUNT_I1D*NAME_MAXSTR
            NPOSN_START=NPOSN_END-NAME_MAXSTR+1     
            wrt_int_state%NAMES_I1D_STRING(NBDL)(NPOSN_START:NPOSN_END)=       &
              ATTRIB_NAME(1:NAME_MAXSTR)                                   !<-- Save the 1D integer names
            NCHAR_I1D=NCHAR_I1D+NAME_MAXSTR                                !<-- Save #of characters in all scalar/1D integer names
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces
!
            DO L=1,LENGTH
              wrt_int_state%ALL_DATA_I1D(LENGTH_SUM_I1D+L,NBDL)=WORK_ARRAY_I1D(L)  !<-- String together the integer data
            ENDDO
!
            LENGTH_SUM_I1D=LENGTH_SUM_I1D+LENGTH                           !<-- Total word sum of scalar/1D integer data
            wrt_int_state%LENGTH_DATA_I1D(KOUNT_I1D,NBDL)=LENGTH        !<-- Store length of each individual integer variable
!
            DEALLOCATE(WORK_ARRAY_I1D)
!
!-----------------------------------------------------------------------
!                  -- SCALAR AND 1D REAL HISTORY DATA --
!-----------------------------------------------------------------------
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN                           ! <-- Extract real data with rank <2
!
            ALLOCATE(WORK_ARRAY_R1D(LENGTH),stat=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Scalar/1-D Real Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(bundle    =FILE_BUNDLE               &  !<-- The write component's history data Bundle
                                  ,name      =ATTRIB_NAME               &  !<-- Name of the Attribute to extract
                                  ,count     =LENGTH                    &  !<-- Length of AttributeME
                                  ,valueList =WORK_ARRAY_R1D            &  !<-- Place the Attribute here
                                  ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_R1D=KOUNT_R1D+1                                          !<-- Count # of real Attributes
!
            NPOSN_END=KOUNT_R1D*NAME_MAXSTR
            NPOSN_START=NPOSN_END-NAME_MAXSTR+1     
            wrt_int_state%NAMES_R1D_STRING(NBDL)(NPOSN_START:NPOSN_END)=  &
               ATTRIB_NAME(1:NAME_MAXSTR)                                  !<-- sclar/1D real names
            NCHAR_R1D=NCHAR_R1D+NAME_MAXSTR                                !<-- Save #of characters in all scalar/1D real names
                                                                           !    Note that each name is being given
                                                                           !    NAME_MAXSTR total spaces
!
            DO L=1,LENGTH
              wrt_int_state%ALL_DATA_R1D(LENGTH_SUM_R1D+L,NBDL)=WORK_ARRAY_R1D(L)  !<-- String together the real data
            ENDDO
!
            LENGTH_SUM_R1D=LENGTH_SUM_R1D+LENGTH                           !<-- Total word sum of scalar/1D real data
            wrt_int_state%LENGTH_DATA_R1D(KOUNT_R1D,NBDL)=LENGTH                !<-- Store length of each individual real variable
!
            DEALLOCATE(WORK_ARRAY_R1D)
!
!
!-----------------------------------------------------------------------
!                  -- SCALAR AND 1D REAL8 HISTORY DATA --
!-----------------------------------------------------------------------
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R8)THEN                           ! <-- Extract real data with rank <2
!
            ALLOCATE(WORK_ARRAY_R1D8(LENGTH),stat=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Scalar/1-D Real8 Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(bundle    =FILE_BUNDLE               &  !<-- The write component's history data Bundle
                                  ,name      =ATTRIB_NAME               &  !<-- Name of the Attribute to extract
                                  ,count     =LENGTH                    &  !<-- Length of AttributeME
                                  ,valueList =WORK_ARRAY_R1D8           &  !<-- Place the Attribute here
                                  ,rc=RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_R1D=KOUNT_R1D+1                                          !<-- Count # of real Attributes
!
            NPOSN_END=KOUNT_R1D*NAME_MAXSTR
            NPOSN_START=NPOSN_END-NAME_MAXSTR+1
            wrt_int_state%NAMES_R1D_STRING(NBDL)(NPOSN_START:NPOSN_END)=  &
               ATTRIB_NAME(1:NAME_MAXSTR)                                  !<-- sclar/1D real names
            NCHAR_R1D=NCHAR_R1D+NAME_MAXSTR                                !<-- Save #of characters in all scalar/1D real names
                                                                           !    Note that each name is being given
                                                                           !    NAME_MAXSTR total spaces
!change r8-->r4
            DO L=1,LENGTH
              wrt_int_state%ALL_DATA_R1D(LENGTH_SUM_R1D+L,NBDL)=WORK_ARRAY_R1D8(L)  !<-- String together the real data
            ENDDO
!
            LENGTH_SUM_R1D=LENGTH_SUM_R1D+LENGTH                           !<-- Total word sum of scalar/1D real data
            wrt_int_state%LENGTH_DATA_R1D(KOUNT_R1D,NBDL)=LENGTH                !<-- Store length of each individual real variable
!
            DEALLOCATE(WORK_ARRAY_R1D8)
!
!-----------------------------------------------------------------------
!                          -- LOGICAL DATA --                       
!-----------------------------------------------------------------------
!
!jw datatype== 9 for logical
!          ELSEIF(DATATYPE==ESMF_TYPEKIND_I1)THEN                              ! <-- Extract logical data
           ELSE
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Get Logical Data from History Bundle"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            CALL ESMF_AttributeGet(bundle =FILE_BUNDLE               &  !<-- The write component's history data Bundle
                                  ,name   =ATTRIB_NAME                  &  !<-- Name of the Attribute to extract
                                  ,value  =WORK_LOGICAL                 &  !<-- Place the Attribute here
                                  ,rc     =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            KOUNT_LOG=KOUNT_LOG+1                                          !<-- Count # of logical Attributes
!
            NPOSN_END=KOUNT_LOG*NAME_MAXSTR
            NPOSN_START=NPOSN_END-NAME_MAXSTR+1     
            wrt_int_state%NAMES_LOG_STRING(NBDL)(NPOSN_START:NPOSN_END)= &
              ATTRIB_NAME(1:NAME_MAXSTR)                                   !<-- Save the logical names
            NCHAR_LOG=NCHAR_LOG+NAME_MAXSTR                                !<-- Save #of characters in all logical names
                                                                           !    Note that each name is being given
                                                                           !    EMSF_MAXSTR total spaces
!
            wrt_int_state%ALL_DATA_LOG(KOUNT_LOG,NBDL)=WORK_LOGICAL             !<-- String together the logical data
!
            LENGTH_SUM_LOG=LENGTH_SUM_LOG+1                                !<-- Total length of all logical data variables
!
          ENDIF
!
!-----------------------------------------------------------------------
!
        ENDDO attribute_loop
        write(0,*)'after attribute,KOUNT_I1D=',KOUNT_I1D,'R1D=',KOUNT_R1D, &
         'LOG=',KOUNT_LOG,LENGTH_SUM_LOG,'NBDL=',NBDL
!
!-----------------------------------------------------------------------
!***  INSERT NUMBER AND LENGTHS OF SCALAR/1D INTEGER AND REAL QUANTITIES
!***  AND LOGICALS INTO THE WRITE COMPONENT'S INTERNAL STATE.
!-----------------------------------------------------------------------
!
        wrt_int_state%KOUNT_I1D(NBDL)=KOUNT_I1D
        wrt_int_state%KOUNT_R1D(NBDL)=KOUNT_R1D
        wrt_int_state%KOUNT_LOG(NBDL)=KOUNT_LOG
!
        wrt_int_state%LENGTH_SUM_I1D(NBDL)=LENGTH_SUM_I1D 
        wrt_int_state%LENGTH_SUM_R1D(NBDL)=LENGTH_SUM_R1D
        wrt_int_state%LENGTH_SUM_LOG(NBDL)=LENGTH_SUM_LOG
!
!-----------------------------------------------------------------------
!***  NOW EXTRACT THE NUMBER OF ESMF Fields IN THE HISTORY DATA Bundle
!***  WRITE COMPONENT'S IMPORT STATE ALONG WITH THEIR NAMES.
!***  SAVE THE Field INFORMATION OF THE 2D HISTORY DATA SINCE 
!***  IT WILL BE NEEDED FOR DATA EXTRACTION FROM THE IMPORT STATE
!***  AND THE TRANSFER TO THE WRITE TASKS.
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIND OUT THE NUMBER OF ESMF Fields IN THE HISTORY DATA Bundle.
!***  IT WAS Fields INTO WHICH THE 2D GRIDDED HISTORY DATA WAS PLACED.
!
!***  THIS INFORMATION WILL BE SAVED IN THE INTERNAL STATE
!***  FOR ALL OUTPUT TIMES AND NOT BE RETRIEVED OVER AND OVER AGAIN.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Get Field Count from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_FieldBundleGet(bundle    =FILE_BUNDLE                   &  !<-- The write component's history data Bundle
                                ,fieldCount=wrt_int_state%NCOUNT_FIELDS(NBDL)   &  !<-- Get total # of Fields in the history data Bundle
                                ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  NOW EXTRACT THE NAMES OF ALL THE Fields IN THE BUNDLE.  
!***  ALSO, THE NUMBER OF Field NAMES RETURNED SHOULD EQUAL 
!***  THE FIELD COUNT IN THE PRECEDING CALL.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Extract Field Names from History Bundle"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        allocate(field_name(5000))
        CALL ESMF_FieldBundleGet(bundle    =FILE_BUNDLE              &  !<-- The write component's history data Bundle
                                ,nameList  =FIELD_NAME               & !<-- Array of ESMF Field names in the Bundle
                                ,nameCount =NUM_FIELD_NAMES             &  !<-- Number of Field names in the Bundle
                                ,rc        =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!        write(0,*)'get num fieldname=',NUM_FIELD_NAMES, 'field_name=',  &
!          trim(field_name(1)),'nfield=',wrt_int_state%NCOUNT_FIELDS(NBDL)
!
        IF(NUM_FIELD_NAMES/=wrt_int_state%NCOUNT_FIELDS(nbdl))THEN
          WRITE(0,*)' WARNING: Number of Fields in Bundle of history'   &
                   ,' output does not equal the number of Field names'
          WRITE(0,*)' They are ',NUM_FIELD_NAMES,' and '                &
                   ,wrt_int_state%NCOUNT_FIELDS(1),', respectively'
        ENDIF
!jw        wrt_int_state%FIELD_NAME(1:NUM_FIELD_NAMES,NBDL)=field_name(1:NUM_FIELD_NAMES)
!jw        deallocate(field_name)
!
!-----------------------------------------------------------------------
!***  DO A PRELIMINARY EXTRACTION OF THE Fields THEMSELVES IN ORDER TO
!***  COUNT THE NUMBER OF REAL AND INTEGER 2D ARRAYS.  
!-----------------------------------------------------------------------
!
        KOUNT_R2D=0
        KOUNT_I2D=0
        NCHAR_I2D(1)=0
        NCHAR_R2D(1)=0
!
        DO N=1,wrt_int_state%NCOUNT_FIELDS(NBDL)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Extract Fields from History Bundle for Counting"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldBundleGet(bundle=FILE_BUNDLE                &  !<-- The write component's history data Bundle
                                  ,name  =FIELD_NAME(N)              &  !<-- The ESMF Field's name
                                  ,field =FIELD_WORK1                   &  !<-- The ESMF Field taken from the Bundle
                                  ,rc    =RC)
         wrt_int_state%FIELD_NAME(N,NBDL)=field_name(N)
!         write(0,*)'fieldname=',wrt_int_state%FIELD_NAME(N,NBDL)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Get Datatype of Fields for Counting Real/Integer"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_FieldGet(field=FIELD_WORK1                          &  !<-- The ESMF 2D Field
                            ,typekind=DATATYPE                          &  !<-- The Field's ESMF Datatype
                            ,rc   =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(DATATYPE==ESMF_TYPEKIND_I4)THEN
            KOUNT_I2D=KOUNT_I2D+1                                          !<-- Add up the total number of integer 2D Fields
!
            IF(KOUNT_I2D>MAX_DATA_I2D)THEN
              WRITE(0,*)' FATAL: YOU HAVE EXCEEDED MAX NUMBER OF INTEGER 2D FIELDS FOR OUTPUT'
              WRITE(0,*)' YOU MUST INCREASE VALUE OF MAX_DATA_I2D WHICH NOW EQUALS ',MAX_DATA_I2D
            ENDIF
!
            NPOSN_END  =KOUNT_I2D*NAME_MAXSTR
            NPOSN_START=NPOSN_END-NAME_MAXSTR+1
            wrt_int_state%NAMES_I2D_STRING(NBDL)(NPOSN_START:NPOSN_END)=  &
              wrt_int_state%FIELD_NAME(N,NBDL)(1:NAME_MAXSTR) !<-- Save the 2D integer Field names 
                                                                                              !<-- in one long string
            NCHAR_I2D(1)=NCHAR_I2D(1)+NAME_MAXSTR
!
          ELSEIF(DATATYPE==ESMF_TYPEKIND_R4)THEN
            KOUNT_R2D=KOUNT_R2D+1                                          !<-- Add up the total number of real 2D Fields
!
            IF(KOUNT_R2D>MAX_DATA_R2D)THEN
              WRITE(0,*)' FATAL: YOU HAVE EXCEEDED MAX NUMBER OF REAL 2D FIELDS FOR OUTPUT'
              WRITE(0,*)' KOUNT_R2D=',KOUNT_R2D
              WRITE(0,*)' YOU MUST INCREASE VALUE OF MAX_DATA_R2D WHICH NOW EQUALS ',MAX_DATA_R2D
            ENDIF
!
            NPOSN_END  =KOUNT_R2D*NAME_MAXSTR
            NPOSN_START=NPOSN_END-NAME_MAXSTR+1
            wrt_int_state%NAMES_R2D_STRING(NBDL)(NPOSN_START:NPOSN_END)=  &
              wrt_int_state%FIELD_NAME(N,NBDL)(1:NAME_MAXSTR) !<-- Save the 2D real Field names 
                                                                                              !<-- in one long string
            NCHAR_R2D(1)=NCHAR_R2D(1)+NAME_MAXSTR
!
          ENDIF
!
        ENDDO
        deallocate(field_name)
!
        wrt_int_state%KOUNT_R2D(NBDL)=KOUNT_R2D
        wrt_int_state%KOUNT_I2D(NBDL)=KOUNT_I2D
        write(0,*)'kount_r2d=',wrt_int_state%KOUNT_R2D(NBDL),'kount_i2d=', &
          wrt_int_state%KOUNT_I2D(NBDL)
!
      ENDIF fcst_tasks
!-----------------------------------------------------------------------
!***  IF THERE ARE NO QUANTITIES SPECIFIED FOR HISTORY OUTPUT,
!***  FORECAST TASK 0 WILL INFORM THE WRITE TASKS AND THEN
!***  EVERYONE WILL RETURN.
!-----------------------------------------------------------------------
!
      NO_FIELDS(1)=ESMF_FALSE
!
!
!-----------------------------------------------------------------------
!#####jw #############     for quilting tasks
!-----------------------------------------------------------------------
!
        if1task2:  if(ntasks>1) then
!##


      IF(MYPE==0)THEN
        IF(wrt_int_state%NCOUNT_FIELDS(1)==0)NO_FIELDS(1)=ESMF_TRUE        !<-- Reset flag saying there are no history quantities
        LAST_WRITE_TASK=NTASKS-1                                           !<-- The last write task in this group
!
        DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK                               !<-- Loop through all the write tasks in the write group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          MESSAGE_CHECK="Fcst Task0 Informs All That There Are No Fields"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                          ,sendData=NO_FIELDS                           &  !<-- Send this data
                          ,count   =1                                   &  !<-- Words sent
                          ,dst     =N                                   &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ENDDO
      ENDIF
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- All write tasks in this group receive
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Told By Fcst Task0 There Are No Fields"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=NO_FIELDS                             &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!
      else   !for 1 tasks
!
        IF(wrt_int_state%NCOUNT_FIELDS(1)==0)NO_FIELDS(1)=ESMF_TRUE        !<-- Reset flag saying there are no history quantities

      endif  if1task2
!-----------------------------------------------------------------------
!
      IF(NO_FIELDS(1)==ESMF_TRUE)THEN
        IF(MYPE==0)THEN
          WRITE(6,*)'WARNING: No Import ESMF quantities for the Write Component'
          WRITE(0,*)'WARNING: No Import ESMF quantities for the Write Component'
        ENDIF
!
        RETURN                                                             !<-- All tasks return if there is no history output
!
      ENDIF
!
!
!-----------------------------------------------------------------------
!#####jw #############     for 1 tasks
!-----------------------------------------------------------------------
!
        if1task3:   if(ntasks >1) then
!-----------------------------------------------------------------------
!***  FORECAST TASK 0 SENDS ALL THE WRITE TASKS THE NUMBER OF
!***  REAL AND INTEGER 2D GRIDDED QUANTITIES PLUS ALL OF THE
!***  LOCAL HORIZONTAL DOMAIN LIMITS IN PREPARATION FOR THE
!***  WRITE TASKS' RECEIVING AND ASSEMBLING THE LOCAL HISTORY
!***  DATA THEY RECEIVE FROM THE FORECAST TASKS.
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Forecast task 0 sends
!
        LAST_WRITE_TASK=NTASKS-1                                           !<-- The last write task in this group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Write Tasks Info for Quilting"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK                               !<-- Loop through all the write tasks in the write group
!
          ITMP(1)=wrt_int_state%NCOUNT_FIELDS(NBDL)
          CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                          ,sendData=ITMP                                &  !<-- Send this data
                          ,count   =1                                   &  !<-- Words sent
                          ,dst     =N                                   &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          ITMP(1)=wrt_int_state%KOUNT_R2D(NBDL)    
          CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                          ,sendData=ITMP                                &  !<-- Send this data
                          ,count   =1                                   &  !<-- Words sent
                          ,dst     =N                                   &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
          ITMP(1)=wrt_int_state%KOUNT_I2D(NBDL)    
          CALL ESMF_VMSend(vm      =VM                                  &  !<-- ESMF Virtual Machine
                          ,sendData=ITMP                                &  !<-- Send this data
                          ,count   =1                                   &  !<-- Words sent
                          ,dst     =N                                   &  !<-- Receiving task in active write group
                          ,rc      =RC)
!
        ENDDO
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!
!-----------------------------------------------------------------------
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- All write tasks in this group receive
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Quilting Info From Fcst Task0"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=ITMP                                  &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%NCOUNT_FIELDS(NBDL)=ITMP(1)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=ITMP                                  &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_R2D(NBDL)=ITMP(1)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=ITMP                                  &  !<-- Recv this data
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_I2D(NBDL)=ITMP(1)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF
!      write(0,*)'all wrt comp,inpes=',wrt_int_state%inpes,'jnpes=',wrt_int_state%jnpes,&
!       'nfield=',wrt_int_state%NCOUNT_FIELDS(NBDL),'KOUNT_R2D(NBDL)=', &
!       wrt_int_state%KOUNT_R2D(NBDL),'wrt_int_state%KOUNT_I2D=', &
!       wrt_int_state%KOUNT_I2D(NBDL)
!
!-----------------------------------------------------------------------
!***  FORECAST TASK 0 SENDS THE 2D DATA NAMES TO ALL LEAD WRITE TASKS.
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN                                                      !<-- Fcst task0 alone can send write task preliminary info
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Write Tasks 2D Data Names"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
       DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK
!
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=NCHAR_I2D                             &  !<-- Send total length of the names of 2D integer data
                        ,count   =1                                     &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
!jw                        ,dst     =LEAD_WRITE_TASK                       &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        if(NCHAR_I2D(1)>0) then
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%NAMES_I2D_STRING(nbdl) &  !<-- Send names of 2D integer history variables
                        ,count   =NCHAR_I2D(1)                          &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
        endif
!
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=NCHAR_R2D                             &  !<-- Send total length of the names of 2D real data
                        ,count   =1                                     &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        if(NCHAR_R2D(1)>0) then
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=wrt_int_state%NAMES_R2D_STRING(NBDL)  &  !<-- Send names of 2D real history variables
                        ,count   =NCHAR_R2D(1)                          &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
         endif
!
          ENDDO
!
!jw      ELSEIF(MYPE==NTASKS-NWTPG)THEN                                       !<-- 1st write task receives 2D preliminary info
      ELSEIF(MYPE>=LEAD_WRITE_TASK) then                                    !<-- 1st write task receives 2D preliminary info
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=NCHAR_I2D                             &  !<-- Recv total length of the names of 2D integer data
                        ,count   =1                                     &  !<-- Words sent
                        ,src     =0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
!
        if(NCHAR_I2D(1)>0) then
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%NAMES_I2D_STRING(NBDL)        &  !<-- Recv names of 2D integer history variables
                        ,count   =NCHAR_I2D(1)                          &  !<-- Words sent
                        ,src     =0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
         endif
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=NCHAR_R2D                             &  !<-- Recv total length of the names of 2D gridded data
                        ,count   =1                                     &  !<-- Words sent
                        ,src     =0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
!
        if(NCHAR_R2D(1)>0) then
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=wrt_int_state%NAMES_R2D_STRING(NBDL)        &  !<-- Recv names of 2D real history variables
                        ,count   =NCHAR_R2D(1)                          &  !<-- Words sent
                        ,src     =0                                     &  !<-- Sending task (fcst task 0)
                        ,rc      =RC)
         endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDIF  !write task
      write(0,*)'write task get names for names_R2D_string=', &
        wrt_int_state%NAMES_R2D_STRING(NBDL)(1:5)
!
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK MUST KNOW THE IDs OF THE FORECAST TASKS
!***  FROM WHICH IT WILL RECEIVE 2D GRIDDED HISTORY DATA.
!########################################################################
!
      IF(MYPE>=LEAD_WRITE_TASK)THEN                                        !<-- The write tasks
!
        MYPE_LOCAL=MOD(MYPE-NUM_PES_FCST,NWTPG)
!
!-----------------------------------------------------------------------
!***  EACH WRITE TASK ALSO MUST KNOW THE NORTH-SOUTH EXTENT OF THE
!***  FULL 2D DOMAIN THAT IT WILL HANDLE.  THIS IS DETERMINED BY
!***  THE COVERAGE OF THE FCST TASKS THAT SEND TO IT.
!-----------------------------------------------------------------------
!
        JSTA_WRITE=wrt_int_state%JSTART_WRITE(mype_local+1)  !<-- JTS of 1st fcst task that sends to this write task
        JEND_WRITE=wrt_int_state%JEND_WRITE(mype_local+1)  !<-- JTE of last fcst task that sends to this write task
         write(0,*)'jsta_write=',JSTA_WRITE,'jend_write=',&
          JEND_WRITE
!
!-----------------------------------------------------------------------
!***  NOW EACH WRITE TASK ALLOCATES ITS OWN SECTION OF THE 2D DOMAIN
!***  FOR ALL THE 2D VARIABLES IT WILL RECEIVE AND ITS 1D EQUIVALENT
!***  USED TO TRANSFER THE DATA TO THE LEAD WRITE TASK.
!-----------------------------------------------------------------------
!
       if(.not.allocated(wrt_int_state%BUFF_INT)) then
        LENGTH=IM*JM
        ALLOCATE(wrt_int_state%BUFF_INT(LENGTH))
        ALLOCATE(wrt_int_state%BUFF_INT_TMP(LENGTH))
        ALLOCATE(wrt_int_state%BUFF_REAL(LENGTH))
        ALLOCATE(wrt_int_state%BUFF_REAL_TMP(LENGTH))
       endif
!
       ENDIF  !MYPE>=LEAD_WRITE_TASK
!
!-----------------------------------------------------------------------
!#####jw #############     for 1 tasks
!-----------------------------------------------------------------------
!
        endif if1task3
!
!-----------------------------------------------------------------------
!***  THE LEAD WRITE TASK ALLOCATES ITS WORKING ARRAYS INTO WHICH
!***  IT WILL ASSEMBLE EACH INDIVIDUAL 2D FIELD THAT WILL BE
!***  WRITTEN TO THE HISTORY FILES.
!-----------------------------------------------------------------------
!
!jw      IF(MYPE==LEAD_WRITE_TASK)THEN
      IF((MYPE>=LEAD_WRITE_TASK .and. MYPE_LOCAL==0) .or. ntasks==1 )THEN
       if(.not.allocated(wrt_int_state%OUTPUT_ARRAY_I2D)) then
        ALLOCATE(wrt_int_state%OUTPUT_ARRAY_I2D(1:IM,1:JM))
        ALLOCATE(wrt_int_state%OUTPUT_ARRAY_R2D(1:IM,1:JM))
       endif
      ENDIF
      write(0,*)'in first pass, output_array_r2d,size=',size(wrt_int_state%OUTPUT_ARRAY_R2D,1), &
       size(wrt_int_state%OUTPUT_ARRAY_R2D,2),im,jm
! 
!-----------------------------------------------------------------------
!#####jw #############     for 1 tasks
!-----------------------------------------------------------------------
!
        if1task4:  if(ntasks >1) then
!
!-----------------------------------------------------------------------
!***  SINCE ALL SCALAR/1D DATA IS IDENTICAL ON ALL FORECAST TASKS,
!***  TASK 0 ALONE CAN SEND THE INFORMATION TO THE LEAD WRITE TASK
!***  THAT WILL LATER WRITE IT TO THE HISTORY FILE.
!--------------------------------------------------------------------
!
!--------------------------------------------------------------------
      task_0_sends: IF(MYPE==0)THEN                                      !<-- Forecast task 0 sends
!--------------------------------------------------------------------
!
      write(0,*)'before task 0 send 1d vars, all processors,lead_write_task=', &
       lead_write_task,'last_write_task=',last_write_task
!------------------------------------------------
!***  SEND SCALAR/1D INTEGER HISTORY INFORMATION. to all the write tasks
!------------------------------------------------
!
       DO N=LEAD_WRITE_TASK,LAST_WRITE_TASK

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Scalar/1D Integer History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ITMP(1)=wrt_int_state%KOUNT_I1D(NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=ITMP                               &  !<-- Send # of scalar/1D integer history variables
                        ,count   =1                                     &  !<-- Words sent
!jw                        ,dst     =LEAD_WRITE_TASK                       &  !<-- Receiving task (1st write task in group)
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
!if there are any I 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_I1D(NBDL)>0 )then
!
        ITMP(1)=wrt_int_state%LENGTH_SUM_I1D(NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=ITMP                                  &  !<-- Send length of string of all such integer history variables
                        ,count   =1                                     &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        ITMP(1:wrt_int_state%KOUNT_I1D(NBDL))=                          &
         wrt_int_state%LENGTH_DATA_I1D(1:wrt_int_state%KOUNT_I1D(NBDL),NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=ITMP                                  &  !<-- Send lengths of each scalar/1D integer history variable
                        ,count   =wrt_int_state%KOUNT_I1D(NBDL)         &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        NAMETMP(1:NCHAR_I1D)=wrt_int_state%NAMES_I1D_STRING(nbdl)(1:NCHAR_I1D)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=NAMETMP                               &  !<-- Send names of each scalar/1D integer history variable
                        ,count   =NCHAR_I1D                                      &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        ITMP(1:wrt_int_state%LENGTH_SUM_I1D(NBDL))=                     &
           wrt_int_state%ALL_DATA_I1D(1:wrt_int_state%LENGTH_SUM_I1D(NBDL),NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=ITMP                                  &  !<-- Send the full string of all scalar/1D integer history data
                        ,count   =wrt_int_state%LENGTH_SUM_I1D(NBDL)       &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        write(0,*)'I1D,kount=',wrt_int_state%KOUNT_I1D(NBDL), 'length=', &
         wrt_int_state%LENGTH_DATA_I1D(1:wrt_int_state%KOUNT_I1D(NBDL),NBDL), &
         'I1D value=',wrt_int_state%ALL_DATA_I1D(1:wrt_int_state%LENGTH_SUM_I1D(NBDL),NBDL) 
         endif
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!---------------------------------------------
!***  SEND SCALAR/1D REAL HISTORY INFORMATION.
!---------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Scalar/1D Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ITMP(1)=wrt_int_state%KOUNT_R1D(NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=ITMP                                  &  !<-- Send # of scalar/1D real history variables
                        ,count   =1                                     &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
!
!if there are any R 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_R1D(NBDL)>0 )then
!
        ITMP(1)=wrt_int_state%LENGTH_SUM_R1D(NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=ITMP                                  &  !<-- Send length of string of all such real history variables
                        ,count   =1                                     &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        ITMP(1:wrt_int_state%KOUNT_R1D(NBDL))=wrt_int_state%LENGTH_DATA_R1D &
          (1:wrt_int_state%KOUNT_R1D(NBDL),NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=ITMP                                  &!<-- Send lengths of each scalar/1D real history variable
                        ,count   =wrt_int_state%KOUNT_R1D(NBDL)         &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        NAMETMP(1:NCHAR_R1D)=wrt_int_state%NAMES_R1D_STRING(NBDL)(1:NCHAR_R1D)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=NAMETMP                               &  !<-- Send names of each scalar/1D real history variable
                        ,count   =NCHAR_R1D                             &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        RTMP(1:wrt_int_state%LENGTH_SUM_R1D(NBDL))=wrt_int_state%ALL_DATA_R1D( &
             1:wrt_int_state%LENGTH_SUM_R1D(NBDL),NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=RTMP                                  &  !<-- Send the full string of all scalar/1D real history data
                        ,count   =wrt_int_state%LENGTH_SUM_R1D(NBDL)       &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
        write(0,*)'R1D,kount=',wrt_int_state%KOUNT_R1D(NBDL), 'length=', &
         wrt_int_state%LENGTH_DATA_R1D(1:wrt_int_state%KOUNT_R1D(NBDL),NBDL)
!
         ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!--------------------------------------
!***  SEND LOGICAL HISTORY INFORMATION.
!--------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Fcst Task0 Sends Logical History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        ITMP(1)=wrt_int_state%KOUNT_LOG(NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=ITMP                                  &  !<-- Send # of logical history variables
                        ,count   =1                                     &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
!
!if there are any R 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_LOG(NBDL)>0 )then
!
        ITMP(1)=wrt_int_state%LENGTH_SUM_LOG(NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=ITMP                                  &  !<-- Send length of string of all logical variables
                        ,count   =1                                     &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        NAMETMP(1:NCHAR_R1D)=wrt_int_state%NAMES_LOG_STRING(NBDL)(1:NCHAR_LOG)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=NAMETMP                               & !<-- Send names of each logical history variable
                        ,count   =NCHAR_LOG                                      &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        write(0,*)'length_sum_log=',wrt_int_state%LENGTH_SUM_LOG(NBDL)
        DO I=1,wrt_int_state%LENGTH_SUM_LOG(NBDL)
          if(wrt_int_state%ALL_DATA_LOG(i,NBDL)==ESMF_TRUE) then
            LTMP(I)=ESMF_TRUE
        write(0,*)'set ltmp for log, for true, I=',I
          else
            LTMP(I)=ESMF_FALSE
        write(0,*)'set ltmp for log, for false, I=',I
          endif
        ENDDO

!        LTMP(1:wrt_int_state%LENGTH_SUM_LOG(NBDL))=                     &
!          wrt_int_state%ALL_DATA_LOG(1:wrt_int_state%LENGTH_SUM_LOG(NBDL),NBDL)
        CALL ESMF_VMSend(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,sendData=LTMP                                  &  !<-- Send the full string of all logical history data
                        ,count   =wrt_int_state%LENGTH_SUM_LOG(NBDL)    &  !<-- Words sent
                        ,dst     =N                                     &  !<-- Receiving task (1st write task in group)
                        ,rc      =RC)
!
        ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      ENDDO
!
!-----------------------------------------------------------------------
      ENDIF task_0_sends
!-----------------------------------------------------------------------
!
      write(0,*)'before write tasks get 1d vars, all processors'
!-----------------------------------------------------------------------
!
!jw      write_task_recvs: IF(MYPE==LEAD_WRITE_TASK)THEN                      !<-- 1st write task in this group receives
      write_task_recvs: IF(MYPE>=LEAD_WRITE_TASK)THEN                      !<-- 1st write task in this group receives
                                                                           !    all of the data just sent to it by
                                                                           !    fcst task 0
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RECEIVE SCALAR/1D INTEGER HISTORY INFORMATION
!***  FROM FORECAST TASK 0.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Scalar/1D Integer History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=ITMP                                  &  !<-- Recv # of integer history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_I1D(NBDL)=ITMP(1)
        write(0,*)'after vm recv. in first_pass,I1D,itmp=',itmp(1)
        write(0,*)'size(length)=',size(wrt_int_state%LENGTH_SUM_I1D),            &
         size( wrt_int_state%LENGTH_DATA_I1D,1)                        
!
!if there are any R 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_I1D(NBDL)>0 )then
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=itmp                                  &  !<-- Recv length of string of all integer history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
         wrt_int_state%LENGTH_SUM_I1D(NBDL)=itmp(1)
        write(0,*)'after vm recv. in first_pass,I1D,length_sum,itmp=',itmp(1)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=itmp                                  &  !<-- Recv lengths of each integer history variable
                        ,count   =wrt_int_state%KOUNT_I1D(NBDL)         &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%LENGTH_DATA_I1D(1:wrt_int_state%KOUNT_I1D(NBDL),NBDL)= &
          itmp(1:wrt_int_state%KOUNT_I1D(NBDL))
!
        NCHAR_I1D=wrt_int_state%KOUNT_I1D(NBDL)*NAME_MAXSTR
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=NAMETMP                               &  !<-- Recv names of integer history variables
                        ,count   =NCHAR_I1D                             &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%NAMES_I1D_STRING(NBDL)=trim(NAMETMP)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=ITMP                                  &  !<-- Recv the string of integer history data
                        ,count   =wrt_int_state%LENGTH_SUM_I1D(NBDL)    &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%ALL_DATA_I1D(1:wrt_int_state%LENGTH_SUM_I1D(NBDL),NBDL)= &
          ITMP(1:wrt_int_state%LENGTH_SUM_I1D(NBDL))
!
        write(0,*)'after vm recv. in first_pass,I1D,length'
        ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RECEIVE SCALAR/1D REAL HISTORY INFORMATION.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Scalar/1D Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=itmp                                  &  !<-- Recv # of scalar/1D real history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_R1D(NBDL)=itmp(1)
        write(0,*)'after vm recv. in first_pass,R1D,itmp=',itmp(1)
!
!if there are any R 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_R1D(NBDL)>0 )then
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=itmp                                  &  !<-- Recv length of string of all such real history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%LENGTH_SUM_R1D(NBDL)=itmp(1)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=itmp                                  &  !<-- Recv lengths of each scalar/1D real history variable
                        ,count   =wrt_int_state%KOUNT_R1D(NBDL)         &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%LENGTH_DATA_R1D(1:wrt_int_state%KOUNT_R1D(NBDL),NBDL)= &
          itmp(1:wrt_int_state%KOUNT_R1D(NBDL))
!
        NCHAR_R1D=wrt_int_state%KOUNT_R1D(NBDL)*ESMF_MAXSTR
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=nametmp                               &  !<-- Recv names of scalar/1D real history variables
                        ,count   =NCHAR_R1D                             &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%NAMES_R1D_STRING(NBDL)=trim(nametmp)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=rtmp                                  &  !<-- Recv the string of all scalar/1D real history data
                        ,count   =wrt_int_state%LENGTH_SUM_R1D(NBDL)       &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%ALL_DATA_R1D(1:wrt_int_state%LENGTH_SUM_R1D(NBDL),NBDL)= &
          rtmp(1:wrt_int_state%LENGTH_SUM_R1D(NBDL))
!
        write(0,*)'after vm recv. in first_pass,R1D length'
        ENDIF
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  RECEIVE LOGICAL HISTORY INFORMATION.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        MESSAGE_CHECK="Write Tasks Recv Logical Real History Data"
!       CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=itmp                                  &  !<-- Recv # of logical history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%KOUNT_LOG(NBDL)=itmp(1)
        write(0,*)'after vm recv. in first_pass,LOG,itmp=',itmp(1)
!
!if there are any R 1D data
!---------------------------------------------------------------------
        if(wrt_int_state%KOUNT_LOG(NBDL)>0 )then

!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=itmp                                  &  !<-- Recv length of string of all logical history variables
                        ,count   =1                                     &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%LENGTH_SUM_LOG(NBDL)=itmp(1)
!
        NCHAR_LOG=wrt_int_state%KOUNT_LOG(NBDL)*NAME_MAXSTR
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=NAMETMP                               & !<-- Recv names of logical history variables
                        ,count   =NCHAR_LOG                             &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%NAMES_LOG_STRING(NBDL)=trim(NAMETMP)
!
        CALL ESMF_VMRecv(vm      =VM                                    &  !<-- ESMF Virtual Machine
                        ,recvData=LTMP                                  &  !<-- Recv the string of all logical history data
                        ,count   =wrt_int_state%LENGTH_SUM_LOG(NBDL)       &  !<-- Words received
                        ,src     =0                                     &  !<-- Sending task (forecast task 0)
                        ,rc      =RC)
        wrt_int_state%ALL_DATA_LOG(1:wrt_int_state%LENGTH_SUM_LOG(NBDL),NBDL)= &
         LTMP(1:wrt_int_state%LENGTH_SUM_LOG(NBDL))

        write(0,*)'after vm recv. in first_pass,LOG,length'
        ENDIF

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        CALL ERR_MSG(RC,MESSAGE_CHECK,RC_WRT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!
      ENDIF write_task_recvs
      write(0,*)'before first_pass_gfs end,KOUNT_I1D=',wrt_int_state%KOUNT_I1D(NBDL), &
       'I1D length=',wrt_int_state%LENGTH_DATA_I1D(1:wrt_int_state%KOUNT_I1D(NBDL),NBDL), &
       'I1D value=',wrt_int_state%ALL_DATA_I1D(1:wrt_int_state%KOUNT_I1D(NBDL),NBDL), &
       'KOUNT_R1D=',wrt_int_state%KOUNT_R1D(NBDL),'KOUNT_LOG=', &
       wrt_int_state%KOUNT_LOG(NBDL)
!
!-----------------------------------------------------------------------
!#####jw #############     for 1 tasks
!-----------------------------------------------------------------------
!
      endif  if1task4
!
!-----------------------------------------------------------------------
!
      DEALLOCATE(INPES,JNPES)
      DEALLOCATE(ITMP)
      DEALLOCATE(RTMP)
      DEALLOCATE(LTMP)
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE FIRST_PASS_GFS
!
!-----------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_INIT_GFS(ATM_GRID_COMP,wrt_comps,imp_state_write, &
        exp_state_write,CLOCK_ATM,WRITE_GROUP_READY_TO_GO)
!jw        exp_state_write,CLOCK_ATM,alarm_output,WRITE_GROUP_READY_TO_GO)

      use module_gfs_mpi_def, only: num_pes_fcst,write_tasks_per_group,    &
          write_groups,petlist_write
! 
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEP OF THE WRITE COMPONENTS.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp),INTENT(INOUT)      :: ATM_GRID_COMP             !<-- The ATM gridded component
      TYPE(ESMF_GRIDComp),INTENT(INOUT)      :: wrt_comps(:)             !<-- The ATM Internal State
      TYPE(ESMF_STATE),INTENT(INOUT)         :: imp_state_write             !<-- The ATM Internal State
      TYPE(ESMF_STATE),INTENT(INOUT)         :: exp_state_write             !<-- The ATM Internal State
      TYPE(ESMF_Clock),INTENT(INOUT)         :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
!jw      TYPE(ESMF_Alarm),INTENT(INOUT)         :: alarm_output
      INTEGER, INTENT(INOUT)                 :: WRITE_GROUP_READY_TO_GO
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config)      :: CF                                        !<-- The config object
      TYPE(ESMF_VM)          :: VM                                        !<-- The ESMF virtual machine.
      TYPE(ESMF_Time)        :: CURRTIME
      TYPE(ESMF_TimeInterval):: TIMEINTERVAL_OUTPUT
!
      INTEGER :: I,INPES,J,JNPES,RC,RC_INIT,NFHOUT
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  RETRIEVE THE CONFIG OBJECT CF FROM THE ATM GRIDDED COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="Get Config Object from ATM Component in Write Init"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_GridCompGet(gridcomp=ATM_GRID_COMP                      &  !<-- The ATM gridded component
                           ,config  =CF                                 &  !<-- The config object (~namelist)
                           ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
!-----------------------------------------------------------------------
!***  SET UP alarm_output
!-----------------------------------------------------------------------
!
!jw      CALL ESMF_ConfigGetAttribute(config=CF                            &  !<-- The ESMF configure object
!jw                                  ,value =nfhout                         &  !<-- # of fcst tasks in I direction
!jw                                  ,label ='nfhout:'                      &  !<-- Give the value of this label to INPES
!jw                                  ,rc    =RC)
!jw      write(0,*)'nfhout=',nfhout 
!
!jw      CALL ESMF_TimeIntervalSet(timeinterval=TIMEINTERVAL_OUTPUT           &  !<-- Time interval between
!jw                               ,h           =nfhout                          &  !<-- Hours between history
!jw                               ,rc          =RC)
!

!jw      CALL ESMF_ClockGet(Clock   =CLOCK_ATM                             &  !<-- Time interval between
!jw                        ,currTime=CURRTIME                              &
!jw                        ,rc      =RC)
!
!jw      write(0,*)'now print out clock infomation'
!jw      call ESMF_ClockPrint(clock=CLOCK_ATM,rc=rc)
!
!jw      ALARM_OUTPUT=ESMF_AlarmCreate(name          ='ALARM_OUTPUT'      &
!jw                                    ,clock            =CLOCK_ATM            &  !<-- ATM Clock
!jw                                    ,ringTime         =CURRTIME             &  !<-- Forecast/Restart start time (ESMF)
!jw                                    ,ringInterval     =TIMEINTERVAL_OUTPUT  &  !<-- Time interval between
!jw                                    ,ringTimeStepCount=1                    &  !<-- The Alarm rings for this many timesteps
!jw                                    ,sticky           =.false.              &  !<-- Alarm does not ring until turned off
!jw                                    ,rc               =RC)
!
!jw      write(0,*)'now print out clock infomation 2'
!jw      call ESMF_ClockPrint(clock=CLOCK_ATM,rc=rc)
!jw      write(0,*)'in wrt_init_gfs, create alarm_output,RC=',RC
!
!-----------------------------------------------------------------------
!***  EXECUTE THE INITIALIZE STEP FOR THE WRITE COMPONENTS.
!***  THESE ARE THE INITIALIZE SUBROUTINES SPECIFIED IN THE
!***  REGISTER ROUTINES CALLED IN ESMF_GridCompSetServices.
!-----------------------------------------------------------------------
!
      DO J=1,WRITE_GROUPS
!
!!!!    CALL ESMF_VMBarrier(VM,rc=RC)    ! Insert barrier since fcst tasks are involved in each iteration of write groups
!
        DO I=1,NUM_PES_FCST+WRITE_TASKS_PER_GROUP
          IF(MYPE==PETLIST_WRITE(I,J))THEN                   !<--  Forecast tasks plus the Write tasks in each write group
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            MESSAGE_CHECK="Execute Initialize Step of Write Component"
!           CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            write(0,*)'before grid comp init,j=',j,'pelist_wrt=',PETLIST_WRITE
            N_GROUP=J
            CALL ESMF_GridCompInitialize(WRT_COMPS(J)                 &  !<-- The Write gridded components
                                        ,importstate=IMP_STATE_WRITE  &  !<-- The Write import state
                                        ,exportstate=EXP_STATE_WRITE  &  !<-- The Write export state
                                        ,clock      =CLOCK_ATM                      &  !<-- The ESMF clock of the ATM component
                                        ,phase      =ESMF_SINGLEPHASE               &
                                        ,rc         =RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            CALL ERR_MSG(RC,MESSAGE_CHECK,RC_INIT)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
            write(0,*)'after wrt_comp, initial'
          ENDIF
        ENDDO
      ENDDO
!
!-----------------------------------------------------------------------
!***  SET THE FIRST WRITE GROUP AS THE FIRST ONE TO ACT.
!-----------------------------------------------------------------------
!
      WRITE_GROUP_READY_TO_GO=1
      N_GROUP=WRITE_GROUP_READY_TO_GO
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_INIT_GFS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_ASYNC_GFS(WRT_COMPS,EXP_STATE,IMP_STATE_WRITE     &
                            ,EXP_STATE_WRITE,CLOCK_ATM,MYPE              &
                            ,WRITE_GROUP_READY_TO_GO)

      use module_gfs_mpi_def, only: num_pes_fcst,write_groups,           &
          write_tasks_per_group,PETLIST_WRITE
!
!-----------------------------------------------------------------------
!***  WRITE OUT A HISTORY FILE USING THE ASYNCHRONOUS QUILTING.
!-----------------------------------------------------------------------
!
      TYPE(ESMF_GridComp)     ,INTENT(INOUT) :: WRT_COMPS(:)               !<-- The write_comp
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: EXP_STATE                  !<-- The export state dyn
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: IMP_STATE_WRITE            !<-- The import state of write
      TYPE(ESMF_STATE)        ,INTENT(INOUT) :: EXP_STATE_WRITE            !<-- The export state of write
      TYPE(ESMF_Clock)        ,INTENT(INOUT) :: CLOCK_ATM                 !<-- The ATM Component's ESMF Clock
      INTEGER,INTENT(IN) :: MYPE
      INTEGER,INTENT(INOUT) :: WRITE_GROUP_READY_TO_GO
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      TYPE(ESMF_Config) :: CF                                             !<-- The configure object (~namelist)
      TYPE(ESMF_Time)   :: CURRTIME                                       !<-- The current forecast time (ESMF)
!
      INTEGER :: YY,MM,DD,H,M,S                                           !<-- Year, Month, Day, Hour, Minute, Second (integer)
!
      INTEGER :: I,INPES,JNPES                                          &
                ,RC,RC_ASYNC
!
      CHARACTER(ESMF_MAXSTR) :: filename                                  !<-- Restart/History label
!jwtest
      TYPE(ESMF_VM) :: myVM
      integer :: mypelocal,mypetcount,mympi
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  WHAT IS THE CURRENT FORECAST TIME?
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Get Current Time from ATM Clock"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_ClockGet(clock   =CLOCK_ATM                             &  !<-- The ATM component's ESMF Clock
                        ,currTime=CURRTIME                              &  !<-- The current forecast time (ESMF)
                        ,rc      =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Convert ESMF Time to Real Time"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_TimeGet (time=CURRTIME                                  &  !<-- The current forecast time (ESMF)
                        ,yy  =YY                                        &  !<-- The current year (integer)
                        ,mm  =MM                                        &  !<-- The current month (integer)
                        ,dd  =DD                                        &  !<-- The current day (integer)
                        ,h   =H                                         &  !<-- The current hour (integer)
                        ,m   =M                                         &  !<-- The current minute (integer)
                        ,s   =S                                         &  !<-- The current second (integer)
                        ,rc  =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      write(0,*)'in async,YY=',YY,'mm=',mm,'dd=',dd,'h=',h,'s=',s
!
!-----------------------------------------------------------------------
!***  THE EXPORT STATE OF THE DYNAMICS COMPONENT LIES WITHIN THE
!***  INTERNAL STATE OF THE ATM GRIDDED COMPONENT AND HOLDS THE
!***  IMPORT STATE OF THE WRITE COMPONENT.
!***  EXTRACT THAT WRITE COMPONENT'S IMPORT STATE SINCE WE ARE 
!***  ABOUT TO EXECUTE THE RUN STEP OF THE WRITE COMPONENT.
!-----------------------------------------------------------------------
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Extract Write Import State from Dyn Export State"
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      CALL ESMF_StateGet(state          =EXP_STATE    &  !<-- The Dyn component's export state
                        ,itemName       ="Write Import State"     &  !<-- Name of state to be extracted
                        ,nestedState    =IMP_STATE_WRITE  &  !<-- The extracted state
                        ,rc             =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      write(0,*)'in async,get imp_state_write,rc=',rc
!
!-----------------------------------------------------------------------
!***  ALL FORECAST TASKS PLUS THOSE WRITE TASKS IN THE APPROPRIATE
!***  WRITE GROUP EXECUTE THE RUN STEP OF A WRITE COMPONENT.
!-----------------------------------------------------------------------
!
      N_GROUP=WRITE_GROUP_READY_TO_GO                          !<-- The active write group
!
!jw
      call esmf_GridCompGet(gridcomp=WRT_COMPS(N_GROUP)     &
                     ,vm=myVM                        &
                     ,rc=rc)
!
     call esmf_vmget(VM=myVM                         &
                     ,localpet=myPElocal                  &
                     ,petcount=mypetcount             &
                     ,rc=rc)

    CALL ESMF_VMBarrier(myVM,rc=RC)    ! Insert barrier since fcst tasks are involved in each iteration of write groups
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
      MESSAGE_CHECK="WRITE_ASYNC: Execute Run Step of Write Components" 
!     CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
      DO I=1,NUM_PES_FCST+WRITE_TASKS_PER_GROUP
        IF(MYPE==PETLIST_WRITE(I,N_GROUP))THEN
          CALL ESMF_GridCompRun(WRT_COMPS(N_GROUP)          &  !<-- The write gridded component
                               ,importState=IMP_STATE_WRITE &  !<-- Its import state
                               ,exportState=EXP_STATE_WRITE &  !<-- Its export state
                               ,clock      =CLOCK_ATM                     &  !<-- The ATM Clock
                               ,phase      =ESMF_SINGLEPHASE              &
                               ,rc         =RC)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ERR_MSG(RC,MESSAGE_CHECK,RC_ASYNC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          IF(I==NUM_PES_FCST+1)THEN                                          !<-- The first write task tells us the history output time
            WRITE(0,101)YY,MM,DD,H,M,S
  101       FORMAT(' Wrote File at ',I4.4,'_',I2.2,'_',I2.2,'_',I2.2,':',I2.2,':',I2.2)
          ENDIF
!
        write(0,*)'mype after grid comp run'
        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  PREPARE TO USE THE NEXT WRITE GROUP AT THE NEXT OUTPUT TIME.
!***  RETURN TO THE 1ST GROUP IF WE HAVE CYCLED THROUGH ALL OF THEM.
!-----------------------------------------------------------------------
!
      IF(WRITE_GROUP_READY_TO_GO==WRITE_GROUPS)THEN
        WRITE_GROUP_READY_TO_GO=1
      ELSE
        WRITE_GROUP_READY_TO_GO=WRITE_GROUP_READY_TO_GO+1
      ENDIF
        write(0,*)'before end asyn'
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_ASYNC_GFS
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
      SUBROUTINE WRITE_NEMSIO_OPEN(WRT_INT_STATE             &
                                        ,NBDL                           &
                                        ,NEMSIOFILE                     &
                                        ,IYEAR_FCST                     &
                                        ,IMONTH_FCST                    &
                                        ,IDAY_FCST                      &
                                        ,IHOUR_FCST                     &
                                        ,IMINUTE_FCST                   &
                                        ,SECOND_FCST                    &
                                        ,NF_HOURS                       &
                                        ,NF_MINUTES                     &
                                        ,NF_SECONDS                     &
                                        ,DIM1,DIM2,NFRAME               &
                                        ,LEAD_WRITE_TASK)
!
!-----------------------------------------------------------------------
!***  WRITE OUT A NEMSIO BINARY RUN HISTORY FILE.
!-----------------------------------------------------------------------
!
      TYPE(WRITE_INTERNAL_STATE_GFS),INTENT(INOUT) :: WRT_INT_STATE             !<-- The Write component's internal state
!
      TYPE(NEMSIO_GFILE),INTENT(INOUT)         :: NEMSIOFILE                !<-- The nemsio file handler
!
      INTEGER,INTENT(IN)  :: IYEAR_FCST                                 &
                            ,IMONTH_FCST                                &
                            ,IDAY_FCST                                  &
                            ,IHOUR_FCST                                 &
                            ,IMINUTE_FCST                               &
                            ,NF_HOURS                                   &
                            ,NF_MINUTES                                 &
                            ,LEAD_WRITE_TASK

      INTEGER,INTENT(IN)  :: NBDL
      INTEGER,INTENT(OUT) :: DIM1,DIM2,NFRAME
!
      REAL,INTENT(IN)     :: NF_SECONDS                                 &
                            ,SECOND_FCST
!
!-----------------------------------------------------------------------
!***  LOCAL VARIABLES
!-----------------------------------------------------------------------
!
      INTEGER :: I,J,N,N1,N2,NPOSN_1,NPOSN_2,LENGTH,MAXLENGTH
!
      INTEGER :: FIELDSIZE,IM,JM,LM,IDATE(7),FCSTDATE(7)                &
                ,INDX_2D,INDX_2D2,INDX_2D3                              &
                ,IRET,IND1,IND2,IND3,IND4,CNT,INI1,INI2                 &  
                ,N2ISCALAR,N2IARY,N2RSCALAR,N2RARY,N2LSCALAR            &
                ,NMETA,TLMETA,VLEV,NSOIL
!
      INTEGER :: NFIELD,RC
!
      INTEGER,DIMENSION(:),POINTER :: ARYILEN                           &
                                     ,ARYRLEN                           &
                                     ,RECLEV                            &
                                     ,VARIVAL
!
      INTEGER,DIMENSION(:,:),POINTER :: ARYIVAL
!
      REAL(4),DIMENSION(:,:,:),POINTER :: VCOORD
!
!jw
      REAL(KIND=kind_io4),DIMENSION(:)  ,POINTER :: VARRVAL
      REAL(KIND=kind_io4),DIMENSION(:,:),POINTER :: ARYRVAL
      REAL(KIND=kind_io4),DIMENSION(:)  ,POINTER :: CPI,RI 
!
      LOGICAL                      :: GLOBAL
      LOGICAL                      :: HYBRID
      LOGICAL                      :: GEN_COORD_HYBRID
      LOGICAL,DIMENSION(:),POINTER :: VARLVAL
!
      CHARACTER(6)                       :: MODEL_LEVEL
      CHARACTER(16)                       :: VLEVTYP
!
      CHARACTER(8),DIMENSION(:) ,POINTER :: ARYINAME                    &
                                           ,ARYRNAME                    &
                                           ,RECNAME                     &
                                           ,VARINAME                    &
                                           ,VARRNAME                    &
                                           ,VARLNAME
!
      CHARACTER(16),DIMENSION(:),POINTER :: RECLEVTYP
!
      CHARACTER(ESMF_MAXSTR) :: NAME,FILENAME
!
      TYPE(ESMF_Logical) :: WORK_LOGICAL
!
      INTEGER :: NREC=0    
!
      INTEGER :: IDVC,IDVM,IDSL,IDRT,NVCOORD
!
      INTEGER :: NAK5,NBK5,NCK5,NSI,NIDVC                     &
                ,NIDVM,NIDSL,NIDRT,Nthermodyn_id              &
                ,Nsfcpress_id,Nvertcoord_id,NTRAC,NCLD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!initialization
!
      NAK5=0;NBK5=0;NCK5=0;NSI=0;NIDVC=0
      NIDVM=0;NIDSL=0;NIDRT=0;Nthermodyn_id=0
      Nsfcpress_id=0;Nvertcoord_id=0;NTRAC=0;NCLD=0
!
      IDVC=0;IDVM=0;IDSL=0;IDRT=0;NVCOORD=0
!
      FCSTDATE(1)=IYEAR_FCST
      FCSTDATE(2)=IMONTH_FCST
      FCSTDATE(3)=IDAY_FCST
      FCSTDATE(4)=IHOUR_FCST
      FCSTDATE(5)=IMINUTE_FCST
      FCSTDATE(6)=nint(SECOND_FCST*100.)
      FCSTDATE(7)=100
      write(0,*)'fcstdate=',fcstdate(1:6),'kount_I1d=',wrt_int_state%KOUNT_I1D(NBDL), &
       'KOUNT_R1D=',wrt_int_state%KOUNT_R1D(NBDL),'LOG=',wrt_int_state%KOUNT_LOG(NBDL),'NAK5=',NAK5
!
!-----------------------------------------------------------------------
!***  INTEGER SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
!-------------------------------------------------------------
!*** Find out the total number of int scalars and int arrays.
!-------------------------------------------------------------
!
      N2ISCALAR=0
      N2IARY=0
      MAXLENGTH=1
!
      DO N=1,wrt_int_state%KOUNT_I1D(NBDL)                                    !<-- Loop through all scalar/1D integer data
        LENGTH=wrt_int_state%LENGTH_DATA_I1D(N,NBDL)
!
        IF(LENGTH==1)THEN
          N2ISCALAR=N2ISCALAR+1
        ELSE
          N2IARY=N2IARY+1
          MAXLENGTH=MAX(LENGTH,MAXLENGTH)
        ENDIF
!
      ENDDO
!
      N2IARY=N2IARY+1
      MAXLENGTH=MAX(MAXLENGTH,7)
      ALLOCATE(VARINAME(N2ISCALAR),VARIVAL(N2ISCALAR))
      ALLOCATE(ARYINAME(N2IARY),ARYILEN(N2IARY),ARYIVAL(MAXLENGTH,N2IARY))
!
!---------------------------------------
!***  SET VALUE TO AVRIVAL and ARYIVAL.
!---------------------------------------
!
      N2=0                                                                 !<-- Word counter for full string of integer scalar/1D data
      N2ISCALAR=0
      N2IARY=0
      IDATE=0
      ARYIVAL=0
!
      DO N=1,wrt_int_state%KOUNT_I1D(NBDL)                                    !<-- Loop through all scalar/1D integer data
!
        NPOSN_1=(N-1)*NAME_MAXSTR+1
        NPOSN_2=N*NAME_MAXSTR
        NAME=wrt_int_state%NAMES_I1D_STRING(NBDL)(NPOSN_1:NPOSN_2)               !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_I1D(N,NBDL)                            !<-- The variable's length in words
         write(0,*)'in I1D ,name=',trim(name),'len=',length
!
        IF(LENGTH==1)THEN
          N2=N2+1
          N2ISCALAR=N2ISCALAR+1
          VARINAME(N2ISCALAR)=TRIM(NAME)
          VARIVAL(N2ISCALAR)=wrt_int_state%ALL_DATA_I1D(N2,NBDL)
          IF(VARINAME(N2ISCALAR)=='IHRST') then
            IDATE(4)=VARIVAL(N2ISCALAR)
!jws
          ELSEIF(VARINAME(N2ISCALAR)=='idrt') then
            IDRT=VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR)=='idsl') then
            IDSL=VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR)=='idvm') then
            IDVM=VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR)=='idvc') then
            IDVC=VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR)=='sfcpress_id') then
            Nsfcpress_id=VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR)=='vertcoord_id') then
            Nvertcoord_id=VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR)=='thermodyn_id') then
            Nthermodyn_id=VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR)=='ntrac') then
            NTRAC=VARIVAL(N2ISCALAR)
          ELSEIF(VARINAME(N2ISCALAR)=='ncld') then
            NCLD=VARIVAL(N2ISCALAR)

          ENDIF
        ELSE
          N2IARY=N2IARY+1
          ARYINAME(N2IARY)=TRIM(NAME)
          ARYILEN(N2IARY)=LENGTH
!            write(0,*)'in I1D array,aryiname=',aryiname(N2IARY),'len=',aryilen(N2IARY),  &
!              wrt_int_state%ALL_DATA_I1D(N2+1:N2+length)
!
          DO N1=1,LENGTH
            N2=N2+1
            ARYIVAL(N1,N2IARY)=wrt_int_state%ALL_DATA_I1D(N2,NBDL)              !<-- Extract the individual data from the data string
          ENDDO

          IF(ARYINAME(N2IARY)=='IDAT') THEN
            IDATE(1)=ARYIVAL(3,N2IARY)
            IDATE(2)=ARYIVAL(2,N2IARY)
            IDATE(3)=ARYIVAL(1,N2IARY)
          ELSEIF(trim(ARYINAME(N2IARY))=='IDATE') THEN
            IDATE(4)=ARYIVAL(1,N2IARY)
            IDATE(2)=ARYIVAL(2,N2IARY)
            IDATE(3)=ARYIVAL(3,N2IARY)
            IDATE(1)=ARYIVAL(4,N2IARY)
          ENDIF
          IDATE(7)=100.
!
            write(0,*)'in I1D array,aryival=',aryival(:,N2IARY)
        ENDIF
!
      ENDDO
!
!-----------------------------
!***  Add fcst_date into ARYI
!-----------------------------
!
      N2IARY=N2IARY+1
      ARYINAME(N2IARY)='FCSTDATE'
      ARYILEN(N2IARY)=7
      ARYIVAL(1:7,N2IARY)=FCSTDATE(1:7)
!
!*** prepare cpi,ri
      allocate(cpi(ntrac+1),ri(ntrac+1))
      cpi=0.
      ri=0.
!
!-----------------------------------------------------------------------
!***  REAL SCALAR/1D HISTORY VARIABLES
!-----------------------------------------------------------------------
!
!------------------------------------------------------------
!***  Find the total number of real scalars and real arrays.
!------------------------------------------------------------
!
      N2RSCALAR=0
      N2RARY=0
      MAXLENGTH=1
!
      DO N=1,wrt_int_state%KOUNT_R1D(NBDL)                                    !<-- Loop through all scalar/1D real data
        LENGTH=wrt_int_state%LENGTH_DATA_R1D(N,NBDL)                       !<-- The variable's length
        IF(LENGTH==1)THEN
           N2RSCALAR=N2RSCALAR+1
        ELSE
          N2RARY=N2RARY+1
          MAXLENGTH=MAX(LENGTH,MAXLENGTH)
        ENDIF
      ENDDO
      write(0,*)'NBDL=',NBDL,'N2RSCALAR=',N2RSCALAR,'N2RARY=',N2RARY
!
      IF(N2RSCALAR>0) &
      ALLOCATE(VARRNAME(N2RSCALAR),VARRVAL(N2RSCALAR))
      IF(N2RARY>0) &
      ALLOCATE(ARYRNAME(N2RARY),ARYRLEN(N2RARY),ARYRVAL(MAXLENGTH,N2RARY))
!
!------------------------------------------------------
!***  Set values for the real scalars and real arrays.
!------------------------------------------------------
      N2=0                                                                 !<-- Word counter for full string of real scalar/1D data
      N2RSCALAR=0
      N2RARY=0
!
      DO N=1,wrt_int_state%KOUNT_R1D(NBDL)                                    !<-- Loop through all scalar/1D real data
!
        NPOSN_1=(N-1)*NAME_MAXSTR+1
        NPOSN_2=N*NAME_MAXSTR
        NAME=wrt_int_state%NAMES_R1D_STRING(NBDL)(NPOSN_1:NPOSN_2)               !<-- The variable's name
        LENGTH=wrt_int_state%LENGTH_DATA_R1D(N,NBDL)                            !<-- The variable's length
        if(NBDL==2) write(0,*)'R1D, N=',N,' NAME=',trim(NAME)
!
        IF(LENGTH==1)THEN
          N2=N2+1
          N2RSCALAR=N2RSCALAR+1
          VARRNAME(N2RSCALAR)=TRIM(NAME)
          VARRVAL(N2RSCALAR)=wrt_int_state%ALL_DATA_R1D(N2,NBDL)
!
        ELSE
!
          N2RARY=N2RARY+1
          ARYRNAME(N2RARY)=TRIM(NAME)
          ARYRLEN(N2RARY)=LENGTH
!
          DO N1=1,LENGTH
            N2=N2+1
            ARYRVAL(N1,N2RARY)=wrt_int_state%ALL_DATA_R1D(N2,NBDL)              !<-- Extract the individual data from the data string
          ENDDO
!
          IF( TRIM(NAME)=='AK5') THEN
            NAK5=N2RARY
          ELSEIF (TRIM(NAME)=='BK5') THEN
            NBK5=N2RARY
          ELSEIF (TRIM(NAME)=='CK5') THEN
            NCK5=N2RARY
          ELSEIF (TRIM(NAME)=='SI') THEN
            NSI=N2RARY
          ELSEIF (TRIM(NAME)=='CPI') THEN
            cpi(1:length)=aryrval(1:length,n2rary)
          ELSEIF (TRIM(NAME)=='RI') THEN
            ri(1:length)=aryrval(1:length,n2rary)
          ENDIF
          write(0,*)'in nemsio, aryrname=',trim(ARYRNAME(N2RARY)),'length=', &
           ARYRLEN(N2RARY),'cpi=',cpi,'ri=',ri

        ENDIF
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  LOGICAL HISTORY VARIABLES
!-----------------------------------------------------------------------
!
      N2LSCALAR=wrt_int_state%KOUNT_LOG(NBDL)                                 !<-- Counter for full string of logical data
!
      IF(N2LSCALAR>0) then
        ALLOCATE(VARLNAME(N2LSCALAR),VARLVAL(N2LSCALAR))
        VARLVAL=.false.
      ENDIF
      N2LSCALAR=0
!
      DO N=1,wrt_int_state%KOUNT_LOG(NBDL)                                    !<-- Loop through all logical data
!
        NPOSN_1=(N-1)*NAME_MAXSTR+1
        NPOSN_2=N*NAME_MAXSTR
        NAME=wrt_int_state%NAMES_LOG_STRING(NBDL)(NPOSN_1:NPOSN_2)               !<-- The variable's name
!
        N2LSCALAR=N2LSCALAR+1
        WORK_LOGICAL=wrt_int_state%ALL_DATA_LOG(N2LSCALAR,NBDL)                 !<-- Extract the individual data from the data string
        VARLNAME(N2LSCALAR)=NAME
        if(WORK_LOGICAL==ESMF_TRUE)VARLVAL(N2LSCALAR)=.true.
        IF(TRIM(NAME)=='GLOBAL') GLOBAL=VARLVAL(N2LSCALAR)
        IF(TRIM(NAME)=='GEN_COORD_HYBRID') GEN_COORD_HYBRID=VARLVAL(N2LSCALAR)
        IF(TRIM(NAME)=='HYBRID') HYBRID=VARLVAL(N2LSCALAR)
        write(0,*)'in nemsio, log name=',trim(name),'GEN_COORD_HYBRID=',GEN_COORD_HYBRID
!
      ENDDO
!
!-----------------------------------------------------------------------
!***  NOW OPEN NEMSIO FILE
!-----------------------------------------------------------------------
!
      write(0,*)' OPEN_NEMSIO_FILE wrt_int_state%IO_NEMSIOFILE=',        &
          trim(wrt_int_state%FILENAME_BASE(NBDL)),'idate=',idate
!
      IF(wrt_int_state%IO_FILE(NBDL)=='DEFERRED')THEN
        N=LEN_TRIM(wrt_int_state%FILENAME_BASE(NBDL))
        WRITE(FILENAME,100)wrt_int_state%FILENAME_BASE(NBDL)(1:n)  &
                          ,wrt_int_state%NFHOUR,'_nemsio'
        write(0,*)'FILENAME=',trim(FILENAME),' n=',n
  100   FORMAT(A4,I3.3,A7)
      ELSE
        FILENAME=wrt_int_state%FILENAME_BASE(NBDL)//'_nemsio'
      ENDIF
!
!----------------------------------------------------
!***  Prepare variables needed by the nemsip header:
!----------------------------------------------------
!
!dimension
      IF(wrt_int_state%core=='NMMB' .and. GLOBAL) THEN
!for global im/jm for data field
        NFRAME=1
      ELSE
!for regional
        NFRAME=0
      ENDIF
      IM=wrt_int_state%im(1)
      JM=wrt_int_state%jm(1)
      DIM1=wrt_int_state%im(1)-2*NFRAME
      DIM2=wrt_int_state%jm(1)-2*NFRAME
!
      LM=wrt_int_state%LM(1)
!
!for nmmb whole domain
      FIELDSIZE=IM*JM
      NREC=wrt_int_state%kount_I2D(NBDL)+wrt_int_state%kount_R2D(NBDL)
      write(0,*)'before vcoord,nbdl=',nbdl,'Kount_i2d=',wrt_int_state%kount_I2D(NBDL) &
       ,'Kount_r2d=',wrt_int_state%kount_R2D(NBDL)
!
!vcoord
      ALLOCATE(VCOORD(LM+1,3,2))
      VCOORD=0.
      if(wrt_int_state%core=='gfs') then
        write(0,*)'in write routine,core=',wrt_int_state%core,'NAK5=',NAK5,   &
         'NBK5=',NBK5,'NCK5=',NCK5
        IF(gen_coord_hybrid) then
          idvc    = Nvertcoord_id
          idvm    = Nthermodyn_id*10 + Nsfcpress_id    ! 1: ln(ps) 2:ps   ! hmhj
          idsl    = 2    ! idsl=2 for middle of layer                   ! hmhj
          nvcoord = 3
          if(NAK5>0 .and. NBK5>0 .and. NCK5>0 ) then
        write(0,*)'in write routine,1,core=',wrt_int_state%core,'NAK5=',NAK5,   &
         'NBK5=',NBK5,'NCK5=',NCK5,'ak5=',ARYRVAL(1:LM+1,NAK5),' bk5=',       &
         ARYRVAL(1:LM+1,NBK5),'ck5=',ARYRVAL(1:LM+1,NCK5)
            VCOORD(1:LM+1,1,1)=ARYRVAL(1:LM+1,NAK5)*1000.
            VCOORD(1:LM+1,2,1)=ARYRVAL(1:LM+1,NBK5)
            VCOORD(1:LM+1,3,1)=ARYRVAL(1:LM+1,NCK5)*1000.
          endif
        ELSEIF(hybrid) then
          idvc    = 2                        ! for hybrid vertical coord.
          nvcoord = 2
          if(NAK5>0 .and. NBK5>0 ) then
            do i=1,LM+1
              vcoord(i,1,1) = ARYRVAL(LM+1+1-i,NAK5)*1000.
              vcoord(i,2,1) = ARYRVAL(LM+1+1-i,NBK5)
            enddo
            VCOORD(1:LM+1,3,1)=0
          endif
        ELSE
          idvc    = 1    ! for sigma vertical coord. (default)
          nvcoord = 1
          if(NSI>0 ) then
            vcoord(:,1,1) = ARYRVAL(:,NSI)
          endif
        ENDIF 
!
      endif 

      write(0,*)'after vcoord,count_I2d=',wrt_int_state%kount_I2D(NBDL),'nrec=',nrec
!
!-----------------------------------------------------------------------
!***  Cut the output I2D array.
!-----------------------------------------------------------------------
!
      ALLOCATE(RECNAME(NREC),RECLEVTYP(NREC),RECLEV(NREC))
      NREC=0
      INI1=0
      INI2=0
!
      DO NFIELD=1,wrt_int_state%KOUNT_I2D(NBDL)
!
        NREC=NREC+1
        NPOSN_1=(NFIELD-1)*NAME_MAXSTR+1
        NPOSN_2=NFIELD*NAME_MAXSTR
        NAME=wrt_int_state%NAMES_I2D_STRING(NBDL)(NPOSN_1:NPOSN_2)               !<-- The name of this 2D integer history quantity
        INDX_2D=index(NAME,"_2D")
!
        IF (INDX_2D > 0) THEN
          MODEL_LEVEL=NAME(INDX_2D-2:INDX_2D-1)
          RECLEV(NREC)=(ICHAR(MODEL_LEVEL(1:1))-48)*10+ICHAR(MODEL_LEVEL(2:2))-48
          INDX_2D2=INDEX(NAME,"_")
          RECNAME(NREC)=NAME(1:INDX_2D2-1)
          CALL LOWERCASE(RECNAME(NREC))
          if(INDX_2D-4>INDX_2D2)  then
             RECLEVTYP(NREC)=NAME(INDX_2D2+1:INDX_2D-4)
          else
             RECLEVTYP(NREC)='mid layer'
          endif
          IF (RECLEV(NREC)==LM+1) RECLEVTYP(NREC-LM:NREC)='layer'

          INI1=INI1+1
        ELSE
          RECLEV(NREC)=1
          INDX_2D2=INDEX(NAME,"_")
          INDX_2D3=INDEX(NAME,"_",BACK=.true.)
          if(INDX_2D2>0) then
            RECNAME(NREC)=trim(NAME(1:INDX_2D2-1))
          else
            RECNAME(NREC)=TRIM(NAME)
          endif
          if(INDX_2D3>0) then
            RECLEVTYP(NREC)=NAME(INDX_2D3+1:)
          else
            RECLEVTYP(NREC)='sfc'
          endif
          CALL LOWERCASE(RECNAME(NREC))

          RECLEV(NREC)=1
          INI2=INI2+1
        ENDIF
!
        IF (RECNAME(NREC)=='ISLTYP') RECNAME(NREC)='sltyp'
        IF (RECNAME(NREC)=='IVGTYP') RECNAME(NREC)='vgtyp'
        IF (RECNAME(NREC)=='NCFRCV') RECNAME(NREC)='cfrcv'
        IF (RECNAME(NREC)=='NCFRST') RECNAME(NREC)='cfrst'
        CALL LOWERCASE(RECNAME(NREC))
      ENDDO
!
      write(0,*)'after I2D,nrec=',nrec
!
!-----------------------------------------------------------------------
!*** Cut the output R2D array.
!-----------------------------------------------------------------------
!
      NSOIL=0
      IND4=0
      IND3=0
      IND2=0
      IND1=0
!
      DO NFIELD=1,wrt_int_state%KOUNT_R2D(NBDL)
!
        NREC=NREC+1
        NPOSN_1=(NFIELD-1)*NAME_MAXSTR+1
        NPOSN_2=NFIELD*NAME_MAXSTR
        NAME=wrt_int_state%NAMES_R2D_STRING(NBDL)(NPOSN_1:NPOSN_2)  !<-- The name of this 2D integer history quantity
        INDX_2D=INDEX(NAME,"_2D")
!
        IF (INDX_2D > 0) THEN
          MODEL_LEVEL=NAME(INDX_2D-2:INDX_2D-1)
          RECLEV(NREC)=(ICHAR(MODEL_LEVEL(1:1))-48)*10+ICHAR(MODEL_LEVEL(2:2))-48
          INDX_2D2=INDEX(NAME,"_")
          RECNAME(NREC)=NAME(1:INDX_2D2-1)
          CALL LOWERCASE(RECNAME(NREC))
          if(INDX_2D-4>INDX_2D2)  then
             RECLEVTYP(NREC)=NAME(INDX_2D2+1:INDX_2D-4)
          else
             RECLEVTYP(NREC)='mid layer'
          endif
          IF (RECLEV(NREC)==LM+1) RECLEVTYP(NREC-LM:NREC)='layer'
!
          IF (RECNAME(NREC)=='smc'.or.RECNAME(NREC)=='smc_soil layer') NSOIL=NSOIL+1
          IF (RECNAME(NREC)=='w') RECNAME(NREC)='vvel'
          IF (RECNAME(NREC)=='cw') RECNAME(NREC)='clwmr'
          IF (RECNAME(NREC)=='u') RECNAME(NREC)='ugrd'
          IF (RECNAME(NREC)=='v') RECNAME(NREC)='vgrd'
          IF (RECNAME(NREC)=='t') RECNAME(NREC)='tmp'
          IF (RECNAME(NREC)=='q') RECNAME(NREC)='spfh'
          IF (RECNAME(NREC)=='pint') THEN
          RECNAME(NREC)='pres'
          IND1=IND1+1
          ELSE IF (RECNAME(NREC)=='smc'.OR.RECNAME(NREC)=='sh2o'.or.RECNAME(NREC)=='stc' &
            .or. RECNAME(NREC)=='slc' ) THEN 
             RECLEVTYP(NREC)='soil layer' 
          IND2=IND2+1
          ELSE
          IND3=IND3+1
          ENDIF
        ELSE
          RECLEV(NREC)=1
          INDX_2D2=INDEX(NAME,"_")
          INDX_2D3=INDEX(NAME,"_",BACK=.true.)
          if(INDX_2D2>0) then
            RECNAME(NREC)=trim(NAME(1:INDX_2D2-1))
          else
            RECNAME(NREC)=TRIM(NAME)
          endif
          if(INDX_2D3>0) then
            RECLEVTYP(NREC)=NAME(INDX_2D3+1:)
          else
            RECLEVTYP(NREC)='sfc'
          endif
          CALL LOWERCASE(RECNAME(NREC))
!
          IF (RECNAME(NREC)=='pd') THEN
            RECNAME(NREC)='dpres'
            RECLEVTYP(NREC)='hybrid sig lev'
          ENDIF
!
          IF (RECNAME(NREC)=='pressfc') RECNAME(NREC)='pres'
          IF (RECNAME(NREC)=='sst') RECNAME(NREC)='tsea'
          IF (RECNAME(NREC)=='fis') RECNAME(NREC)='hgt'
          IF (RECNAME(NREC)=='ustar') RECNAME(NREC)='uustar'
          IF (RECNAME(NREC)=='z0') RECNAME(NREC)='zorl'
          IND4=IND4+1
        ENDIF
!
        write(0,*)'nrec=',nrec,'recname(nrec)=',recname(nrec),'reclev=',reclev(nrec), &
          'reclevtyp=',reclevtyp(nrec)
      ENDDO
!
!set meta data record
      NMETA=12
      IF (NAK5==0.and.NBK5==0.and.NCK5==0) NMETA=5
!
!-----------------------------------------------------------------------
!                      SET UP NEMSIO WRITE
!-----------------------------------------------------------------------
!
      CALL NEMSIO_INIT(IRET=IRET)
      write(0,*)'after nemsio_init, iret=',iret,'dim1=',dim1,'dim2=',dim2, &
       'nsoil=',nsoil,'ntrac=',ntrac,'ncldt=',ncld,'nrec=',nrec,'cpi=',cpi, &
       'ri=',ri,'nmetavarl=',n2lscalar,'nmetaaryr=',n2rary
!
!-----------------------------------------------------------------------
!***  OPEN NEMSIO FILE
!-----------------------------------------------------------------------
!
      CALL NEMSIO_OPEN(NEMSIOFILE,trim(FILENAME),'write',iret,           &
        modelname="GFS", gdatatype=wrt_int_state%io_form(NBDL),          &
        idate=IDATE,nfhour=NF_HOURS, &
        nfminute=NF_MINUTES,nfsecondn=nint(NF_SECONDS*100),              &
        nfsecondd=100,dimx=DIM1,dimy=DIM2,dimz=LM,nframe=NFRAME,         &
        nmeta=NMETA,                                                     &
        nsoil=NSOIL,ntrac=ntrac,nrec=nrec, ncldt=ncld,                   &
        vcoord=vcoord,cpi=cpi,ri=ri,idrt=idrt,                           &
        extrameta=.true.,nmetavari=N2ISCALAR,                            &
        nmetavarr=N2RSCALAR,nmetavarl=N2LSCALAR,nmetaaryi=N2IARY,        &
        nmetaaryr=N2RARY,variname=VARINAME,varival=VARIVAL,              &
        varrname=VARRNAME,varrval=VARRVAL,varlname=VARLNAME,             &
        varlval=VARLVAL,aryiname=ARYINAME,aryilen=ARYILEN,               &
        aryival=ARYIVAL,aryrname=ARYRNAME,aryrlen=ARYRLEN,               &
        aryrval=ARYRVAL,recname=RECNAME,reclevtyp=RECLEVTYP,reclev=RECLEV)
      write(0,*)'after nemsio_open, iret=',iret
!
!-----------------------------------------------------------------------
!***  CLEAN UP
!-----------------------------------------------------------------------
!
      DEALLOCATE(VCOORD,CPI,RI)
      if(associated(VARINAME))  DEALLOCATE(VARINAME,VARIVAL)
      if(associated(ARYINAME))  DEALLOCATE(ARYINAME,ARYILEN,ARYIVAL)
      if(associated(VARRNAME))  DEALLOCATE(VARRNAME,VARRVAL)
      if(associated(ARYRNAME))  DEALLOCATE(ARYRNAME,ARYRLEN,ARYRVAL)
      if(associated(varlname))  DEALLOCATE(VARLNAME,VARLVAL)
!
      write(0,*)'end of write_nemsio_open'
!-----------------------------------------------------------------------
!
      END SUBROUTINE WRITE_NEMSIO_OPEN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!
      elemental subroutine lowercase(word)
!
!-----------------------------------------------------------------------
!***  convert a word to lower case
!-----------------------------------------------------------------------
!
      character (len=*) , intent(inout) :: word
      integer :: i,ic,nlen
      nlen = len(word)
!
      do i=1,nlen
        ic = ichar(word(i:i))
        if (ic >= 65 .and. ic < 91) word(i:i) = char(ic+32)
      end do
!
!
!-----------------------------------------------------------------------
!
      end subroutine lowercase
!
!-----------------------------------------------------------------------
!
      SUBROUTINE CMONTH(IMON,CMON)
!
!-----------------------------------------------------------------------
!***  CONVERT MONTH
!-----------------------------------------------------------------------
!
      INTEGER,INTENT(IN) :: IMON
      CHARACTER(LEN=3)   :: CMON
!
!-----------------------------------------------------------------------
!
      SELECT CASE (IMON)
        CASE(1)
            CMON='Jan'
        CASE(2)
            CMON='Feb'
        CASE(3)
            CMON='Mar'
        CASE(4)
            CMON='Apr'
        CASE(5)
            CMON='May'
        CASE(6)
            CMON='Jun'
        CASE(7)
            CMON='Jul'
        CASE(8)
            CMON='Aug'
        CASE(9)
            CMON='Sep'
        CASE(10)
            CMON='Oct'
        CASE(11)
            CMON='Nov'
        CASE(12)
            CMON='Dec'
      END SELECT
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE CMONTH
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_WRITE_ROUTINES_GFS
!
!-----------------------------------------------------------------------
