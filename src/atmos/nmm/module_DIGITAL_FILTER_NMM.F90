#include "../../ESMFVersionDefine.h"

#if (ESMF_MAJOR_VERSION < 5 || ESMF_MINOR_VERSION < 2)
#undef ESMF_520r
#else
#define ESMF_520r
#endif

!-----------------------------------------------------------------------
!
      MODULE module_DIGITAL_FILTER_NMM
!
!-----------------------------------------------------------------------
!
! a generic digital filter for any model under ESMF 
!
!-----------------------------------------------------------------------
! March    2007	Hann-Ming Henry Juang
! February 2008 Weiyu Yang, updated to use the ESMF 3.1.0 library.
! February 2011 Weiyu Yang, Updated to use both the ESMF 4.0.0rp2 library,
!                           ESMF 5 library and the the ESMF 3.1.0rp2 library.
! May      2011 Weiyu Yang, Modified for using the ESMF 5.2.0r_beta_snapshot_07.
! September2011 Weiyu Yang, Modified for using the ESMF 5.2.0r library.
!----------------------------------------------------------------------------
!
      USE esmf_mod
      use module_include
      use module_exchange,only: halo_exch
      use module_dm_parallel, only: its,ite,jts,jte, lm

!      type(esmf_config),save :: cf_1                                !<-- The config object

      implicit none

! ---------
! dynamics
! ---------
      real, allocatable, save :: array_save_2d(:,:,:)
      real, allocatable, save :: array_save_3d(:,:,:,:)
      real, allocatable, save :: array_save_4d(:,:,:,:,:)
      real, allocatable, save :: dolph_wgts(:)
      real             , save :: totalsum
      character(20), allocatable, save :: name_save_2d(:)
      character(20), allocatable, save :: name_save_3d(:)
      character(20), allocatable, save :: name_save_4d(:)
      integer,                    save :: tot_rank_2d=0 
      integer,                    save :: tot_rank_3d=0
      integer,                    save :: tot_rank_4d=0
      integer,                    save :: kstep, nstep
      
! ---------
! physics
! ---------
      character(20), allocatable, save :: phy_name(:)
      real, allocatable, save :: array_save_2d_phys(:,:,:)
      real, allocatable, save :: array_save_3d_phys(:,:,:,:)
      character(20), allocatable, save :: name_save_2d_phys(:)
      character(20), allocatable, save :: name_save_3d_phys(:)
      integer,                    save :: tot_rank_2d_phys=0 
      integer,                    save :: tot_rank_3d_phys=0


      contains

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_dyn_init_nmm(filt_bundle                &
                                            ,ndfistep                   &
                                            ,dt_int,dt_num,dt_den       &
                                            ,num_water                  &
                                            ,num_tracers)
!-----------------------------------------------------------------------
!
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE
      integer, intent(in)          :: ndfistep
      integer, intent(in)          :: num_water,num_tracers
      integer, intent(in)          :: dt_int,dt_num,dt_den
!
!     type(esmf_field)             :: tmp_field 
      integer                      :: tmp_rank,dyn_items                &
                                     ,dfihr
      integer                      :: spec_max,rc,n,m,NUM_FIELDS
      character(20), allocatable   :: dyn_name(:)
      character(20)                :: state_name
      real                         :: taus, dt
      TYPE(ESMF_Field) :: tmpfield
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      nstep = ndfistep 
      kstep = - nstep -1

      if (allocated(dolph_wgts)) deallocate(dolph_wgts)
      allocate(dolph_wgts(-nstep:nstep))

      dt=float(dt_int)+float(dt_num)/float(dt_den)

!     hardwiring cutoff frequency based on length of filtering window

      taus=float(2*ndfistep)*dt
!
      call dolph(dt,taus,nstep,dolph_wgts)
!

!     Retrieve the dynamical fields to be filtered from the bundle

        CALL ESMF_FieldBundleGet(bundle    =FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,fieldCount=NUM_FIELDS          &  !<-- # of Fields in the Bundle
                                ,rc        =RC)

      tot_rank_2d=0
      tot_rank_3d=0
      tot_rank_4d=0

      if (.not. allocated(dyn_name))                                    &
      allocate(dyn_name(NUM_FIELDS))

      if (.not. allocated(name_save_2d))                                &
      allocate(name_save_2d(NUM_FIELDS))

      if (.not. allocated(name_save_3d))                                &
      allocate(name_save_3d(NUM_FIELDS))

      if (.not. allocated(name_save_4d))                                &
      allocate(name_save_4d(NUM_FIELDS))

      DO N=1,NUM_FIELDS

        CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,N                               &
                                ,tmpfield                        &
                                ,rc)

        CALL ESMF_FieldGet(field=tmpfield, name=dyn_name(N), dimCount=tmp_rank, rc=rc)
!
        IF (tmp_rank == 2) THEN
          tot_rank_2d=tot_rank_2d+1
          name_save_2d(tot_rank_2d)=dyn_name(N)
        ENDIF
!
        IF (tmp_rank == 3) THEN
          tot_rank_3d=tot_rank_3d+1
          name_save_3d(tot_rank_3d)=dyn_name(N)
        ENDIF
!
        IF (tmp_rank ==4 ) THEN
          tot_rank_4d=tot_rank_4d+1  
          name_save_4d(tot_rank_4d)=dyn_name(N)
        ENDIF
!
      ENDDO 
!
!      SPEC_MAX=MAX(NUM_TRACERS,NUM_WATER)
!     
      IF (tot_rank_2d > 0 .and. .not. allocated(array_save_2d)) THEN
      	allocate(array_save_2d(ITS:ITE,JTS:JTE,tot_rank_2d))
      ENDIF
!
      IF (tot_rank_3d > 0 .and. .not. allocated(array_save_3d)) THEN
      	allocate(array_save_3d(ITS:ITE,JTS:JTE,LM,tot_rank_3d))
      ENDIF
!
      IF (tot_rank_4d > 0 .and. .not. allocated(array_save_4d)) THEN
      	allocate(array_save_4d(ITS:ITE,JTS:JTE,LM,SPEC_MAX,tot_rank_4d))
      ENDIF
!
      IF (allocated(array_save_2d)) THEN
        array_save_2d=0.
      ENDIF

      IF (allocated(array_save_3d)) THEN
        array_save_3d=0.
      ENDIF

      IF (allocated(array_save_4d)) THEN
        array_save_4d=0.
      ENDIF
     
      deallocate(dyn_name)
      totalsum=0.
!-----------------------------------------------------------------------

      end subroutine digital_filter_dyn_init_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_dyn_sum_nmm(filt_bundle               &
                                           ,mean_on                     &
                                           ,num_water                   &
                                           ,num_tracers)
!-----------------------------------------------------------------------

      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE
      integer(kind=kint),intent(in) :: mean_on,num_water,num_tracers
!
      integer(kind=kint) :: i,ii,j,jj,l,n,num_spec,p,rc,rc_upd
      real(kind=kfpt) :: digfil,prod,sx,wx
      real(kind=kfpt),dimension(:,:)    ,pointer :: hold_2d
      real(kind=kfpt),dimension(:,:,:)  ,pointer :: hold_3d
      real(kind=kfpt),dimension(:,:,:,:),pointer :: hold_4d
      logical :: dolph
!
      character(20) :: field_name
!
      type(ESMF_Field) :: hold_field
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      rc    =esmf_success
      rc_upd=esmf_success
!
      kstep = kstep + 1

!-----------------------------------------------------------------------
!     Future task: make this dolph switch logical a configure file item
!-----------------------------------------------------------------------

      dolph=.true.

      IF (dolph) THEN
        digfil=dolph_wgts(kstep)

      ELSE 
        sx = acos(-1.)*kstep/nstep
        wx = acos(-1.)*kstep/(nstep+1)
        if( kstep/=0)then
          digfil = sin(wx)/wx*sin(sx)/sx
        else
          digfil=1.
        endif 

        if(mean_on>0)then
          digfil=1.
        endif

      ENDIF
!
      totalsum = totalsum + digfil
!
      if(tot_rank_2d>0) then
!
        do n=1,tot_rank_2d
          field_name=name_save_2d(N)
          nullify(hold_2d) 

          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)

          call ESMF_FieldGet(field     =HOLD_FIELD                      &  !<-- Field that holds the data pointer
                            ,localDe   =0                               &
                            ,farrayPtr =HOLD_2D                         &  !<-- Put the pointer here
                            ,rc        =RC)

          do j=jts,jte
          do i=its,ite
            array_save_2d(i,j,n)=array_save_2d(i,j,n)+digfil*hold_2d(i,j)
          enddo
          enddo
        enddo
      endif
!
      if(tot_rank_3d>0)then
        do n=1,tot_rank_3d
          field_name=name_save_3d(N)
          nullify(hold_3d)

          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)

          call ESMF_FieldGet(field     =HOLD_FIELD                      &  !<-- Field that holds the data pointer
                            ,localDe   =0                               &
                            ,farrayPtr =HOLD_3D                         &  !<-- Put the pointer here
                            ,rc        =RC)

          do l=1,lm  
            do j=jts,jte
            do i=its,ite
              array_save_3d(i,j,l,n)=array_save_3d(i,j,l,n)+digfil*hold_3d(i,j,l)
            enddo
            enddo
         enddo
        enddo
!
      endif
!
     if(tot_rank_4d>0)then
        do n=1,tot_rank_4d
          field_name=name_save_4d(N)
          nullify(hold_4d)
          if (field_name == 'TRACERS') then
            num_spec=num_tracers
          else if (field_name == 'WATER') then
            num_spec=num_water
          endif

          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)

          call ESMF_FieldGet(field     =HOLD_FIELD                      &  !<-- Field that holds the data pointer
                            ,localDe   =0                               &
                            ,farrayPtr =HOLD_4D                         &  !<-- Put the pointer here
                            ,rc        =RC)

          do p=1,num_spec
            do l=1,lm
              do j=jts,jte
              do i=its,ite
                array_save_4d(i,j,l,p,n)=array_save_4d(i,j,l,p,n)+digfil*hold_4d(i,j,l,p)
              enddo
              enddo
            enddo
          enddo
        enddo

      endif


      if(rc_upd==esmf_success)then
!       write(0,*)'DYNAMICS UPDATE SUCCEEDED'
      else
        write(0,*)'DYNAMICS UPDATE FAILED RC_UPD=',rc_upd
      endif
!-----------------------------------------------------------------------

      end subroutine digital_filter_dyn_sum_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_dyn_average_nmm(filt_bundle               &
                                               ,num_water               &
                                               ,num_tracers)
!-----------------------------------------------------------------------
!
      USE MODULE_DM_PARALLEL
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE
      INTEGER, intent(in)          :: NUM_WATER,NUM_TRACERS
!
      INTEGER(KIND=KINT) :: I,II,J,JJ,L,N,P,RC,RC_UPD
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
      REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: HOLD_4D
!
      CHARACTER(20)    :: FIELD_NAME
!
      TYPE(ESMF_Field) :: HOLD_FIELD

!     TYPE(ESMF_Field)                :: tmp_field
!
      CHARACTER(ESMF_Maxstr)          :: name
      real, dimension(:,:), pointer   :: tmp_ptr
      real                            :: totalsumi
      integer                         :: NUM_SPEC, NUM_FIELDS

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC    =ESMF_SUCCESS
      RC_UPD=ESMF_SUCCESS
!
      totalsumi = 1.0 / totalsum   

      IF (tot_rank_2d > 0) THEN
        DO N=1,tot_rank_2d
          DO J=JTS,JTE
          DO I=ITS,ITE
            array_save_2d(I,J,N)=totalsumi*array_save_2d(I,J,N)
          ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      IF (tot_rank_3d > 0) THEN
        DO N=1,tot_rank_3d
          DO L=1,LM
            DO J=JTS,JTE
            DO I=ITS,ITE
              array_save_3d(I,J,L,N)=totalsumi*array_save_3d(I,J,L,N)
            ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      IF (tot_rank_4d > 0) THEN
        DO N=1,tot_rank_4d
          FIELD_NAME=name_save_4d(N)
          IF (FIELD_NAME == 'TRACERS') THEN
            NUM_SPEC=NUM_TRACERS
          ELSE IF (FIELD_NAME == 'WATER') THEN
            NUM_SPEC=NUM_WATER
          ENDIF
          DO P=1,NUM_SPEC
            DO L=1,LM
              DO J=JTS,JTE
              DO I=ITS,ITE
                array_save_4d(I,J,L,P,N)=totalsumi*array_save_4d(I,J,L,P,N)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
!
      IF (tot_rank_2d > 0) THEN
        DO N=1,tot_rank_2d
          FIELD_NAME=name_save_2d(N)
          NULLIFY(HOLD_2D)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

          call ESMF_FieldGet(field     =HOLD_FIELD                      &  !<-- Field that holds the data pointer
                            ,localDe   =0                               &
                            ,farrayPtr =HOLD_2D                         &  !<-- Put the pointer here
                            ,rc        =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DO J=JTS,JTE
          DO I=ITS,ITE
            HOLD_2D(I,J)=array_save_2d(I,J,N)
          ENDDO
          ENDDO

          CALL HALO_EXCH(hold_2d,1,2,2)
        ENDDO
      ENDIF
!
      IF (tot_rank_3d > 0) THEN
        DO N=1,tot_rank_3d
          FIELD_NAME=name_save_3d(N)
          NULLIFY(HOLD_3D)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

          call ESMF_FieldGet(field     =HOLD_FIELD                      &  !<-- Field that holds the data pointer
                            ,localDe   =0                               &
                            ,farrayPtr =HOLD_3D                         &  !<-- Put the pointer here
                            ,rc        =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DO L=1,LM
            DO J=JTS,JTE
            DO I=ITS,ITE
              HOLD_3D(I,J,L)=array_save_3d(I,J,L,N)
            ENDDO
            ENDDO
          ENDDO

          CALL HALO_EXCH(hold_3d,LM,2,2)
        ENDDO
      ENDIF
!
      IF (tot_rank_4d > 0) THEN 
        DO N=1,tot_rank_4d 
          FIELD_NAME=name_save_4d(N)
          IF (FIELD_NAME == 'TRACERS') THEN
            NUM_SPEC=NUM_TRACERS
         ELSE IF (FIELD_NAME == 'WATER') THEN
            NUM_SPEC=NUM_WATER
          ENDIF
          NULLIFY(HOLD_4D)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         MESSAGE_CHECK="Dyn Update: Extract Temperature Pointer from Field"
!         CALL ESMF_LogWrite(MESSAGE_CHECK,ESMF_LOGMSG_INFO,rc=RC)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

          call ESMF_FieldGet(field     =HOLD_FIELD                      &  !<-- Field that holds the data pointer
                            ,localDe   =0                               &
                            ,farrayPtr =HOLD_4D                         &  !<-- Put the pointer here
                            ,rc        =RC)

! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!         CALL ERR_MSG(RC,MESSAGE_CHECK,RC_UPD)
! ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
!
          DO P=1,NUM_SPEC
            DO L=1,LM
              DO J=JTS,JTE
              DO I=ITS,ITE
                HOLD_4D(I,J,L,P)=array_save_4d(I,J,L,P,N)
              ENDDO
              ENDDO
            ENDDO
          ENDDO
          CALL ESMF_FieldBundleAdd(bundle=FILT_BUNDLE            &  !<-- The Filt Bundle for Filtered variables
                                  ,field =HOLD_FIELD                    &  !<-- Add this Field to the Bundle
                                  ,rc    =RC )
        ENDDO
      ENDIF


      IF (tot_rank_2d > 0) THEN 
        deallocate(name_save_2d)
        deallocate(array_save_2d)
      ENDIF

      IF (tot_rank_3d > 0) THEN 
        deallocate(name_save_3d)
        deallocate(array_save_3d)
      ENDIF

      IF (tot_rank_4d > 0) THEN 
        deallocate(name_save_4d)
        deallocate(array_save_4d)
      ENDIF

      tot_rank_2d=0
      tot_rank_3d=0
      tot_rank_4d=0
      kstep=0
      nstep=0

!-----------------------------------------------------------------------
      end subroutine digital_filter_dyn_average_nmm
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_phy_init_nmm(filt_bundle)
!-----------------------------------------------------------------------
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE
!
      TYPE(ESMF_Field) :: tmpfield
      integer :: rc, NUM_FIELDS, n, tmp_rank
      character(len=20) :: field_name
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

!     Retrieve the dynamical fields to be filtered from the bundle

        CALL ESMF_FieldBundleGet(bundle    =FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,fieldCount=NUM_FIELDS          &  !<-- # of Fields in the Bundle
                                ,rc        =RC)
!
      tot_rank_2d_phys=0
      tot_rank_3d_phys=0

      if (.not. allocated(phy_name))                                    &
      allocate(phy_name(NUM_FIELDS))

      if (.not. allocated(name_save_2d_phys))                                &
      allocate(name_save_2d_phys(NUM_FIELDS))

      if (.not. allocated(name_save_3d_phys))                                &
      allocate(name_save_3d_phys(NUM_FIELDS))

      DO N=1,NUM_FIELDS

        CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,N                               &
                                ,tmpfield                        &
                                ,rc)

        CALL ESMF_FieldGet(field=tmpfield, name=phy_name(N), dimCount=tmp_rank, rc=rc)

        IF (tmp_rank == 2) THEN
          tot_rank_2d_phys=tot_rank_2d_phys+1
          name_save_2d_phys(tot_rank_2d_phys)=phy_name(N)
        ENDIF
!
        IF (tmp_rank == 3) THEN
          tot_rank_3d_phys=tot_rank_3d_phys+1
          name_save_3d_phys(tot_rank_3d_phys)=phy_name(N)
        ENDIF
!
      ENDDO 
!
      IF (tot_rank_2d_phys > 0 .and. .not. allocated(array_save_2d_phys)) THEN
      	allocate(array_save_2d_phys(ITS:ITE,JTS:JTE,tot_rank_2d_phys))
      ENDIF
!
      IF (tot_rank_3d_phys > 0 .and. .not. allocated(array_save_3d_phys)) THEN
      	allocate(array_save_3d_phys(ITS:ITE,JTS:JTE,LM,tot_rank_3d_phys))
      ENDIF
!
      deallocate(phy_name)
      array_save_2d_phys=0.
      array_save_3d_phys=0.

!-----------------------------------------------------------------------
!
      end subroutine digital_filter_phy_init_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_phy_save_nmm(filt_bundle)
!-----------------------------------------------------------------------
!      use module_dm_parallel, only: its,ite,jts,jte
!
      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE
!
      TYPE(ESMF_Field)             :: hold_field
      integer                      :: n, rc, i,j,l, LDIM
      character(len=20)            :: field_name
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
!      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,ALLOCATABLE :: TMP_T, TMP_PINT, TMP_Q


!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!

      if(tot_rank_2d_phys>0) then
!
        do n=1,tot_rank_2d_phys
          field_name=name_save_2d_phys(N)
          nullify(hold_2d) 

          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)

          call ESMF_FieldGet(field     =HOLD_FIELD                      &  !<-- Field that holds the data pointer
                            ,localDe   =0                               &
                            ,farrayPtr =HOLD_2D                         &  !<-- Put the pointer here
                            ,rc        =RC)

          do j=jts,jte
          do i=its,ite
            array_save_2d_phys(i,j,n)=hold_2d(i,j)
          enddo
          enddo

        enddo
      endif


      if(tot_rank_3d_phys>0)then
        do n=1,tot_rank_3d_phys
          field_name=name_save_3d_phys(N)
          nullify(hold_3d)

          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)

          call ESMF_FieldGet(field     =HOLD_FIELD                      &  !<-- Field that holds the data pointer
                            ,localDe   =0                               &
                            ,farrayPtr =HOLD_3D                         &  !<-- Put the pointer here
                            ,rc        =RC)

	LDIM=size(HOLD_3D,dim=3)

          do l=1, LDIM
            do j=jts,jte
            do i=its,ite
              array_save_3d_phys(i,j,l,n)=hold_3d(i,j,l)
            enddo
            enddo
         enddo
        enddo
!
      endif
!-----------------------------------------------------------------------
!
      end subroutine digital_filter_phy_save_nmm

!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
      subroutine digital_filter_phy_restore_nmm( FILT_BUNDLE )
!-----------------------------------------------------------------------
!

      TYPE(ESMF_FieldBundle),INTENT(INOUT) :: FILT_BUNDLE

      INTEGER(KIND=KINT) :: I,II,J,JJ,L,N,P,RC,RC_UPD,NUM_FIELDS,LDIM
      REAL(KIND=KFPT),DIMENSION(:,:)    ,POINTER :: HOLD_2D
      REAL(KIND=KFPT),DIMENSION(:,:,:)  ,POINTER :: HOLD_3D
!     REAL(KIND=KFPT),DIMENSION(:,:,:,:),POINTER :: HOLD_4D

!
      type(ESMF_Field) :: hold_field

      character(len=20):: vars_list(150)
      CHARACTER(20)    :: FIELD_NAME
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      IF (tot_rank_2d_phys > 0) THEN
        DO N=1,tot_rank_2d_phys
          FIELD_NAME=name_save_2d_phys(N)
          NULLIFY(HOLD_2D)

          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)

          call ESMF_FieldGet(field     =HOLD_FIELD                      &  !<-- Field that holds the data pointer
                            ,localDe   =0                               &
                            ,farrayPtr =HOLD_2D                         &  !<-- Put the pointer here
                            ,rc        =RC)

          DO J=JTS,JTE
          DO I=ITS,ITE
            HOLD_2D(I,J)=array_save_2d_phys(I,J,N)
          ENDDO
          ENDDO

          CALL HALO_EXCH(hold_2d,1,2,2)

        ENDDO
      ENDIF

      IF (tot_rank_3d_phys > 0) THEN
        DO N=1,tot_rank_3d_phys
          FIELD_NAME=name_save_3d_phys(N)

          NULLIFY(HOLD_3D)


          CALL ESMF_FieldBundleGet(FILT_BUNDLE         &  !<-- The ESMF Bundle of arrays to be filtered
                                ,name=field_name       &
                                ,field=HOLD_FIELD                        &
                                ,rc=rc)

          call ESMF_FieldGet(field           =HOLD_FIELD                &  !<-- Field that holds the data pointer
                            ,localDe         =0                         &
                            ,farrayPtr       =HOLD_3D                   &  !<-- Put the pointer here
                            ,rc              =RC)
!
         LDIM=size(HOLD_3D,dim=3)
!
         DO L=1,LDIM
          DO J=JTS,JTE
          DO I=ITS,ITE
            HOLD_3D(I,J,L)=array_save_3d_phys(I,J,L,N)
          ENDDO
          ENDDO
        ENDDO

          CALL HALO_EXCH(hold_3d,LDIM,2,2)

        ENDDO
      ENDIF

!-----------------------------------------------------------------------

      end subroutine digital_filter_phy_restore_nmm

!-----------------------------------------------------------------------

      end module module_digital_filter_nmm
!
!-----------------------------------------------------------------------

   SUBROUTINE dolph(deltat, taus, m, window)

!     calculation of dolph-chebyshev window or, for short,
!     dolph window, using the expression in the reference:
!
!     antoniou, andreas, 1993: digital filters: analysis,
!     design and applications. mcgraw-hill, inc., 689pp.
!
!     the dolph window is optimal in the following sense:
!     for a given main-lobe width, the stop-band attenuation
!     is minimal; for a given stop-band level, the main-lobe
!     width is minimal.

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)                  ::  m
      REAL, DIMENSION(0:2*M), INTENT(OUT)    ::  window
      REAL, INTENT(IN)                     :: deltat, taus

      ! local data
      integer, PARAMETER        :: NMAX = 5000
      REAL, dimension(0:NMAX)   :: t, w, time
      real, dimension(0:2*nmax) :: w2
      INTEGER                   :: NPRPE=0        ! no of pe
      CHARACTER*80              :: MES

      real    :: pi, thetas, x0, term1, term2, rr, r,db, sum, arg, sumw
      integer :: n, nm1, i, nt

      PI = 4*ATAN(1.D0)

!      print *, 'in dfcoef, deltat = ', deltat, 'taus=',taus

      N = 2*M+1
      NM1 = N-1

      THETAS = 2*PI*ABS(DELTAT/TAUS)
      X0 = 1/COS(THETAS/2)
      TERM1 = (X0 + SQRT(X0**2-1))**(FLOAT(N-1))
      TERM2 = (X0 - SQRT(X0**2-1))**(FLOAT(N-1))
      RR = 0.5*(TERM1+TERM2)
      R = 1/RR
      DB = 20*LOG10(R)


!      WRITE(0,'(1X,''DOLPH: M,N='',2I8)')M,N
!      WRITE(0,'(1X,''DOLPH: THETAS (STOP-BAND EDGE)='',F10.3)')THETAS
!      WRITE(0,'(1X,''DOLPH: R,DB='',2F10.3)')R, DB

      DO NT=0,M
         SUM = 1
         DO I=1,M
            ARG = X0*COS(I*PI/N)
            CALL CHEBY(T,NM1,ARG)
            TERM1 = T(NM1)
            TERM2 = COS(2*NT*PI*I/N)
            SUM = SUM + R*2*TERM1*TERM2
         ENDDO
         W(NT) = SUM/N
         TIME(NT) = NT
      ENDDO
!     fill in the negative-time values by symmetry.
      DO NT=0,M
         W2(M+NT) = W(NT)
         W2(M-NT) = W(NT)
      ENDDO

!     fill up the array for return
      SUMW = 0.
      DO NT=0,2*M
         SUMW = SUMW + W2(NT)
      ENDDO
!      WRITE(0,'(1X,''DOLPH: SUM OF WEIGHTS W2='',F10.4)')SUMW

      DO NT=0,2*M
         WINDOW(NT) = W2(NT)
      ENDDO

      RETURN

   END SUBROUTINE dolph


   SUBROUTINE cheby(t, n, x)

!     calculate all chebyshev polynomials up to order n
!     for the argument value x.

!     reference: numerical recipes, page 184, recurrence
!         t_n(x) = 2xt_{n-1}(x) - t_{n-2}(x) ,  n>=2.

      IMPLICIT NONE

      ! Arguments
      INTEGER, INTENT(IN)  :: n
      REAL, INTENT(IN)     :: x
      REAL, DIMENSION(0:N) :: t

      integer  :: nn

      T(0) = 1
      T(1) = X
      IF(N.LT.2) RETURN
      DO NN=2,N
         T(NN) = 2*X*T(NN-1) - T(NN-2)
      ENDDO

      RETURN

   END SUBROUTINE cheby
