!-----------------------------------------------------------------------
      MODULE MODULE_WRITE_INTERNAL_STATE
!-----------------------------------------------------------------------
!***  THE INTERNAL STATE OF THE WRITE COMPONENT.
!-----------------------------------------------------------------------
!***
!***  HISTORY
!***
!       xx Feb 2007:  W. Yang - Originator
!       14 Jun 2007:  T. Black - Name revisions
!       14 Aug 2007:  T. Black - Some pointers changed to arrays
!       11 Sep 2007:  T. Black - Updates for quilting
!       15 Aug 2008:  J. Wang  - Add NEMSIO variables
!       16 Sep 2008:  J. Wang  - 3-D output arrays revert to 2-D
!       04 Sep 2009:  T. Black - Add the 1-D boundary restart arrays
!       22 Apr 2010:  T. Black - Add minutes and seconds to elapsed
!                                forecast time.
!
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
      USE MODULE_INCLUDE
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: WRITE_INTERNAL_STATE,WRITE_WRAP                         &
               ,MAX_DATA_I1D,MAX_DATA_I2D                               &
               ,MAX_DATA_R1D,MAX_DATA_R2D                               &
               ,MAX_DATA_LOG
!
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_I1D=50                      !<-- Max # of 1D integer arrays
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_I2D=50                      !<-- Max # of 2D integer arrays
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_R1D=50                      !<-- Max # of 1D real arrays
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_R2D=10000                   !<-- Max # of 2D real arrays and layers
                                                                           !    of all real 3D arrays combined
      INTEGER(kind=KINT),PARAMETER :: MAX_DATA_LOG=10
!
!-----------------------------------------------------------------------

      TYPE WRITE_INTERNAL_STATE

!--------------------------------
! PE INFORMATION AND TASK LAYOUT
!--------------------------------
!
      INTEGER(kind=KINT) :: MYPE
      INTEGER(kind=KINT) :: INPES,JNPES
      INTEGER(kind=KINT) :: IHALO,JHALO
      INTEGER(kind=KINT) :: NTASKS
      INTEGER(kind=KINT) :: WRITE_GROUPS,WRITE_TASKS_PER_GROUP
!
!-----------------------------
!***  FULL DOMAIN INFORMATION
!-----------------------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: IM
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: IDS
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: IDE
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JM
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JDS
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: JDE
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: LM
!
      TYPE(ESMF_Logical) :: GLOBAL
!
!--------------------
!***  SUBDOMAIN SIZE
!--------------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: LOCAL_ISTART       &
                                                    ,LOCAL_IEND         &
                                                    ,LOCAL_JSTART       &
                                                    ,LOCAL_JEND
!
!----------------------------------------------------
!***  IDs OF FCST TASKS THAT SEND TO EACH WRITE TASK
!----------------------------------------------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: ID_FTASK_RECV_STA  &
                                                    ,ID_FTASK_RECV_END
!
!----------------------------------------------------
!***  # OF WORDS SENT BY EACH FORECAST TASK 
!***  TO ITS DESIGNATED WRITE TASK.
!----------------------------------------------------
!
      INTEGER(kind=KINT) :: NUM_WORDS_SEND_I2D_HST                      &
                           ,NUM_WORDS_SEND_R2D_HST                      &
                           ,NUM_WORDS_SEND_I2D_RST                      & 
                           ,NUM_WORDS_SEND_R2D_RST
!
!----------------------------------------------------
!***  # OF WORDS RECEIVED BY EACH WRITE TASK FROM
!***  ALL OF ITS DESIGNATED FORECAST TASKS.
!----------------------------------------------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: NUM_WORDS_RECV_I2D_HST  &
                                                    ,NUM_WORDS_RECV_R2D_HST  &
                                                    ,NUM_WORDS_RECV_I2D_RST  &
                                                    ,NUM_WORDS_RECV_R2D_RST
!
!------------------------------
!***  HISTORY DATA INFORMATION
!------------------------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: KOUNT_I1D          &
                                                    ,KOUNT_I2D          &
                                                    ,KOUNT_R1D          &
                                                    ,KOUNT_R2D          &
                                                    ,KOUNT_LOG
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: RST_KOUNT_I1D      &
                                                    ,RST_KOUNT_I2D      &
                                                    ,RST_KOUNT_R1D      &
                                                    ,RST_KOUNT_R2D      &
                                                    ,RST_KOUNT_LOG
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: LENGTH_DATA_I1D    &
                                                    ,LENGTH_DATA_R1D    &
                                                    ,LENGTH_DATA_R2D
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: RST_LENGTH_DATA_I1D  &
                                                    ,RST_LENGTH_DATA_R1D  &
                                                    ,RST_LENGTH_DATA_R2D
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: LENGTH_SUM_I1D     &
                                                    ,LENGTH_SUM_R1D     &
                                                    ,LENGTH_SUM_R2D     &
                                                    ,LENGTH_SUM_LOG
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: RST_LENGTH_SUM_I1D  &
                                                    ,RST_LENGTH_SUM_R1D  &
                                                    ,RST_LENGTH_SUM_R2D  &
                                                    ,RST_LENGTH_SUM_LOG
!
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: NCOUNT_FIELDS
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: RST_NCOUNT_FIELDS
!
      INTEGER(kind=KINT),DIMENSION(:)  ,ALLOCATABLE :: ALL_DATA_I1D
      INTEGER(kind=KINT),DIMENSION(:)  ,ALLOCATABLE :: ALL_DATA_I2D
      INTEGER(kind=KINT),DIMENSION(:,:),ALLOCATABLE :: OUTPUT_ARRAY_I2D
!
      INTEGER(kind=KINT),DIMENSION(:)  ,ALLOCATABLE :: RST_ALL_DATA_I1D
      INTEGER(kind=KINT),DIMENSION(:)  ,ALLOCATABLE :: RST_ALL_DATA_I2D
      INTEGER(kind=KINT),DIMENSION(:,:),ALLOCATABLE :: RST_OUTPUT_ARRAY_I2D
!
      REAL(kind=KFPT),DIMENSION(:)  ,ALLOCATABLE :: ALL_DATA_R1D
      REAL(kind=KFPT),DIMENSION(:)  ,ALLOCATABLE :: ALL_DATA_R2D
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE :: OUTPUT_ARRAY_R2D
!
      REAL(kind=KFPT),DIMENSION(:)  ,ALLOCATABLE :: RST_ALL_DATA_R1D
      REAL(kind=KFPT),DIMENSION(:)  ,ALLOCATABLE :: RST_ALL_DATA_R2D
      REAL(kind=KFPT),DIMENSION(:,:),ALLOCATABLE :: RST_OUTPUT_ARRAY_R2D
!
!-----------------------------------------------------------------------
!***  BOUNDARY RESTART 
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT) :: LNSV                                           !<-- # of V bndry rows obtained from Dynamics
!
!-----------
!***  Local
!-----------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: NUM_WORDS_BC_SOUTH &  !<-- Word counts of 1-D boundary data strings
                                                    ,NUM_WORDS_BC_NORTH &  !    for each side of the domain.
                                                    ,NUM_WORDS_BC_WEST  &  !
                                                    ,NUM_WORDS_BC_EAST                !<--
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: RST_BC_DATA_SOUTH     &  !<-- 1-D strings of boundary data 
                                                 ,RST_BC_DATA_NORTH     &  !    for each side of the domain.
                                                 ,RST_BC_DATA_WEST      &  !
                                                 ,RST_BC_DATA_EAST         !<--
!
!-----------------
!***  Full-domain
!-----------------
!
      INTEGER(kind=KINT),DIMENSION(:),ALLOCATABLE :: NUM_WORDS_SEND_BC     !<-- Word count of full-domain 1-D boundary data string
!
      REAL(kind=KFPT),DIMENSION(:,:,:,:),ALLOCATABLE :: UBS,UBN,UBW,UBE &  !<-- 1-D string for U,V on each side of domain
                                                       ,VBS,VBN,VBW,VBE
!
      REAL(kind=KFPT),DIMENSION(:),ALLOCATABLE :: RST_ALL_BC_DATA          !<-- 1-D string of full-domain boundary data
!
!-----------------------------------------------------------------------
!*** STORAGE ARRAYS
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT),DIMENSION(:)    ,ALLOCATABLE :: BUFF_INT
      INTEGER(kind=KINT),DIMENSION(:,:,:),ALLOCATABLE :: WRITE_SUBSET_I
      REAL(kind=KFPT)   ,DIMENSION(:)    ,ALLOCATABLE :: BUFF_REAL
      REAL(kind=KFPT)   ,DIMENSION(:,:,:),ALLOCATABLE :: WRITE_SUBSET_R
!
      INTEGER(kind=KINT),DIMENSION(:)    ,ALLOCATABLE :: RST_BUFF_INT
      INTEGER(kind=KINT),DIMENSION(:,:,:),ALLOCATABLE :: RST_WRITE_SUBSET_I
      REAL(kind=KFPT)   ,DIMENSION(:)    ,ALLOCATABLE :: RST_BUFF_REAL
      REAL(kind=KFPT)   ,DIMENSION(:,:,:),ALLOCATABLE :: RST_WRITE_SUBSET_R
!
      TYPE(ESMF_Logical),DIMENSION(:),ALLOCATABLE :: ALL_DATA_LOG 
      TYPE(ESMF_Logical),DIMENSION(:),ALLOCATABLE :: RST_ALL_DATA_LOG
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(5000) :: FIELD_NAME
      CHARACTER(ESMF_MAXSTR),DIMENSION(5000) :: RST_FIELD_NAME
!
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I1D) :: NAMES_I1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I2D) :: NAMES_I2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R1D) :: NAMES_R1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R2D) :: NAMES_R2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_LOG) :: NAMES_LOG_STRING
!
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I1D) :: RST_NAMES_I1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_I2D) :: RST_NAMES_I2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R1D) :: RST_NAMES_R1D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_R2D) :: RST_NAMES_R2D_STRING
      CHARACTER(ESMF_MAXSTR*MAX_DATA_LOG) :: RST_NAMES_LOG_STRING
!
!-----------------------------------------------------------------------
!***  THE OUTPUT FILE
!-----------------------------------------------------------------------
!
      INTEGER(kind=KINT) :: IO_HST_UNIT,IO_RST_UNIT
      INTEGER(kind=KINT) :: IO_RECL
      INTEGER(kind=KINT) :: NFHOURS
      INTEGER(kind=KINT) :: NFMINUTES
!
      REAL(kind=KFPT)    :: NFSECONDS
!
      CHARACTER(ESMF_MAXSTR) :: IO_HST_FILE,IO_RST_FILE
      CHARACTER(ESMF_MAXSTR) :: HST_NAME_BASE,RST_NAME_BASE
      CHARACTER(ESMF_MAXSTR) :: IO_STATUS
      CHARACTER(ESMF_MAXSTR) :: IO_ACCESS
      CHARACTER(ESMF_MAXSTR) :: IO_FORM
      CHARACTER(ESMF_MAXSTR) :: IO_POSITION
      CHARACTER(ESMF_MAXSTR) :: IO_ACTION
      CHARACTER(ESMF_MAXSTR) :: IO_DELIM
      CHARACTER(ESMF_MAXSTR) :: IO_PAD
!
!-------------------------------------
!***  Times used in history filenames
!-------------------------------------
!
      TYPE(ESMF_Time)         :: IO_BASETIME
      TYPE(ESMF_TimeInterval) :: IO_CURRTIMEDIFF
!
!-----------------------------------------
!***  I/O direction flags (Read or Write)
!-----------------------------------------
!
      LOGICAL(kind=KLOG) :: WRITE_HST_FLAG,WRITE_RST_FLAG
      LOGICAL(kind=KLOG) :: WRITE_NEMSIOFLAG
      LOGICAL(kind=KLOG) :: WRITE_NEMSIOCTL
      LOGICAL(kind=KLOG) :: WRITE_DONEFILEFLAG
 
!-----------------------------------------------------------------------
!
      END TYPE WRITE_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!***  THIS STATE IS SUPPORTED BY C POINTERS BUT NOT F90 POINTERS
!***  THEREFORE WE NEED THIS WRAP.
!-----------------------------------------------------------
!
      TYPE WRITE_WRAP
        TYPE(WRITE_INTERNAL_STATE),POINTER :: WRITE_INT_STATE
      END TYPE WRITE_WRAP

!-----------------------------------------------------------
!
      END MODULE MODULE_WRITE_INTERNAL_STATE
!
!-----------------------------------------------------------
