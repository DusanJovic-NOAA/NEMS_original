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
!
!-----------------------------------------------------------------------
!
      USE ESMF_Mod
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
      INTEGER,PARAMETER :: MAX_DATA_I1D=50             !<-- Max # of 1D integer arrays
      INTEGER,PARAMETER :: MAX_DATA_I2D=50             !<-- Max # of 2D integer arrays
      INTEGER,PARAMETER :: MAX_DATA_R1D=50             !<-- Max # of 1D real arrays
      INTEGER,PARAMETER :: MAX_DATA_R2D=3000           !<-- Max # of 2D real arrays and layers of all real 3D arrays combined
      INTEGER,PARAMETER :: MAX_DATA_LOG=10
!
!-----------------------------------------------------------------------

      TYPE WRITE_INTERNAL_STATE

!--------------------------------
! PE INFORMATION AND TASK LAYOUT
!--------------------------------
!
      INTEGER :: MYPE
      INTEGER :: INPES,JNPES
      INTEGER :: IHALO,JHALO
      INTEGER :: NTASKS
      INTEGER :: WRITE_GROUPS,WRITE_TASKS_PER_GROUP
!
!-----------------------------
!***  FULL DOMAIN INFORMATION
!-----------------------------
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: IM
      INTEGER,DIMENSION(:),ALLOCATABLE :: JM
      INTEGER,DIMENSION(:),ALLOCATABLE :: LM
      TYPE(ESMF_Logical)           :: GLOBAL
!
!--------------------
!***  SUBDOMAIN SIZE
!--------------------
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: LOCAL_ISTART                  &
                                         ,LOCAL_IEND                    &
                                         ,LOCAL_JSTART                  &
                                         ,LOCAL_JEND
!
!----------------------------------------------------
!***  IDs OF FCST TASKS THAT SEND TO EACH WRITE TASK
!----------------------------------------------------
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: ID_FTASK_RECV_STA             &
                                         ,ID_FTASK_RECV_END
!
!------------------------------
!***  HISTORY DATA INFORMATION
!------------------------------
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: KOUNT_I1D                     &
                                         ,KOUNT_I2D                     &
                                         ,KOUNT_R1D                     &
                                         ,KOUNT_R2D                     &
                                         ,KOUNT_LOG
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: RST_KOUNT_I1D                     &
                                         ,RST_KOUNT_I2D                     &
                                         ,RST_KOUNT_R1D                     &
                                         ,RST_KOUNT_R2D                     &
                                         ,RST_KOUNT_LOG
!
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: LENGTH_DATA_I1D               &
                                         ,LENGTH_DATA_R1D               &
                                         ,LENGTH_DATA_R2D
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: RST_LENGTH_DATA_I1D               &
                                         ,RST_LENGTH_DATA_R1D               &
                                         ,RST_LENGTH_DATA_R2D
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: LENGTH_SUM_I1D                &
                                         ,LENGTH_SUM_R1D                &
                                         ,LENGTH_SUM_R2D                &
                                         ,LENGTH_SUM_LOG
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: RST_LENGTH_SUM_I1D                &
                                         ,RST_LENGTH_SUM_R1D                &
                                         ,RST_LENGTH_SUM_R2D                &
                                         ,RST_LENGTH_SUM_LOG
!
!
      INTEGER,DIMENSION(:),ALLOCATABLE :: NCOUNT_FIELDS
      INTEGER,DIMENSION(:),ALLOCATABLE :: RST_NCOUNT_FIELDS
!
      INTEGER,DIMENSION(:)  ,ALLOCATABLE :: ALL_DATA_I1D
      INTEGER,DIMENSION(:)  ,ALLOCATABLE :: ALL_DATA_I2D
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: OUTPUT_ARRAY_I2D
!
      INTEGER,DIMENSION(:)  ,ALLOCATABLE :: RST_ALL_DATA_I1D
      INTEGER,DIMENSION(:)  ,ALLOCATABLE :: RST_ALL_DATA_I2D
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: RST_OUTPUT_ARRAY_I2D
!
      REAL   ,DIMENSION(:)  ,ALLOCATABLE :: ALL_DATA_R1D
      REAL   ,DIMENSION(:)  ,ALLOCATABLE :: ALL_DATA_R2D
      REAL   ,DIMENSION(:,:),ALLOCATABLE :: OUTPUT_ARRAY_R2D
!
      REAL   ,DIMENSION(:)  ,ALLOCATABLE :: RST_ALL_DATA_R1D
      REAL   ,DIMENSION(:)  ,ALLOCATABLE :: RST_ALL_DATA_R2D
      REAL   ,DIMENSION(:,:),ALLOCATABLE :: RST_OUTPUT_ARRAY_R2D
!
!-----------------------------------------------------------------------
!*** STORAGE ARRAYS
!-----------------------------------------------------------------------
!
      INTEGER,DIMENSION(:)    ,ALLOCATABLE :: BUFF_INT
      INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: WRITE_SUBSET_I
      REAL   ,DIMENSION(:)    ,ALLOCATABLE :: BUFF_REAL
      REAL   ,DIMENSION(:,:,:),ALLOCATABLE :: WRITE_SUBSET_R
!
      INTEGER,DIMENSION(:)    ,ALLOCATABLE :: RST_BUFF_INT
      INTEGER,DIMENSION(:,:,:),ALLOCATABLE :: RST_WRITE_SUBSET_I
      REAL   ,DIMENSION(:)    ,ALLOCATABLE :: RST_BUFF_REAL
      REAL   ,DIMENSION(:,:,:),ALLOCATABLE :: RST_WRITE_SUBSET_R
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
      INTEGER                :: IO_HST_UNIT,IO_RST_UNIT
      INTEGER                :: IO_RECL
      INTEGER                :: NFHOUR
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
      LOGICAL :: WRITE_HST_FLAG,WRITE_RST_FLAG
      LOGICAL :: WRITE_NEMSIOFLAG
      LOGICAL :: WRITE_NEMSIOCTL
 
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
