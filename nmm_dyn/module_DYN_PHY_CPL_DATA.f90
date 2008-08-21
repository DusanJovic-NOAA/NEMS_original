!-----------------------------------------------------------------------
!
      MODULE MODULE_DYN_PHY_CPL_DATA
!
!-----------------------------------------------------------------------
!
!***  LIST THE ROOTS OF THE FIELD NAMES OF THE ARRAYS THAT WILL BE
!***  TRANSFERRED BETWEEN EXPORT AND IMPORT STATES OF THE DYNAMICS-
!***  PHYSICS COUPLER DURING THE INTEGRATION.
!***  THESE LISTS ARE THEN USED IN SIMPLE DO LOOPS WITHIN THE
!***  COUPLER TO DIRECT THE POINTERS AT THE DATA LOCATIONS.
!
!-----------------------------------------------------------------------
!
      USE ESMF_MOD
!
!-----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!-----------------------------------------------------------------------
!***  SET THE NUMBER OF 2- and 3-DIMENSIONAL FIELD NAMES
!***  AND HOW MANY TO BE TRANSFERRED BETWEEN COMPONENTS.
!-----------------------------------------------------------------------
!
!---------------------------------
!***  TOTAL NUMBER OF FIELD NAMES
!---------------------------------
!
      INTEGER,PARAMETER :: NDATA_2D=1       ! Number of 2D variables handled by the coupler
      INTEGER,PARAMETER :: NDATA_3D=5       ! Number of 3D variables handled by the coupler
!
!-----------------------------------------------------
!***  NUMBER OF FIELDS MOVED FROM DYNAMICS TO PHYSICS
!-----------------------------------------------------
!
      INTEGER,PARAMETER :: NDATA_2D_DYN_TO_PHY=NDATA_2D
      INTEGER,PARAMETER :: NDATA_3D_DYN_TO_PHY=NDATA_3D
!
!-----------------------------------------------------
!***  NUMBER OF FIELDS MOVED FROM PHYSICS TO DYNAMICS
!-----------------------------------------------------
!
      INTEGER,PARAMETER :: NDATA_2D_PHY_TO_DYN=NDATA_2D
      INTEGER,PARAMETER :: NDATA_3D_PHY_TO_DYN=NDATA_3D-1        ! We know a priori that one 3D Field (OMGALF) is not changed by Physics
!
!----------------------------------------------------------------
!***  THE NAMES OF THE FIELDS THAT WILL MOVE THROUGH THE
!***  DYN-PHY COUPLER.
!----------------------------------------------------------------
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(NDATA_2D) :: DATANAMES_2D        &
                                                     =(/'PD'/)
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(NDATA_3D) :: DATANAMES_3D        &
                                                     =(/'T     '        &
                                                       ,'U     '        &
                                                       ,'V     '        &
                                                       ,'Q2    '        &
                                                       ,'OMGALF'        &
                                                            /)
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYN_PHY_CPL_DATA
!
!-----------------------------------------------------------------------
