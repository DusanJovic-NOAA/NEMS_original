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
!***  Set the number of 2- AND 3-dimensional Field names
!***  and how many to be transferred between components.
!-----------------------------------------------------------------------
!
!----------------------------
!***  Total number of Fields
!----------------------------
!
      INTEGER,PARAMETER :: NDATA_2D=1       ! Number of 2D variables handled by the coupler
      INTEGER,PARAMETER :: NDATA_3D=7       ! Number of 3D variables handled by the coupler
!
!-------------------------------------------------
!***  Number of Fields moved from Dynamics export
!-------------------------------------------------
!
      INTEGER,PARAMETER :: NDATA_2D_FROM_DYN=NDATA_2D
      INTEGER,PARAMETER :: NDATA_3D_FROM_DYN=NDATA_3D
!
!------------------------------------------------
!***  Number of Fields moved from Physics export
!------------------------------------------------
!
      INTEGER,PARAMETER :: NDATA_2D_FROM_PHY=NDATA_2D
      INTEGER,PARAMETER :: NDATA_3D_FROM_PHY=NDATA_3D-1        ! We know a priori that one 3D Field (OMGALF) is not changed by Physics
!
!----------------------------------------------------------------
!***  The names of the Fields that will move through the
!***  Dyn-Phy Coupler.
!----------------------------------------------------------------
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(NDATA_2D) :: DATANAMES_2D        &
                                                     =(/'PD'/)
!
      CHARACTER(ESMF_MAXSTR),DIMENSION(NDATA_3D) :: DATANAMES_3D        &
                                                     =(/'T     '        &
                                                       ,'U     '        &
                                                       ,'V     '        &
                                                       ,'W     '        &
                                                       ,'Z     '        &
                                                       ,'Q2    '        &
                                                       ,'OMGALF'        &
                                                            /)
!
!-----------------------------------------------------------------------
!
      END MODULE MODULE_DYN_PHY_CPL_DATA
!
!-----------------------------------------------------------------------
