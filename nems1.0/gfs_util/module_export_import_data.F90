!-----------------------------------------------------------------------
!
      module module_export_import_data
!
!-----------------------------------------------------------------------
!
!***  list the roots of the field names of the arrays that will be
!***  transferred between export and import states during the
!***  integration.  these lists can then be used in simple do loops
!***  to redirect the data's pointers.
!
!-----------------------------------------------------------------------
!
      use esmf_mod
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!***  set the number of 2- and 3-dimensional field names
!***  and how many to be transferred between components.
!-----------------------------------------------------------------------
!
!---------------------------------
!***  total number of field names
!---------------------------------
!
      integer,parameter :: ndata_1d=1
      integer,parameter :: ndata_2d=2
      integer,parameter :: ndata_3d=9
!
!-----------------------------------------------------
!***  number of fields moved from dynamics 
!-----------------------------------------------------
!
      integer,parameter :: ndata_1d_dyn_imp=1
      integer,parameter :: ndata_2d_dyn_imp=2
      integer,parameter :: ndata_3d_dyn_imp=6
      integer,parameter :: ndata_1d_dyn_exp=1
      integer,parameter :: ndata_2d_dyn_exp=2
      integer,parameter :: ndata_3d_dyn_exp=9
!
!-----------------------------------------------------
!***  number of fields moved from physics 
!-----------------------------------------------------
!
      integer,parameter :: ndata_1d_phy_imp=1
      integer,parameter :: ndata_2d_phy_imp=2
      integer,parameter :: ndata_3d_phy_imp=9
      integer,parameter :: ndata_1d_phy_exp=1
      integer,parameter :: ndata_2d_phy_exp=2
      integer,parameter :: ndata_3d_phy_exp=6
!
!----------------------------------------------------------------
!***  the names of the fields that will move through the coupler
!----------------------------------------------------------------
!
      character(esmf_maxstr),dimension(ndata_1d) :: datanames_1d        &
                                                     =(/'date'/)
!
      character(esmf_maxstr),dimension(ndata_2d) :: datanames_2d        &
                                                     =(/'hs'            &
                                                       ,'ps'            &
                                                            /)
!
      character(esmf_maxstr),dimension(ndata_3d) :: datanames_3d        &
                                                     =(/'t     '        &
                                                       ,'u     '        &
                                                       ,'v     '        &
                                                       ,'shum  '        &
                                                       ,'oz    '        &
                                                       ,'cld   '        &
                                                       ,'p     '        &
                                                       ,'dp    '        &
                                                       ,'dpdt  '        &
                                                            /)
!
!-----------------------------------------------------------------------
!
      end module module_export_import_data
!
!-----------------------------------------------------------------------
