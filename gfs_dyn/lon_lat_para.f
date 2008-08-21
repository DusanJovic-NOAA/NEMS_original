      SUBROUTINE LONLAT_PARA(global_lats_a,XLON,XLAT,lonsperlat)
!
c***********************************************************************
!
      USE gfs_dyn_MACHINE , ONLY : kind_grid

      use gfs_dyn_resol_def
      use gfs_dyn_layout1
      use gfs_dyn_gg_def
      use gfs_dyn_physcons, pi => con_pi
      implicit none
      integer i,j,lat
      integer                 lonsperlat(latg)
      real (kind=kind_grid) tpi,hpi,bphi
      PARAMETER (TPI=2.E0*PI,HPI=0.5E0*PI)
      integer              global_lats_a(latg)
      real (kind=kind_grid) XLON(lonf,lats_node_a)
      real (kind=kind_grid) XLAT(lonf,lats_node_a)
!
      xlon=0.
      xlat=0.
 
      DO j=1,lats_node_a
        lat = global_lats_a(ipt_lats_node_a-1+j)
        BPHI = TPI/lonsperlat(lat)
        if (lat.le.latg2) then
          DO i=1,lonsperlat(lat)
            XLON(I,J) = (i-1) * BPHI
            XLAT(I,J) = HPI - colrad_a(lat)
          ENDDO
        else
          DO i=1,lonsperlat(lat)
            XLON(I,J) =  (i-1) * BPHI
!           XLAT(I,J) = colrad_a(lat)-HPI
            XLAT(I,J) = colrad_a(latg+1-lat)-HPI
          ENDDO
        endif
      ENDDO
 
      RETURN
      END
 
