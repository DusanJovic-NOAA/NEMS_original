      subroutine sig2press(njeff,nsize_ar,pgr,prsl,dprs)
 
      use gfs_dyn_machine , only : kind_grid
      use gfs_dyn_resol_def
      use gfs_dyn_vert_def
      use gfs_dyn_coordinate_def
      implicit none
 
      integer njeff,nsize_ar
      real(kind=kind_grid)  pgr(nsize_ar)
      real(kind=kind_grid) prsl(nsize_ar,levs)
      real(kind=kind_grid) dprs(nsize_ar,levs)

      integer i,k
 
      do k=1,levs
        do i=1,njeff
          prsl(i,k)  = sl(k)*pgr(i)
          dprs(i,k)  = (si(k)-si(k+1))*pgr(i)
        enddo
      enddo
 
      return
      end
