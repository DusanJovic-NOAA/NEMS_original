       module gfs_dyn_dfi_mod
!
!*** jw: set up pointers for digital filter variables
!
        use gfs_dyn_machine, only:kind_evod
        implicit none
!
        type gfs_dfi_grid_gr

          integer z_imp,ps_imp,temp_imp,u_imp,v_imp,tracer_imp,         &
     &          p_imp,dp_imp,dpdt_imp
          real(kind=kind_evod),pointer:: hs(:,:,:)=>null()
          real(kind=kind_evod),pointer:: ps(:,:,:)=>null()
          real(kind=kind_evod),pointer:: t(:,:,:)=>null()
          real(kind=kind_evod),pointer:: u(:,:,:)=>null()
          real(kind=kind_evod),pointer:: v(:,:,:)=>null()
          real(kind=kind_evod),pointer:: tracer(:,:,:)=>null()
          real(kind=kind_evod),pointer:: p(:,:,:)=>null()
          real(kind=kind_evod),pointer:: dp(:,:,:)=>null()
          real(kind=kind_evod),pointer:: dpdt(:,:,:)=>null()

        end type gfs_dfi_grid_gr
!
       end module gfs_dyn_dfi_mod
