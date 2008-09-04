#### this is elementary grid component or independent grid component ####
 F77     = mpxlf_r
 F90     = mpxlf95_r
 EXEC	 = gfs_phys.a
#
 FINC0    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
 FINCS    = $(FINC0) -I/nwprod/lib/incmod/w3_d

 DEBUG    = -qcheck -qsource -qlist
#---------------
# For non-debug.
#---------------
FFLAGM_DEBUG  = -O2

 OPTS     = -qsuffix=cpp=f -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave $(FFLAGM_DEBUG)
 OPTS90   = -qsuffix=cpp=f -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave $(FFLAGM_DEBUG)
 OPTS90A  = -qsuffix=cpp=f -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave $(FFLAGM_DEBUG)

 FFLAG90  = $(OPTS90) $(FINCS) -qfree -NS2048
 FFLAG90A = $(OPTS90A) $(FINCS) -qfree -NS2048
 FFLAGS   = $(OPTS) 
 FFLAGX   = $(OPTS)
 FFLAGIO  = $(OPTS) 
 FFLAGY   = $(OPTS)
 FFLAGM   = $(OPTS) $(FINCS) $(FFLAGM_DEBUG) -NS2048
 FFLAGSF  = -qnosave -q free=f90 -qcheck $(FFLAGM_DEBUG)
 FFLAGSI  = -qnosave -q free=f90 $(FFLAGM_DEBUG)
 FFLAGB   = -qnosave -qfixed $(FFLAGM_DEBUG)

 AR = ar
 ARFLAGS = -r
 
.SUFFIXES: .o .f .F .h
#
# *****************************************************************
#
OBJS0   = gfs_physics_err_msg_mod.o            \
	  gfs_physics_sfc_flx_mod.o	       \
	  gfs_physics_sfc_flx_set_mod.o	       \
	  gfs_physics_namelist_mod.o           \
	  gfs_physics_start_time_get_mod.o     \
	  gfs_physics_internal_state_mod.o     \
	  gfs_physics_grid_create_mod.o        \
	  gfs_physics_getcf_mod.o              \
 	  gfs_physics_initialize_mod.o         \
	  gfs_physics_run_mod.o                \
	  gfs_physics_finalize_mod.o           \
	  gfs_physics_add_get_state_mod.o      \
	  gfs_physics_states_mod.o             \
	  gfs_physics_grid_comp_mod.o

OBJ_MOD	= machine.o		\
	  iounitdef.o           \
	  resol_def.o           \
	  layout1.o             \
	  vert_def.o            \
	  gg_def.o              \
	  date_def.o            \
	  physcons_v.o          \
	  funcphys_v.o          \
	  progtm_module.o       \
	  rascnvv2_v.o          \
	  namelist_physics_def.o        \
	  namelist_soilveg.o    \
	  coordinate_def.o      \
	  tracer_const_h.o  	\
	  reduce_lons_grid_module.o	\
	  module_bfmicrophysics.o	\
	  mersenne_twister.o	\
	  mpi_def.o		\
	  d3d_def.o        	\
	  sfcio_module.o        

OBJS	= \
	  glats_physics.o               \
	  getcon_physics.o		\
	  lon_lat_para.o		\
	  compns_physics.o		\
	  lats_shuff.o			\
	  setlats.o			\
	  get_lats_node_shuff_r_fix.o	\
          common_to_physics_vars.o	\
	  do_physics_one_step.o		\
          physics_to_common_vars.o	\
	  gcycle_v.o           		\
	  hpmdummy.o

OBJS_IO= \
	  wrtout_physics.o   \
          fix_fields.o       \
	  read_fix.o         \
	  gribit.o           \
	  wrt3d.o            \
	  wrt3d_hyb.o        \
	  para_fixio_w.o     \
	  wrtsfc.o           \
	  sfcsub.o           

OBJS_RAD	=              \
	radlw_param.o          \
	radlw_datatb.o         \
	radlw_main.o           \
	radsw_param.o          \
	radsw_datatb.o         \
	radsw_main.o           \
	radiation_astronomy.o  \
	radiation_aerosols.o   \
	radiation_gases.o      \
	radiation_clouds.o     \
	radiation_surface.o    \
	grrad.o                \
	gloopr.o


OBJS_PHY= \
	ozinterp.o          \
	ozphys_v.o          \
	gloopb.o            \
	gbphys_v.o          \
	gbphys_adv_hyb_gc.o \
	gbphys_adv_hyb_gc_h.o \
	dcyc2_v.o           \
	dcyc2_v.pre.rad.o   \
	set_soilveg.o       \
	sfc_drv.o           \
	sfc_land.o          \
	progt2_v.o          \
	sfc_sice.o          \
	sfc_ocean.o         \
	sfc_diff.o          \
	sfc_diag.o          \
	sflx.o              \
	moninp_v.o          \
	moninq_v.o          \
	gwdps_v.o           \
	gwdc.o           \
	sascnv_v.o          \
	sascnvn_v.o          \
	cnvcld_v.o          \
	shalcv_v.o          \
	shalcnv_v.o          \
	lrgsclr_v.o         \
	gscond_v.o          \
	precpd_v.o          \
	mstadb_v.o          \
	mstcnv_v.o          \
	get_prs_v.o         \
	gsmddrive_v.o 

SRC= 				\
	$(OBJS0:.o=.f) 		\
	$(OBJ_MOD:.o=.f) 	\
	$(OBJS:.o=.f) 		\
	$(OBJS_RAD:.o=.f) 	\
	$(OBJS_PHY:.o=.f) 	\
	$(OBJS_IO:.o=.f)

#
# *****************************************************************
#
                                                                                
$(EXEC):  $(OBJ_MOD) $(OBJS0) $(OBJS)  $(OBJS_RAD) $(OBJS_PHY) $(OBJS_IO)
	$(AR) $(ARFLAGS) $(EXEC) $(OBJ_MOD) $(OBJS0) $(OBJS)  $(OBJS_RAD) $(OBJS_PHY) $(OBJS_IO)
	echo " $(EXEC) created!"

clean:
	rm -f $(OBJ_MOD) $(OBJS0) $(OBJS) $(OBJS_RAD) $(OBJS_PHY) $(OBJS_IO) *.mod *.a

.F.o:
	$(F77) $(FFLAGS) -c -d $< 
	#$(F77) $(FFLAGS) -c -d -WF,-DCLR:${RASV} $< 
.f.o:
	$(F77) $(FFLAGS) -c $< 


#
# *****************************************************************
#
cnvcld_v.o:	cnvcld_v.f
		$(F77) $(FFLAGM) -c cnvcld_v.f

dcyc2_v.o:	dcyc2_v.f
		$(F77) $(FFLAGM) -c dcyc2_v.f

dcyc2_v.pre.rad.o:	dcyc2_v.pre.rad.f
		$(F77) $(FFLAGM) -c dcyc2_v.pre.rad.f

funcphys_subsx.o:	funcphys_subsx.f
		$(F77) $(FFLAGM) -c funcphys_subsx.f

gbphys_adv_hyb_gc.o:	gbphys_adv_hyb_gc.f
		$(F77) $(FFLAGM) -c gbphys_adv_hyb_gc.f

gbphys_adv_hyb_gc_h.o:	gbphys_adv_hyb_gc_h.f
		$(F77) $(FFLAGM) -c gbphys_adv_hyb_gc_h.f

gbphys_v.o:	gbphys_v.f
		$(F77) $(FFLAGM) -c gbphys_v.f

get_prs_v.o:	get_prs_v.f
		$(F77) $(FFLAGM) -c get_prs_v.f

gscond_v.o:	gscond_v.f
		$(F77) $(FFLAGM) -c gscond_v.f

gsmddrive_v.o:	gsmddrive_v.f
		$(F77) $(FFLAGM) -c gsmddrive_v.f

module_bfmicrophysics.o:	module_bfmicrophysics.f
		$(F77) $(FFLAGM) -c module_bfmicrophysics.f

gwdps_v.o:	gwdps_v.f
		$(F77) $(FFLAGM) -c gwdps_v.f

gwdc.o:		gwdc.f
		$(F77) $(FFLAGM) -c gwdc.f

lrgsclr_v.o:	lrgsclr_v.f
		$(F77) $(FFLAGM) -c lrgsclr_v.f

moninp_v.o:	moninp_v.f
		$(F77) $(FFLAGM) -c moninp_v.f

moninq_v.o:	moninq_v.f
		$(F77) $(FFLAGM) -c moninq_v.f

moninp_han_v.o:	moninp_han_v.f
		$(F77) $(FFLAGM) -c moninp_han_v.f

mstadb_v.o:	mstadb_v.f
		$(F77) $(FFLAGM) -c mstadb_v.f

mstcnv_v.o:	mstcnv_v.f
		$(F77) $(FFLAGM) -c mstcnv_v.f

ozinterp.o:	ozinterp.f
		$(F77) $(FFLAGM) -c ozinterp.f

ozphys_v.o:	ozphys_v.f
		$(F77) $(FFLAGM) -c ozphys_v.f

precpd_v.o:	precpd_v.f
		$(F77) $(FFLAGM) -c precpd_v.f

rascnvv2_v.o:	rascnvv2_v.f
		$(F77) $(FFLAGM) -c rascnvv2_v.f

sfc_sice.o:	sfc_sice.f
		$(F77) $(FFLAGM) -c sfc_sice.f

set_soilveg.o:	set_soilveg.f
		$(F77) $(FFLAGM) -c set_soilveg.f

namelist_soilveg.o:	namelist_soilveg.f
		$(F77) $(FFLAGM) -c namelist_soilveg.f

sfc_land.o:	sfc_land.f
		$(F77) $(FFLAGM) -c sfc_land.f
                                                                            
progt2_v.o:	progt2_v.f
		$(F77) $(FFLAGM) -c progt2_v.f

sfc_drv.o:	sfc_drv.f
		$(F77) $(FFLAGM) -c sfc_drv.f

sflx.o: 	sflx.f
		$(F77) $(FFLAGM) -c sflx.f

sfc_ocean.o:	sfc_ocean.f
		$(F77) $(FFLAGM) -c sfc_ocean.f

sfc_diff.o:	sfc_diff.f
		$(F77) $(FFLAGM) -c sfc_diff.f

sfc_diag.o:	sfc_diag.f
		$(F77) $(FFLAGM) -c sfc_diag.f

sascnv_v.o:	sascnv_v.f
		$(F77) $(FFLAGM) -c sascnv_v.f

sascnvn_v.o:	sascnvn_v.f
		$(F77) $(FFLAGM) -c sascnvn_v.f

shalcv_v.o:	shalcv_v.f
		$(F77) $(FFLAGM) -c shalcv_v.f

shalcnv_v.o:	shalcnv_v.f
		$(F77) $(FFLAGM) -c shalcnv_v.f

# *****************************************************************

radlw_param.o:	radlw_param.f
		$(F77) $(FFLAGM) -c radlw_param.f


radlw_datatb.o:	radlw_datatb.f
		$(F77) $(FFLAGM) -c radlw_datatb.f


radlw_main.o:	radlw_main.f
		$(F77) $(FFLAGM) -c radlw_main.f

radsw_param.o:	radsw_param.f
		$(F77) $(FFLAGM) -c radsw_param.f

radsw_datatb.o:	radsw_datatb.f
		$(F77) $(FFLAGM) -c radsw_datatb.f

radsw_main.o:	radsw_main.f
		$(F77) $(FFLAGM) -c radsw_main.f

radiation_astronomy.o:	radiation_astronomy.f
		$(F77) $(FFLAGM) -c radiation_astronomy.f

radiation_aerosols.o:	radiation_aerosols.f
		$(F77) $(FFLAGM) -c radiation_aerosols.f

radiation_gases.o:	radiation_gases.f
		$(F77) $(FFLAGM) -c radiation_gases.f

radiation_clouds.o:	radiation_clouds.f
		$(F77) $(FFLAGM) -c radiation_clouds.f

radiation_surface.o:	radiation_surface.f
		$(F77) $(FFLAGM) -c radiation_surface.f

grrad.o:	grrad.f
		$(F77) $(FFLAGM) -c grrad.f

progtm_module.o:	progtm_module.f
		$(F77) $(FFLAGM) -c progtm_module.f

machine.o:	machine.f
		$(F77) $(FFLAGM) -c machine.f

#
# *****************************************************************
#
physcons_v.o:	physcons_v.f
		$(F77) $(FFLAG90) -c physcons_v.f

iounitdef.o:	iounitdef.f
		$(F77) $(FFLAG90) -c iounitdef.f

funcphys_v.o:	funcphys_v.f
		$(F77) $(FFLAG90) -c funcphys_v.f

sfcio_module.o:	sfcio_module.f
		$(F77) $(FFLAGSF) -c sfcio_module.f

sigio_module.o:	sigio_module.f
		$(F77) $(FFLAGSI) -c sigio_module.f

gfsio_def.o:	gfsio_def.f
		$(F77) $(FFLAGSI) -c gfsio_def.f

gfsio_module.o:	gfsio_module.f
		$(F77) $(FFLAGSI) -c gfsio_module.f

sigio_r_module.o:	sigio_r_module.f
		$(F77) $(FFLAGSI) -c sigio_r_module.f

gloopr.o:	gloopr.f
		$(F77) $(FFLAGS) -c gloopr.f

gloopb.o:	gloopb.f
		$(F77) $(FFLAGM) -c gloopb.f

sfcsub.o:	sfcsub.f
		$(F77) $(FFLAGM) -c sfcsub.f

gcycle_v.o:	gcycle_v.f
		$(F77) $(FFLAGM) -c gcycle_v.f

getaer.o:	getaer.f
		$(F77) $(FFLAGIO) -c getaer.f

wrt3d.o:	wrt3d.f
		$(F77) $(FFLAGIO) -c wrt3d.f

wrt3d_hyb.o:	wrt3d_hyb.f
		$(F77) $(FFLAGIO) -c wrt3d_hyb.f

gribit.o:	gribit.f
		$(F77) $(FFLAGIO) -c gribit.f

wrtsfc.o:	wrtsfc.f
		$(F77) $(FFLAGIO) -c wrtsfc.f

para_fixio_w.o:		para_fixio_w.f
		$(F77) $(FFLAGIO) -c para_fixio_w.f

wrtout_physics.o:	wrtout_physics.f
		$(F77) $(FFLAGIO) -c wrtout_physics.f

gfs_physics_err_msg_mod.o:      gfs_physics_err_msg_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_err_msg_mod.f

gfs_physics_sfc_flx_mod.o:      gfs_physics_sfc_flx_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_sfc_flx_mod.f
                                                                                
gfs_physics_sfc_flx_set_mod.o:      gfs_physics_sfc_flx_set_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_sfc_flx_set_mod.f
                                                                                
gfs_physics_namelist_mod.o:     gfs_physics_namelist_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_namelist_mod.f
                                                                                
gfs_physics_start_time_get_mod.o:       gfs_physics_start_time_get_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_start_time_get_mod.f
                                                                                
gfs_physics_grid_create_mod.o:  gfs_physics_grid_create_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_grid_create_mod.f
                                                                                
gfs_physics_getcf_mod.o:        gfs_physics_getcf_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_getcf_mod.f
                
gfs_physics_internal_state_mod.o:       gfs_physics_internal_state_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_internal_state_mod.f
                
gfs_physics_add_get_state_mod.o:        gfs_physics_add_get_state_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_add_get_state_mod.f

gfs_physics_states_mod.o:       gfs_physics_states_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_states_mod.f
                                                                                
gfs_physics_initialize_mod.o:   gfs_physics_initialize_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_initialize_mod.f
                                                                                
gfs_physics_run_mod.o:  gfs_physics_run_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_run_mod.f
                                                                                
gfs_physics_finalize_mod.o:     gfs_physics_finalize_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_finalize_mod.f
                                                                                
gfs_physics_grid_comp_mod.o:    gfs_physics_grid_comp_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_grid_comp_mod.f
                                                                                

