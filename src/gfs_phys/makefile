#### this is elementary grid component or independent grid component ####
#---------------
# For non-debug.
#---------------
#
include ../../ush/compile_settings
MAKEFILE=makefile
TARGET= gfs_phys.a
# *****************************************************************
#
OBJS0   = \
	  gfs_physics_err_msg_mod.o            \
	  gfs_physics_sfc_flx_mod.o	       \
	  gfs_physics_sfc_flx_set_mod.o	       \
	  gfs_physics_gridgr_mod.o      	\
	  gfs_physics_namelist_mod.o           \
	  gfs_physics_start_time_get_mod.o     \
	  gfs_physics_internal_state_mod.o     \
	  gfs_physics_grid_create_mod.o        \
	  gfs_physics_output.o        \
	  gfs_physics_getcf_mod.o              \
 	  gfs_physics_initialize_mod.o         \
	  gfs_physics_run_mod.o                \
	  gfs_physics_finalize_mod.o           \
	  gfs_physics_add_get_state_mod.o      \
	  gfs_phy_states_mod.o             \
	  gfs_physics_grid_comp_mod.o

OBJ_MOD	= \
	  machine.o		\
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
	  gfs_phy_tracer_config.o  	\
	  reduce_lons_grid_module.o	\
	  module_bfmicrophysics.o	\
	  mersenne_twister.o	\
	  mpi_def.o		\
	  d3d_def.o        	\
	  sfcio_module.o        \
          mod_state.o    

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
	  fix_fields.o       \
	  read_fix.o         \
	  wrtout_physics.o   \
	  gribit.o           \
	  para_fixio_w.o     \
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


OBJS_PHY = \
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

SRC = 	\
	$(OBJS0:.o=.f) 		\
	$(OBJ_MOD:.o=.f) 	\
	$(OBJS:.o=.f) 		\
	$(OBJS_RAD:.o=.f) 	\
	$(OBJS_PHY:.o=.f) 	\
	$(OBJS_IO:.o=.f)

#
# *****************************************************************
#

.SUFFIXES:    .o .f

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)


$(LIBDIR)/lib/$(TARGET):  $(OBJ_MOD) $(OBJS0) $(OBJS) $(OBJS_RAD) $(OBJS_PHY) $(OBJS_IO)
	$(AR) $(ARFLAGS) $@ $?

.f.o: %.o: %.f
	$(FC) $(FFLAGS) -c $*.f

clean:
	/bin/rm -f $(OBJ_MOD) 
	/bin/rm -f $(OBJS0) 
	/bin/rm -f $(OBJS) 
	/bin/rm -f $(OBJS_RAD) 
	/bin/rm -f $(OBJS_PHY) 
	/bin/rm -f *.mod 
	/bin/rm -f *.o
	/bin/rm -f $(LIBDIR)/lib/$(TARGET)

#
# *****************************************************************
#
cnvcld_v.o:	cnvcld_v.f
		$(FC) $(FFLAGM) -c cnvcld_v.f

dcyc2_v.o:	dcyc2_v.f
		$(FC) $(FFLAGM) -c dcyc2_v.f

dcyc2_v.pre.rad.o:	dcyc2_v.pre.rad.f
		$(FC) $(FFLAGM) -c dcyc2_v.pre.rad.f

funcphys_subsx.o:	funcphys_subsx.f
		$(FC) $(FFLAGM) -c funcphys_subsx.f

gbphys_adv_hyb_gc.o:	gbphys_adv_hyb_gc.f
		$(FC) $(FFLAGM) -c gbphys_adv_hyb_gc.f

gbphys_adv_hyb_gc_h.o:	gbphys_adv_hyb_gc_h.f
		$(FC) $(FFLAGM) -c gbphys_adv_hyb_gc_h.f

gbphys_v.o:	gbphys_v.f
		$(FC) $(FFLAGM) -c gbphys_v.f

get_prs_v.o:	get_prs_v.f
		$(FC) $(FFLAGM) -c get_prs_v.f

gscond_v.o:	gscond_v.f
		$(FC) $(FFLAGM) -c gscond_v.f

gsmddrive_v.o:	gsmddrive_v.f
		$(FC) $(FFLAGM) -c gsmddrive_v.f

module_bfmicrophysics.o:	module_bfmicrophysics.f
		$(FC) $(FFLAGM) -c module_bfmicrophysics.f

gwdps_v.o:	gwdps_v.f
		$(FC) $(FFLAGM) -c gwdps_v.f

gwdc.o:		gwdc.f
		$(FC) $(FFLAGM) -c gwdc.f

lrgsclr_v.o:	lrgsclr_v.f
		$(FC) $(FFLAGM) -c lrgsclr_v.f

moninp_v.o:	moninp_v.f
		$(FC) $(FFLAGM) -c moninp_v.f

moninq_v.o:	moninq_v.f
		$(FC) $(FFLAGM) -c moninq_v.f

moninp_han_v.o:	moninp_han_v.f
		$(FC) $(FFLAGM) -c moninp_han_v.f

mstadb_v.o:	mstadb_v.f
		$(FC) $(FFLAGM) -c mstadb_v.f

mstcnv_v.o:	mstcnv_v.f
		$(FC) $(FFLAGM) -c mstcnv_v.f

ozinterp.o:	ozinterp.f
		$(FC) $(FFLAGM) -c ozinterp.f

ozphys_v.o:	ozphys_v.f
		$(FC) $(FFLAGM) -c ozphys_v.f

precpd_v.o:	precpd_v.f
		$(FC) $(FFLAGM) -c precpd_v.f

rascnvv2_v.o:	rascnvv2_v.f
		$(FC) $(FFLAGM) -c rascnvv2_v.f

sfc_sice.o:	sfc_sice.f
		$(FC) $(FFLAGM) -c sfc_sice.f

set_soilveg.o:	set_soilveg.f
		$(FC) $(FFLAGM) -c set_soilveg.f

namelist_soilveg.o:	namelist_soilveg.f
		$(FC) $(FFLAGM) -c namelist_soilveg.f

sfc_land.o:	sfc_land.f
		$(FC) $(FFLAGM) -c sfc_land.f
                                                                            
progt2_v.o:	progt2_v.f
		$(FC) $(FFLAGM) -c progt2_v.f

sfc_drv.o:	sfc_drv.f
		$(FC) $(FFLAGM) -c sfc_drv.f

sflx.o: 	sflx.f
		$(FC) $(FFLAGM) -c sflx.f

sfc_ocean.o:	sfc_ocean.f
		$(FC) $(FFLAGM) -c sfc_ocean.f

sfc_diff.o:	sfc_diff.f
		$(FC) $(FFLAGM) -c sfc_diff.f

sfc_diag.o:	sfc_diag.f
		$(FC) $(FFLAGM) -c sfc_diag.f

sascnv_v.o:	sascnv_v.f
		$(FC) $(FFLAGM) -c sascnv_v.f

sascnvn_v.o:	sascnvn_v.f
		$(FC) $(FFLAGM) -c sascnvn_v.f

shalcv_v.o:	shalcv_v.f
		$(FC) $(FFLAGM) -c shalcv_v.f

shalcnv_v.o:	shalcnv_v.f
		$(FC) $(FFLAGM) -c shalcnv_v.f

# *****************************************************************

radlw_param.o:	radlw_param.f
		$(FC) $(FFLAGM) -c radlw_param.f

radlw_datatb.o:	radlw_datatb.f
		$(FC) $(FFLAGM) -c radlw_datatb.f

radlw_main.o:	radlw_main.f
		$(FC) $(FFLAGM) -c radlw_main.f

radsw_param.o:	radsw_param.f
		$(FC) $(FFLAGM) -c radsw_param.f

radsw_datatb.o:	radsw_datatb.f
		$(FC) $(FFLAGM) -c radsw_datatb.f

radsw_main.o:	radsw_main.f
		$(FC) $(FFLAGM) -c radsw_main.f

radiation_astronomy.o:	radiation_astronomy.f
		$(FC) $(FFLAGM) -c radiation_astronomy.f

radiation_aerosols.o:	radiation_aerosols.f
		$(FC) $(FFLAGM) -c radiation_aerosols.f

radiation_gases.o:	radiation_gases.f
		$(FC) $(FFLAGM) -c radiation_gases.f

radiation_clouds.o:	radiation_clouds.f
		$(FC) $(FFLAGM) -c radiation_clouds.f

radiation_surface.o:	radiation_surface.f
		$(FC) $(FFLAGM) -c radiation_surface.f

grrad.o:	grrad.f
		$(FC) $(FFLAGM) -c grrad.f

progtm_module.o:	progtm_module.f
		$(FC) $(FFLAGM) -c progtm_module.f

machine.o:	machine.f
		$(FC) $(FFLAGM) -c machine.f

mpi_def.o:	mpi_def.f
		$(FC) $(FFLAG_GFS_PHY) -c mpi_def.f

#
# *****************************************************************
#
physcons_v.o:	physcons_v.f
		$(FC) $(FFLAG90) -c physcons_v.f

iounitdef.o:	iounitdef.f
		$(FC) $(FFLAG90) -c iounitdef.f

funcphys_v.o:	funcphys_v.f
		$(FC) $(FFLAG90) -c funcphys_v.f

sfcio_module.o:	sfcio_module.f
		$(FC) $(FFLAGSF) -c sfcio_module.f

sigio_module.o:	sigio_module.f
		$(FC) $(FFLAGSI) -c sigio_module.f

gfsio_def.o:	gfsio_def.f
		$(FC) $(FFLAGSI) -c gfsio_def.f

gfsio_module.o:	gfsio_module.f
		$(FC) $(FFLAGSI) -c gfsio_module.f

sigio_r_module.o:	sigio_r_module.f
		$(FC) $(FFLAGSI) -c sigio_r_module.f

gloopr.o:	gloopr.f
		$(FC) $(FFLAGS) -c gloopr.f

gloopb.o:	gloopb.f
		$(FC) $(FFLAGM) -c gloopb.f

sfcsub.o:	sfcsub.f
		$(FC) $(FFLAGM) -c sfcsub.f

gcycle_v.o:	gcycle_v.f
		$(FC) $(FFLAGM) -c gcycle_v.f

getaer.o:	getaer.f
		$(FC) $(FFLAGIO) -c getaer.f

gribit.o:	gribit.f
		$(FC) $(FFLAGIO) -c gribit.f

para_fixio_w.o:  para_fixio_w.f
		$(FC) $(FFLAGIO) -c para_fixio_w.f

wrtout_physics.o: wrtout_physics.f
		$(FC) $(FFLAG_GFS_PHY) -c wrtout_physics.f

gfs_physics_err_msg_mod.o:      gfs_physics_err_msg_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_err_msg_mod.f

gfs_physics_gridgr_mod.o:	gfs_physics_gridgr_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_gridgr_mod.f

gfs_physics_sfc_flx_mod.o:      gfs_physics_sfc_flx_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_sfc_flx_mod.f
                                                                                
gfs_physics_sfc_flx_set_mod.o:      gfs_physics_sfc_flx_set_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_sfc_flx_set_mod.f
                                                                                
gfs_physics_namelist_mod.o:     gfs_physics_namelist_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_namelist_mod.f
                                                                                
gfs_physics_start_time_get_mod.o:       gfs_physics_start_time_get_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_start_time_get_mod.f
                                                                                
gfs_physics_grid_create_mod.o:  gfs_physics_grid_create_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_grid_create_mod.f
                                                                                
gfs_physics_output.o:   gfs_physics_output.f
		$(FC) $(FFLAG90) -c gfs_physics_output.f

gfs_physics_getcf_mod.o:        gfs_physics_getcf_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_getcf_mod.f
                
gfs_physics_internal_state_mod.o:       gfs_physics_internal_state_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_internal_state_mod.f
                
gfs_physics_add_get_state_mod.o:        gfs_physics_add_get_state_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_add_get_state_mod.f

gfs_phy_states_mod.o:		gfs_phy_states_mod.f
		$(FC) $(FFLAG90) -c gfs_phy_states_mod.f

gfs_physics_initialize_mod.o:   gfs_physics_initialize_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_initialize_mod.f
                                                                                
gfs_physics_run_mod.o:  gfs_physics_run_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_run_mod.f
                                                                                
gfs_physics_finalize_mod.o:     gfs_physics_finalize_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_finalize_mod.f
                                                                                
gfs_physics_grid_comp_mod.o:    gfs_physics_grid_comp_mod.f
		$(FC) $(FFLAG90) -c gfs_physics_grid_comp_mod.f

### use FFLAGSDB compile option
gfs_phy_tracer_config.o:	gfs_phy_tracer_config.f
		$(FC)  $(FFLAGSDB) -c gfs_phy_tracer_config.f
               
