##  This is an elementary grid component or say independent grid component ##
 F77     = mpxlf95_r
 F90     = mpxlf95_r
 EXEC    = gfs_dyn.a
#
 FINCS    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
 OPTS     = -qsuffix=cpp=f -O0 -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave
 OPTS90   = -qsuffix=cpp=f -O0 -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave 
 OPTS90A  = -qsuffix=cpp=f -O0 -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave 
#TRAPS    = -qinitauto=FF911299 -qflttrap=ov:zero:inv:en -qsigtrap
#DEBUG    = -g -qfullpath -qcheck -qextchk
#
#FFLAGM  = -g -qsuffix=cpp=f -O0 -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave -qinitauto=7FBFFFFF -qflttrap=enable:overflow:zerodivide:invalid -qsigtrap=xl__trcedump -qcheck -qextchk -NS2048

 FFLAG90  = $(OPTS90) $(FINCS) -qfree -NS2048
 FFLAG90A = $(OPTS90A) $(FINCS) -qfree -NS2048
 FFLAGS   = $(OPTS) $(TRAPS) -qfixed
 FFLAGX   = $(OPTS) $(TRAPS) -qfixed
 FFLAGIO  = $(OPTS) $(TRAPS) -qfixed
 FFLAGY   = $(OPTS) -qfixed
 FFLAGM   = $(OPTS) $(FINCS) $(TRAPS) $(DEBUG) -NS2048 -qfixed
 FFLAGSF  = -O0 -qnosave -q free=f90 -qcheck
 FFLAGSI  = -qnosave -O0 -q free=f90
 FFLAGB   = -qnosave -O0 -q fixed

 ESMFLIB  = /meso/save/wx20rv/ESMF_libs/esmf_test_310r/lib/libO/AIX.default.64.mpi.default
 ESMFLIBS = ${ESMFLIB}/libesmf.a ${ESMFLIB}/libnetcdf_stubs.a
#----------------------------------
# For haze machine.
# LDR     = ncepmpxlf -qsmp=noauto
#----------------------------------
 LDR     = mpxlf95_r -qsmp=noauto 
 LDFLAGS = -lessl_r -lmass -qsmp=noauto ${ESMFLIBS}
# LIBS = -lC -L /nwprod/lib/ -l w3_d -l bacio_4
  LIBS = -lC -L /nwprod/lib/ -l bacio_4 -L /global/save/wx23lu/NEMS/gfs/lib/ -l w3_d
 AR = ar
 ARFLAGS = -r

.SUFFIXES: .o .f .F .h
#
# *****************************************************************
#
OBJS0	= gfs_dynamics_err_msg_mod.o		\
	  gfs_dynamics_namelist_mod.o		\
	  gfs_dynamics_start_time_get_mod.o	\
	  gfs_dynamics_internal_state_mod.o	\
	  gfs_dynamics_grid_create_mod.o	\
	  gfs_dynamics_getcf_mod.o		\
	  gfs_dynamics_initialize_mod.o		\
	  gfs_dynamics_run_mod.o		\
	  gfs_dynamics_finalize_mod.o		\
	  gfs_dynamics_add_get_state_mod.o	\
	  gfs_dynamics_states_mod.o		\
	  gfs_dynamics_grid_comp_mod.o
	  


OBJ_MOD	= gfs_dyn_machine.o             \
	  gfs_dyn_physcons.o          	\
	  gfs_dyn_resol_def.o           \
	  gfs_dyn_gg_def.o              \
	  gfs_dyn_vert_def.o            \
	  sig_io.o              	\
	  gfs_dyn_date_def.o            \
	  gfs_dyn_layout1.o             \
	  namelist_dynamics_def.o       \
	  gfs_dyn_coordinate_def.o      \
	  gfs_dyn_tracer_const.o  	\
	  gfs_dyn_reduce_lons_grid_module.o \
	  gfs_dyn_matrix_sig_def.o      \
	  gfs_dyn_mpi_def.o             \
	  semi_lag_def.o        	\
	  gfs_dyn_deldifs_def.o         \
	  gfs_dyn_bfilt_def.o   	\
	  gfsio_module.o        	\
	  gfsio_def.o           

#	  sigio_module.o        	\
#	  sigio_r_module.o      	\

OBJS	= \
	get_cd_hyb.o get_am_bm_hyb.o\
	get_cd_hyb_gc.o get_am_bm_hyb_gc.o\
	setlats.o\
	getysk.o\
	getysk_gc.o\
	redist_lats_a.o\
	ifshuff.o\
	shuffle_grid.o\
	wrtout_dynamics.o\
	grid_collect.o\
	excha.o\
	info.o\
	four2grid_thread.o\
	hpmdummy.o\
	coundummy.o\
	damp_speed.o\
	deldifs.o\
	input_fields.o\
	delnpe.o\
	delnpo.o \
	dezouv.o\
	impadj_sig.locl.o\
	impadj_hyb.locl.o\
	impadj_hyb.locl_gc.o\
	dozeuv.o\
	epslon.o\
	filtereo.o \
	filter1eo.o \
	filter2eo.o \
	get_ls_node.o\
	get_lats_node_shuff_a_fix.o\
	getcon_dynamics.o \
	gfidi_sig.o\
	gfidi_hyb.o\
	gfidi_hyb_gc.o\
	gfidi_hyb_gc_h.o\
	glats.o\
	grid_to_spect.o\
	spect_to_grid.o\
	spect_to_gridxy.o\
	lon_lat_para.o\
	compns_dynamics.o\
	gozrineo.o\
	pln2eo.o \
	rms_spect.o \
	lats_shuff.o\
	uvoedz.o                     \
	mpi_quit.o                   \
	sicdif_sig.o                 \
	sicdif_hyb.o                 \
	sicdif_hyb_gc.o              \
	common_to_model_vars.o      \
        sig2press.o          	\
	hyb2press.o          	\
	gch2press.o          	\
	omega_sig.o          	\
	omega_hyb.o          	\
	omega_gch.o          	\
	do_dynamics_mod.o          \
	do_dynamics_one_loop.o       \
	do_dynamics_two_loop.o       \
	model_to_common_vars.o      \
	gather_times.o               \
	four2fln.o                   \
	sortrx.o                     \
	sumfln.o                     \
	sumder2.o                    \
	triseori.o                   \
	uveodz.o                     \
	ysminv.o                     \
	ver_sig.o       

OBJS_PORT	= \
	fftpack.o \
	four2grid_thread.fftpack.o \
	noblas.o\

OBJS_IO= \
	grid_fix.o         \
	gribit.o           \
	treadeo.gfsio.o    \
	grid_to_spect_inp.o    \
	spect_to_grid_inp.o    \
	setsig.o           \
	bafrio.o        


#
INCS = f_hpm.h mpi_inc.h function2

#
# *****************************************************************
#
gfs_dynamics : $(EXEC)

$(EXEC):  $(OBJ_MOD) $(OBJS0) $(OBJS) $(OBJS_IO) 
	$(AR) $(ARFLAGS) $(EXEC) $(OBJ_MOD) $(OBJS0) $(OBJS) $(OBJS_IO) 
	echo " $(EXEC) created!"

clean:
	rm -f $(OBJ_MOD) $(OBJS0) $(OBJS) $(OBJS_IO) *.mod *.a

.F.o:
	$(F77) $(FFLAGS) -c -d $< 
.f.o:
	$(F77) $(FFLAGS) -c $< 


#
# *****************************************************************
#

gfs_dyn_machine.o:	gfs_dyn_machine.f
		$(F77) $(FFLAGM) -c gfs_dyn_machine.f

#
# *****************************************************************
#
gfs_dyn_physcons.o:	gfs_dyn_physcons.f
		$(F77) $(FFLAG90) -c gfs_dyn_physcons.f

gfsio_def.o:	gfsio_def.f
		$(F77) $(FFLAGSI) -c gfsio_def.f

gfsio_module.o:	gfsio_module.f
		$(F77) $(FFLAGSI) -c gfsio_module.f

bafrio.o:	bafrio.f
		$(F77) $(FFLAGB) -c bafrio.f

read_fix.o:	read_fix.f
		$(F77) $(FFLAGX) -c read_fix.f

softcount.o:	softcount.f
		$(F77) $(FFLAGX) -c softcount.f

getaer.o:	getaer.f
		$(F77) $(FFLAGIO) -c getaer.f

gribit.o:	gribit.f
		$(F77) $(FFLAGIO) -c gribit.f

treadeo.gfsio.o:	treadeo.gfsio.f
		$(F77) $(FFLAGIO) -c treadeo.gfsio.f

spect_to_grid_inp.o:	spect_to_grid_inp.f
		$(F77) $(FFLAGIO) -c spect_to_grid_inp.f

grid_to_spec_inp.o:	grid_to_spec_inp.f
		$(F77) $(FFLAGIO) -c grid_to_spec_inp.f

wrtout_dynamics.o:	wrtout_dynamics.f
		$(F77) $(FFLAGIO) -c wrtout_dynamics.f

gfs_dynamics_err_msg_mod.o:	 gfs_dynamics_err_msg_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_err_msg_mod.f

gfs_dynamics_namelist_mod.o:	 gfs_dynamics_namelist_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_namelist_mod.f

gfs_dynamics_start_time_get_mod.o:	 gfs_dynamics_start_time_get_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_start_time_get_mod.f

gfs_dynamics_grid_create_mod.o:	 gfs_dynamics_grid_create_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_grid_create_mod.f

gfs_dynamics_getcf_mod.o:	 gfs_dynamics_getcf_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_getcf_mod.f

gfs_dynamics_internal_state_mod.o:	 gfs_dynamics_internal_state_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_internal_state_mod.f

gfs_dynamics_add_get_state_mod.o:	 gfs_dynamics_add_get_state_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_add_get_state_mod.f

gfs_dynamics_states_mod.o:	 gfs_dynamics_states_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_states_mod.f

gfs_dynamics_initialize_mod.o:	 gfs_dynamics_initialize_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_initialize_mod.f

gfs_dynamics_run_mod.o:	 gfs_dynamics_run_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_run_mod.f

gfs_dynamics_finalize_mod.o:	 gfs_dynamics_finalize_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_finalize_mod.f

gfs_dynamics_grid_comp_mod.o:	 gfs_dynamics_grid_comp_mod.f
		$(F77) $(FFLAG90) -c gfs_dynamics_grid_comp_mod.f

