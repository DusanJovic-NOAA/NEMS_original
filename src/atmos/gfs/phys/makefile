include ../../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../../share -I../libutil -I../../phys
ARCH=$(shell unmae -s)
GOCARTINCS = -I../../../chem/gocart/${ARCH}/include/Chem_Base \
             -I../../../chem/gocart/${ARCH}/include/GMAO_mpeu \
             -I../../../chem/gocart/${ARCH}/include/GEOSchem_GridComp


LIBRARY   = gfs_physics.a

DEBUG    = -qcheck -qsource -qlist
#FFLAGM_DEBUG  = -O0 -qcheck
FFLAGM_DEBUG  = -O3
OPTS     = -qsuffix=cpp=f -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto  -qnosave $(FFLAGM_DEBUG)

FFLAG90  = $(FFLAGS) $(OPTS) -qfree -NS2048 $(UTILINCS)
FFLAGS1  = $(FFLAGS) $(OPTS) $(TRAPS) -qfixed $(UTILINCS) $(GOCARTINCS)
FFLAGIO  = $(OPTS) $(TRAPS) $(UTILINCS) -qfixed
FFLAGM   = $(OPTS) $(FINCS) -qfixed $(TRAPS) $(FFLAGM_DEBUG) $(UTILINCS) -NS2048
FFLAGSF=-qnosave -q free=f90 -qcheck $(FFLAGM_DEBUG)
FFLAGSI  = -qrealsize=8 -qnosave -q free=f90 $(FFLAGM_DEBUG)
FFLAG_GFS_PHY=$(FFLAGS) $(OPTS) $(TRAPS) -qfixed    $(UTILINCS) $(GOCARTINCS)
FFLAG_GFS_PHY90=$(FFLAG90) $(TRAPS)     $(UTILINCS) $(GOCARTINCS)


OBJS0   =                                      \
	  gfs_physics_err_msg_mod.o            \
	  gfs_physics_sfc_flx_mod.o	       \
	  gfs_physics_nst_var_mod.o            \
	  gfs_physics_sfc_flx_set_mod.o	       \
	  gfs_physics_gridgr_mod.o             \
	  gfs_physics_g2d_mod.o                \
	  gfs_physics_g3d_mod.o      	       \
	  gfs_physics_namelist_mod.o           \
	  gfs_physics_start_time_get_mod.o     \
	  gfs_physics_internal_state_mod.o     \
	  gfs_physics_grid_create_mod.o        \
	  gfs_physics_output.o                 \
	  gfs_physics_getcf_mod.o              \
 	  gfs_physics_initialize_mod.o         \
	  gfs_physics_run_mod.o                \
	  gfs_physics_finalize_mod.o           \
	  gfs_physics_add_get_state_ESMFField_mod.o      \
	  gfs_phy_states_mod.o                 \
	  gfs_physics_grid_comp_mod.o

OBJ_MOD	=                           \
	  resol_def.o               \
	  layout1.o                 \
	  vert_def.o                \
	  gg_def.o                  \
	  namelist_physics_def.o    \
	  coordinate_def.o          \
	  reduce_lons_grid_module.o \
	  mpi_def.o		    \
	  d3d_def.o        	    \
	  sfcio_module.o            \
          mod_state.o               \
	  nstio_module.o

OBJS    =                               \
          glats_physics.o               \
          getcon_physics.o              \
          lon_lat_para.o                \
          compns_physics.o              \
          info.o                        \
          mymaxmin.o                    \
          lats_shuff.o                  \
          setlats.o                     \
          coundummy.o                   \
          sortrx.o                      \
          mpi_quit.o                    \
          get_lats_node_shuff_r_fix.o   \
          do_physics_one_step.o         \
          gcycle.o                      \
          dfi_fixwr.o                  \
          hpmdummy.o

OBJS_IO=                       \
	  fix_fields.o         \
	  read_fix.o           \
	  GFS_simple_scatter.o \
	  wrtout_physics.o     \
	  fld_collect.o        \
	  fld_collect_nst.o    \
	  para_fixio_w.o       \
	  para_nst_w.o

OBJS_RAD =                    \
	gloopr.o


OBJS_PHY = \
	gloopb.o              \
	gbphys_adv_hyb_gc.o   \
	gbphys_adv_hyb_gc_h.o

OBJS_CC=            \
	mpi_more.o  \
	cmp.comm.o  \
	atm.comm.o  \
	tiles.o

#
# *****************************************************************
#

DEPEND_FILES = ${OBJS0:.o=.f} ${OBJ_MOD:.o=.f} ${OBJS:.o=.f} ${OBJS_IO:.o=.f} ${OBJS_RAD:.o=.f} ${OBJS_PHY:.o=.f} $(OBJS_CC:.o=.f)

all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY):  $(OBJ_MOD) $(OBJS_CC) $(OBJS0) $(OBJS) $(OBJS_RAD) $(OBJS_PHY) $(OBJS_IO)
	$(AR) $(ARFLAGS) $@ $?

.f.o: %.o: %.f
	$(F77) $(FFLAGS1) -c $*.f

clean:
	$(RM) -f $(LIBRARY) *.o *.mod *.lst lm map depend

#
# *****************************************************************
#

gbphys_adv_hyb_gc.o:	gbphys_adv_hyb_gc.f
		$(F77) $(FFLAGM) -c gbphys_adv_hyb_gc.f

gbphys_adv_hyb_gc_h.o:	gbphys_adv_hyb_gc_h.f
		$(F77) $(FFLAGM) -c gbphys_adv_hyb_gc_h.f

# *****************************************************************

mpi_def.o:	mpi_def.f
		$(F77) $(FFLAG_GFS_PHY) -c mpi_def.f

#
# *****************************************************************
#

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
		$(F77) $(FFLAGS1) -c gloopr.f

gloopb.o:	gloopb.f
		$(F77) $(FFLAGM) -c gloopb.f

gcycle.o:	gcycle.f
		$(F77) $(FFLAGM) -c gcycle.f

getaer.o:	getaer.f
		$(F77) $(FFLAGIO) -c getaer.f

gribit.o:	gribit.f
		$(F77) $(FFLAGIO) -c gribit.f

fld_collect.o:  fld_collect.f
	        $(F77) $(FFLAG_GFS_PHY) -c fld_collect.f

fld_collect_nst.o:  fld_collect_nst.f
	        $(F77) $(FFLAG_GFS_PHY) -c fld_collect_nst.f

read_fix.o: read_fix.f
	        $(F77) $(FFLAG_GFS_PHY) -c read_fix.f

para_fixio_w.o:  para_fixio_w.f
		$(F77) $(FFLAG_GFS_PHY) -c para_fixio_w.f

para_nst_w.o:  para_nst_w.f
		$(F77) $(FFLAG_GFS_PHY) -c para_nst_w.f

mpi_more.o:	mpi_more.f
		$(F77) $(FFLAGxSM) -c mpi_more.f

cmp.comm.o:	cmp.comm.f
		$(F77) $(FFLAGM) -c cmp.comm.f

atm.comm.o:	atm.comm.f
		$(F77) $(FFLAGM) -c atm.comm.f

tiles.o:	tiles.f
		$(F77) $(FFLAGM) -c tiles.f

para_nstio_w.o:  para_nstio_w.f
		$(F77) $(FFLAG_GFS_PHY) -c para_nstio_w.f

wrtout_physics.o: wrtout_physics.f
		$(F77) $(FFLAG_GFS_PHY) -c wrtout_physics.f

nstio_module.o:	  nstio_module.f
		$(F77) $(FFLAG90) -c nstio_module.f

gfs_physics_err_msg_mod.o:      gfs_physics_err_msg_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_err_msg_mod.f

gfs_physics_gridgr_mod.o:	gfs_physics_gridgr_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_gridgr_mod.f

gfs_physics_g2d_mod.o:          gfs_physics_g2d_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_g2d_mod.f

gfs_physics_g3d_mod.o:		gfs_physics_g3d_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_g3d_mod.f

gfs_physics_sfc_flx_mod.o:      gfs_physics_sfc_flx_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_sfc_flx_mod.f

gfs_physics_nst_var_mod.o:      gfs_physics_nst_var_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_nst_var_mod.f
                                                                                
gfs_physics_sfc_flx_set_mod.o:      gfs_physics_sfc_flx_set_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_sfc_flx_set_mod.f

gfs_physics_nst_flx_var.o:	gfs_physics_nst_var_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_nst_var_mod.f
                                                                                
gfs_physics_namelist_mod.o:     gfs_physics_namelist_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_namelist_mod.f
                                                                                
gfs_physics_start_time_get_mod.o:       gfs_physics_start_time_get_mod.f
		$(F77) $(FFLAG_GFS_PHY90) -c gfs_physics_start_time_get_mod.f
                                                                                
gfs_physics_grid_create_mod.o:  gfs_physics_grid_create_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_grid_create_mod.f
                                                                                
gfs_physics_output.o:   gfs_physics_output.f
		$(F77) $(FFLAG90) -c gfs_physics_output.f

gfs_physics_getcf_mod.o:        gfs_physics_getcf_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_getcf_mod.f
                
gfs_physics_internal_state_mod.o:       gfs_physics_internal_state_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_internal_state_mod.f
                
gfs_physics_add_get_state_ESMFField_mod.o:        gfs_physics_add_get_state_ESMFField_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_add_get_state_ESMFField_mod.f

gfs_phy_states_mod.o:		gfs_phy_states_mod.f
		$(F77) $(FFLAG90) -c gfs_phy_states_mod.f

gfs_physics_initialize_mod.o:   gfs_physics_initialize_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_initialize_mod.f
                                                                                
gfs_physics_run_mod.o:  gfs_physics_run_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_run_mod.f
                                                                                
gfs_physics_finalize_mod.o:     gfs_physics_finalize_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_finalize_mod.f
                                                                                
gfs_physics_grid_comp_mod.o:    gfs_physics_grid_comp_mod.f
		$(F77) $(FFLAG90) -c gfs_physics_grid_comp_mod.f

#
#
#

include ../../../conf/make.rules

include depend
