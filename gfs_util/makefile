####  this is dependent grd component, needs elemetary grid components ###
 F77     = mpxlf_r
 EXEC    = gfs_util.a
#
 DIR_DYN = ../gfs_dyn
 DIR_PHY = ../gfs_phys
 GCINCS  = -I$(DIR_DYN) -I$(DIR_PHY) 
#
 FINC0    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
 FINCS    = $(FINC0) $(GCINCS)
 DEBUG    = -qcheck -qsource -qlist
#FFLAGM_DEBUG  = -g -fullpath -O0 -qinitauto=FF911299 -qflttrap=enable:overflow:zerodivide:invalid -qsigtrap=xl__trcedump -NS2048 $(DEBUG)

#---------------
# For non-debug.
#---------------
FFLAGM_DEBUG  = -O0

 OPTS     = -qsuffix=cpp=f -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave $(FFLAGM_DEBUG)
 OPTS90   = -qsuffix=cpp=f -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave $(FFLAGM_DEBUG)
 OPTS90A  = -qsuffix=cpp=f -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave $(FFLAGM_DEBUG)
#TRAPS    = -qinitauto=FF911299 -qflttrap=ov:zero:inv:en -qsigtrap -qsource -qlist

#
 FFLAG90  = $(OPTS90) $(FINCS) -qfree -NS2048
 FFLAG90A = $(OPTS90A) $(FINCS) -qfree -NS2048
 FFLAGS   = $(OPTS) $(FINCS) $(TRAPS)
 FFLAGX   = $(OPTS) $(TRAPS)
 FFLAGIO  = $(OPTS) $(TRAPS)
 FFLAGY   = $(OPTS)
 FFLAGM   = $(OPTS) $(FINCS) $(TRAPS) $(FFLAGM_DEBUG) -NS2048
 FFLAGSF  = -qnosave -q free=f90 -qcheck $(FFLAGM_DEBUG)
 FFLAGSI  = -qnosave -q free=f90 $(FFLAGM_DEBUG)
 FFLAGB   = -qnosave -qfixed $(FFLAGM_DEBUG)

 AR = ar
 ARFLAGS = -r

.SUFFIXES: .o .f .F90 .h
#
# *****************************************************************
#
OBJS0	= atmos_err_msg_mod.o 		\
	  module_digital_filter_gfs.o 	\
	  module_include.o 		\
	  module_export_import_data.o 	\
	  atmos_dyn_phy_cpl_comp_mod.o \
	  atmos_grid_comp_mod.o 

#
#
# *****************************************************************
#

$(EXEC): $(OBJS0) 
	$(AR) $(ARFLAGS) $(EXEC) $(OBJS0) 
	echo " $(EXEC) is created! "
clean:
	rm -f *.o *.mod *.lst *.a $(EXEC)

.F90.o:
	$(F77) $(FFLAGS) -c -d $< 
.f.o:
	$(F77) $(FFLAGS) -c $< 

#
# *****************************************************************
atmos_err_msg_mod.o:     atmos_err_msg_mod.F90
	$(F77) $(FFLAG90) -c atmos_err_msg_mod.F90
                                                                                  
module_digital_filter_gfs.o:     module_digital_filter_gfs.F90
	$(F77) $(FFLAG90) -c module_digital_filter_gfs.F90
                                                                                  
module_include.o:     module_include.F90
	$(F77) $(FFLAG90) -c module_include.F90
                                                                                  
module_export_import_data.o:     module_export_import_data.F90
	$(F77) $(FFLAG90) -c module_export_import_data.F90
                                                                                  
atmos_dyn_phy_cpl_comp_mod.o:     atmos_dyn_phy_cpl_comp_mod.F90
	$(F77) $(FFLAG90) -c atmos_dyn_phy_cpl_comp_mod.F90
                                                                                  
atmos_grid_comp_mod.o:     atmos_grid_comp_mod.F90
	$(F77) $(FFLAG90) -c atmos_grid_comp_mod.F90
