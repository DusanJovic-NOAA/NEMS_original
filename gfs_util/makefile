####  this is dependent grd component, needs elemetary grid components ###
 F77     = mpxlf95_r
 EXEC    = gfs_util.a
#
 DIR_DYN = ../gfs_dyn
 DIR_PHY = ../gfs_phys
 GCINCS  = -I$(DIR_DYN) -I$(DIR_PHY) 
 FINC0    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
 FINCS    = $(FINC0) $(GCINCS)
 OPTS     = -qsuffix=cpp=f -O0 -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave 
 OPTS90   = -qsuffix=cpp=f -O0 -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave 
 OPTS90A  = -qsuffix=cpp=f -O0 -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave 
#TRAPS    = -qinitauto=FF911299 -qflttrap=ov:zero:inv:en -qsigtrap -qsource -qlist
#DEBUG    = -g -qfullpath -qcheck -qextchk -qsource -qlist
#
#FFLAGM  = -g -qsuffix=cpp=f -O0 -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave -qinitauto=7FBFFFFF -qflttrap=enable:overflow:zerodivide:invalid -qsigtrap=xl__trcedump -qcheck -qextchk -NS2048
 FFLAG90  = $(OPTS90) $(FINCS) -qfree -NS2048
 FFLAG90A = $(OPTS90A) $(FINCS) -qfree -NS2048
 FFLAGS   = $(OPTS) $(FINCS) $(TRAPS) -qfixed
 FFLAGX   = $(OPTS) $(TRAPS) -qfixed
 FFLAGIO  = $(OPTS) $(TRAPS) -qfixed
 FFLAGY   = $(OPTS) -qfixed
 FFLAGM   = $(OPTS) $(FINCS) $(TRAPS) $(DEBUG) -NS2048 -qfixed
 FFLAGSF  = -O0 -qnosave -q free=f90 -qcheck
 FFLAGSI  = -qnosave -O0 -q free=f90
 FFLAGB   = -qnosave -O0 -q fixed
 ESMFLIB  = /meso/save/wx20rv/ESMF_libs/esmf_test_310r/lib/libO/AIX.default.64.mpi.default
 ESMFLIBS = ${ESMFLIB}/libesmf.a 
# --------------------------------
# for haze machine
# LDR     = ncepmpxlf -qsmp=noauto 
# --------------------------------
 LDR     = mpxlf95_r -qsmp=noauto  
 LDFLAGS = -lessl_r -lmass -qsmp=noauto ${ESMFLIBS}
 LIBS = -lC -L /nwprod/lib/ -l w3_d -l bacio_4 

 
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
gfs_atmosphere: $(EXEC)

$(EXEC): $(OBJS0) 
	$(AR) $(ARFLAGS) $(EXEC) $(OBJS0) 
	echo " $(EXEC) is created! "
clean:
	rm -f *.o *.mod *.lst  $(EXEC)

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
