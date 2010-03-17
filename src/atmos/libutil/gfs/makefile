####  this is dependent grd component, needs elemetary grid components ###
#
#
# *****************************************************************
#
include ../../../../build/compile_settings
MAKEFILE=makefile
TARGET= gfs_libutil.a

OBJS01  = \
	module_gfs_machine.o          \
	module_include_gfs.o            \
	module_gfs_mpi_def.o

OBJS02=\
	module_NEMSIO.o                 \
	module_gfs_physcons.o           \
	module_gfs_funcphys.o           \
	module_gfs_tropp.o

OBJS03  = \
	module_ERR_MSG.o           \
	module_DM_PARALLEL_GFS.o      \
	module_digital_filter_gfs.o   \
	module_export_import_data.o   \
	atmos_dyn_phy_cpl_comp_mod.o \
	atmos_phy_chem_cpl_comp_mod.o \
	atmos_chem_phy_cpl_comp_mod.o



OBJS01_F90 = ${OBJS01:.o=.F90}
OBJS02_F90 = ${OBJS02:.o=.F90}
OBJS03_F90 = ${OBJS03:.o=.F90}

#
#
# *****************************************************************
#
.SUFFIXES:  .o .F90 .h

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)

$(LIBDIR)/lib/$(TARGET): $(OBJS01) $(OBJS02) $(OBJS03)
	$(AR) $(ARFLAGS) $@ $?

$(OBJS01): %.o: %.F90
	$(F77) $(FFLAG_GFS_UTIL) -c $*.F90

$(OBJS02): %.o: %.F90
	$(F90) $(FFLAG_GFS_UTIL) -c $*.F90

$(OBJS03): %.o: %.F90
	$(F77) $(FFLAG_GFS_UTIL90) -c $*.F90


depend: ${OBJS01_F90} ${OBJS02_F90} ${OBJS03_F90}
	$(MAKEDEPEXEC)

clean:
	/bin/rm -f *.o 
	/bin/rm -f *.mod 
	/bin/rm -f *.lst 
	/bin/rm -f $(LIBDIR)/lib/$(TARGET)
	/bin/rm -f depend

#
# *****************************************************************
#
include depend
