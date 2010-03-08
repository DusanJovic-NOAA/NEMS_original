####  this is dependent grd component, needs elemetary grid components ###
#
#
# *****************************************************************
#
include ../../ush/compile_settings
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
	atmos_phy_chem_cpl_comp_mod.o



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

.SUFFIXES:      .F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

$(OBJS01): %.o: %.f90 $(DEPS)
	$(F77) $(FFLAG_GFS_UTIL) -c  $*.f90

$(OBJS02): %.o: %.f90 $(DEPS)
	$(F90) $(FFLAG_GFS_UTIL90) -c $*.f90

$(OBJS03): %.o: %.f90 $(DEPS)
	$(F90) $(FFLAG_GFS_UTIL90) -c $*.f90

clean:
	/bin/rm -f *.o 
	/bin/rm -f *.f90
	/bin/rm -f *.mod 
	/bin/rm -f *.lst 
	/bin/rm -f $(LIBDIR)/lib/$(TARGET)

#
# *****************************************************************
