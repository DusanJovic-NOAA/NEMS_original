####  this is dependent grd component, needs elemetary grid components ###
#
#
# *****************************************************************
#
include ../../ush/compile_settings
MAKEFILE=makefile
TARGET= gfs_libutil.a

OBJ=\
	module_nemsio.o

OBJS01  = \
	module_gfs_machine.o          \
	module_include_gfs.o            \
	module_gfs_mpi_def.o          

OBJS02  = \
	module_ERR_MSG.o           \
	module_DM_PARALLEL_GFS.o      \
	module_digital_filter_gfs.o   \
	module_export_import_data.o   \
	atmos_dyn_phy_cpl_comp_mod.o \
	atmos_dyn_chem_cpl_comp_mod.o \
	atmos_phy_chem_cpl_comp_mod.o



#
#
# *****************************************************************
#
.SUFFIXES:  .o .F90 .h

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)

$(LIBDIR)/lib/$(TARGET): $(OBJ) $(OBJS01) $(OBJS02)
	$(AR) $(ARFLAGS) $@ $?

$(OBJ): %.o: %.F90
	$(F90) $(FFLAG_GFS_UTIL) -c  $*.F90

$(OBJS01): %.o: %.F90
	$(FC) $(FFLAG_GFS_UTIL) -c $*.F90

$(OBJS02): %.o: %.F90
	$(FC) $(FFLAG_GFS_UTIL90) -c $*.F90

clean:
	/bin/rm -f *.o 
	/bin/rm -f *.mod 
	/bin/rm -f *.lst 
	/bin/rm -f $(LIBDIR)/lib/$(TARGET)

#
# *****************************************************************
