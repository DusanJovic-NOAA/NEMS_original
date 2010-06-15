include ../../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../../share
GOCARTINCS = -I../../../chem/gocart/AIX/include/Chem_Base \
             -I../../../chem/gocart/AIX/include/GMAO_mpeu \
             -I../../../chem/gocart/AIX/include/GEOSchem_GridComp

LIBRARY  = gfs_libutil.a

OPTS     = -O3

FFLAG_GFS_UTIL90 = $(FFLAGS_GFS) $(OPTS) $(R8) $(FREE)  $(UTILINCS) $(GOCARTINCS)

OBJS02  = \
	module_gfs_machine.o          \
	module_include_gfs.o          \
	module_gfs_mpi_def.o          \
	module_gfs_physcons.o         \
	module_gfs_funcphys.o         \
	module_gfs_tropp.o            \
	module_DM_PARALLEL_GFS.o      \
	module_digital_filter_gfs.o   \
	module_export_import_data.o   \
	gfs_dyn_phy_cpl_comp_mod.o    \
	atmos_phy_chem_cpl_comp_mod.o \
	atmos_chem_phy_cpl_comp_mod.o


OBJS02_F90 = ${OBJS02:.o=.F90}

#
# *****************************************************************
#

.SUFFIXES:  .o .F90 .h

default:
	@gmake -f makefile all

all: depend $(LIBRARY)

$(LIBRARY): $(OBJS02)
	$(AR) $(ARFLAGS) $@ $?

$(OBJS02): %.o: %.F90
	$(FC) $(FFLAG_GFS_UTIL90) -c $*.F90

depend: ${OBJS02_F90}
	$(MAKEDEPF90) -nosrc $^ > depend

clean:
	$(RM) -f *.o 
	$(RM) -f *.mod 
	$(RM) -f *.lst 
	$(RM) -f $(LIBRARY)
	$(RM) -f depend

#
# DEPENDENCIES :
#
#
include depend
