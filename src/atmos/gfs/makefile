include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../share -Ilibutil -Iio

GOCARTINCS = -I../../chem/gocart/AIX/include/Chem_Base \
             -I../../chem/gocart/AIX/include/GMAO_mpeu \
             -I../../chem/gocart/AIX/include/GEOSchem_GridComp

LIBRARY = libgfs.a

DIR_LIBUTIL = libutil
DIR_IO      = io
DIR_DYN     = dyn
DIR_PHY     = phys
GCINCS      = -I$(DIR_DYN) -I$(DIR_PHY)


OPTS     = $(FFLAGS_GFS) $(R8) -O3  $(UTILINCS) $(GCINCS) $(GOCARTINCS)
FFLAG90  = $(OPTS) $(FREE)

OBJS   = module_GFS_CORE_SETUP.o module_GFS_WRITE.o module_GOCART_ROUTINES.o \
	 module_GFS_INTERNAL_STATE.o module_GFS_INTEGRATE.o module_GFS_GRID_COMP.o

DEPEND_FILES = ${OBJS:.o=.F90}

MODULES_STUB  = module_GFS_GRID_COMP_stub.o

all default: depend
	cd $(DIR_LIBUTIL) && gmake && cd ..
	cd $(DIR_IO) && gmake && cd ..
	cd $(DIR_DYN) && gmake && cd ..
	cd $(DIR_PHY) && gmake && cd ..
	@gmake -f makefile $(LIBRARY)

$(LIBRARY): $(OBJS) 
	$(AR) $(ARFLAGS) $@ $?

stub: $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_STUB)

$(OBJS): %.o: %.f90
	$(FC) $(FFLAG90) -c $*.f90

$(MODULES_STUB): %.o: %.F90
	$(FC) $(FFLAGS_GFS) $(UTILINCS) -c $*.F90

clean:
	cd $(DIR_LIBUTIL) && gmake clean && cd ..
	cd $(DIR_IO) && gmake clean && cd ..
	cd $(DIR_DYN) && gmake clean && cd ..
	cd $(DIR_PHY) && gmake clean && cd ..
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst lm map depend

include ../../conf/make.rules

include depend
