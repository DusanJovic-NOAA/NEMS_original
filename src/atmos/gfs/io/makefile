include ../../../conf/configure.nems

MAKEFILE=makefile

UTILINCS = -I../../share -I../libutil

LIBRARY = gfs_io.a

OBJ_IO	= \
	  module_WRITE_INTERNAL_STATE_GFS.o   \
	  module_GET_CONFIG_WRITE_GFS.o       \
	  module_WRITE_ROUTINES_GFS.o         \
	  module_WRITE_GRID_COMP_GFS.o 

DEPEND_FILES = ${OBJ_IO:.o=.F90}

all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(OBJ_IO)
	$(AR) $(ARFLAGS) $@ $?

$(OBJ_IO): %.o: %.f90
	$(FC) $(FFLAGS_GFS) $(UTILINCS) -c $*.f90

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst lm map depend

include ../../../conf/make.rules

include depend
