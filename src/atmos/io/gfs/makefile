##  This is an elementary grid component or say independent grid component ##
#---------------
# For non-debug.
#---------------
#
include ../../../../build/compile_settings
MAKEFILE=makefile
TARGET= gfs_io.a
# *****************************************************************
#
OBJ_IO	= \
	  module_WRITE_INTERNAL_STATE_GFS.o   \
	  module_GET_CONFIG_WRITE_GFS.o       \
	  module_WRITE_ROUTINES_GFS.o         \
	  module_WRITE_GRID_COMP_GFS.o 

# *****************************************************************
#
OBJ_IO_f90 = ${OBJ_IO:.o=.f90}

.SUFFIXES: .F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)

# *****************************************************************
#

$(LIBDIR)/lib/$(TARGET):  $(OBJ_IO)
	$(AR) $(ARFLAGS) $@ $?

$(OBJ_IO): %.o: %.f90  $(DEPS)
	$(F77) $(FFLAG_GFS_IO) -c $*.f90

depend: $(OBJ_IO_f90)
	$(MAKEDEPEXEC)

clean:
	/bin/rm -f $(OBJS_IO) 
	/bin/rm -f *.mod 
	/bin/rm -f *.f90
	/bin/rm -f *.o
	/bin/rm -f $(LIBDIR)/lib/$(TARGET)
	/bin/rm -f depend

#
# DEPENDENCIES :
#
#
include depend
