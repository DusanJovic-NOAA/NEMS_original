include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../../share

LIBRARY       = libfim.a

MODULES  = fim_internal_state.o fim_grid_comp.o
MODULES_STUB  = fim_grid_comp_stub.o

MODULES_f90 = ${MODULES:.o=.f90}

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

all: depend $(MODULES)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES)
	
stub: depend $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_STUB)
	
$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS_FIM) $(UTILINCS) -c $*.f90

$(MODULES_STUB): %.o: %.f90
	$(FC) $(FFLAGS_FIM) $(UTILINCS) -c $*.f90

depend: $(MODULES_f90)
	makedepf90 -nosrc $^ > depend

clean:
	$(RM) -f $(LIBRARY) *.o *.mod *.lst lm map

#
# DEPENDENCIES :
#
include depend
