include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../../share

LIBRARY  = libfim.a

MODULES  = fim_internal_state.o fim_grid_comp.o
MODULES_STUB  = fim_grid_comp_stub.o

DEPEND_FILES = ${MODULES:.o=.F90}

all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $(MODULES)

stub: $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_STUB)
	
$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS_FIM) $(UTILINCS) -c $*.f90

$(MODULES_STUB): %.o: %.f90
	$(FC) $(FFLAGS_FIM) $(UTILINCS) -c $*.f90

clean:
	$(RM) -f $(LIBRARY) *.o *.mod fim_grid_comp.f90 fim_grid_comp_stub.f90 fim_internal_state.f90 *.lst lm map depend

include ../../conf/make.rules

include depend
