include ../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../atmos/share  -Ihycom  -Igenocn

LIBRARY  = libocean.a

MODULES  = module_OCN_INTERNAL_STATE.o module_OCN_GRID_COMP.o

MODULES_stub  = module_OCN_GRID_COMP_stub.o

DEPEND_FILES = ${MODULES:.o=.F90}

DEPEND_FILES_stub = ${MODULES_stub:.o=.F90}

all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

stub: depend_stub
	@gmake -f $(MAKEFILE) $(LIBRARY_stub)

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $?
	
$(LIBRARY_stub): $(MODULES_stub)
	$(AR) $(ARFLAGS) $@ $?
	
$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

$(MODULES_stub): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

clean:
	cd genocn        && gmake clean
	cd hycom         && gmake clean
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst depend

include ../conf/make.rules

include depend
