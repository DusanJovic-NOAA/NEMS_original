include ../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../atmos/share  -Ihycom  -Igenocn

LIBRARY  = libocean.a

MODULES  = module_OCN_INTERNAL_STATE.o module_OCN_GRID_COMP.o

DEPEND_FILES = ${MODULES:.o=.F90}


all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $?
	
$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

clean:
	cd genocn        && gmake clean
	cd hycom         && gmake clean
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst depend

include ../conf/make.rules

include depend
