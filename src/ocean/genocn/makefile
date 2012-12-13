include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../../atmos/share

LIBRARY  = libgenocn.a

MODULES  = module_GENOCN_INTERNAL_STATE.o module_GENOCN_GRID_COMP.o

MODULES_STUB  = module_GENOCN_GRID_COMP_stub.o

DEPEND_FILES = ${MODULES:.o=.F90}


all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $?
	
$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

$(MODULES_STUB): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

stub: $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_STUB)

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst depend

MKDEPENDS = ../../../exe/mkDepends.pl

include ../../conf/make.rules

include depend
