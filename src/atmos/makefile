include ../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -Ishare -Inmm -Igfs -Ifim

LIBRARY  = libatmos.a

MODULES  = module_ATM_INTERNAL_STATE.o module_ATM_GRID_COMP.o

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

default:
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $?
	
$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

clean:
	cd share        && gmake clean
	cd phys         && gmake clean
	cd nmm          && gmake clean
	cd gfs          && gmake clean
	cd fim          && gmake clean
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst
