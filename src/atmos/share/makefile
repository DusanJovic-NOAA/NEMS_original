include ../../conf/configure.nems

MAKEFILE = makefile

LIBRARY  = libshare.a

MODULES  = module_ERR_MSG.o \
           module_INCLUDE.o \
           module_CONSTANTS.o \
           module_NEMSIO_MPI.o

MODULES_F90 = ${MODULES:.o=.F90}
DEPEND_FILES = $(MODULES_F90)

all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES)
	$(AR) $(ARFLAGS) $@ $(MODULES)

$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst lm map depend

include ../../conf/make.rules

include depend
