include ../../conf/configure.nems

MAKEFILE = makefile

LIBRARY  = libshare.a

MODULES  = module_ERR_MSG.o \
           module_INCLUDE.o \
           module_CONSTANTS.o \
           module_NEMSIO.o \
           module_NEMSIO_MPI.o

MODULES_f90 = ${MODULES:.o=.f90}

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

default:
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): depend $(MODULES)
	$(AR) $(ARFLAGS) $@ $(MODULES)

$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

depend: $(MODULES_f90)
	$(MAKEDEPF90) -nosrc $^ > depend

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst lm map

#
# DEPENDENCIES :
#
include depend
