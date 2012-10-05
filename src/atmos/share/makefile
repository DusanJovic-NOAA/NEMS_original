include ../../conf/configure.nems

MAKEFILE = makefile

LIBRARY  = libshare.a

MODULES_F90  = module_ERR_MSG.o \
               module_INCLUDE.o \
               module_CONSTANTS.o \
	       ffsync.o 

MODULES_ALL = $(MODULES_F90)

DEPEND_FILES = ${MODULES_F90:.o=.F90}


all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES_ALL)
	$(AR) $(ARFLAGS) $@ $(MODULES_ALL)

$(MODULES_F90): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst *.a lm map depend

include ../../conf/make.rules

include depend
