################################################################################
# 
#     Makefile for Global NMM Model 
#
#     Use:
#     make         -  build the executable
#     make clean   -  start with a clean slate
#
#     The following macros will be of interest:
#
#         FC       - name of Fortran compiler
#         CPP      - name of CPP
#         ARCH     - architecture
#         CPPFLAGS - CPP flags
#         OPTS     - compiler code optimizations
#         LIST     - source listing
#         SMP      - threading
#         TRAPS    - runtime traps for floating point exceptions
#         PROFILE  - source code profiling ( -pg )
#         DEBUG    - -g
#         MEM      - user data area and stack size
#         MAP      - load map
#         W3LIB    - w3lib
#         ESSL     - ESSL library
#         MASS     - MASS library
#         TRACE    - MPI tracing
#
#################################################################################
#
# Define the name of the executable
#
#
include ../../../build/compile_settings
MAKEFILE = makefile
TARGET = share.a
 
MODULES=	\
		module_NEMSIO.o \
		module_INCLUDE.o \
		module_ERR_MSG.o 

MODULES_f90 = ${MODULES:.o=.f90}

.SUFFIXES:	.F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)

$(LIBDIR)/lib/$(TARGET): $(MODULES) 
	$(AR) $(ARFLAGS) $@ $?

$(MODULES): %.o: %.f90	$(DEPS)
	$(F90) $(FFLAG_NMM_UTIL) -c $*.f90

depend: $(MODULES_f90) 
	$(MAKEDEPEXEC)

clean:	
	/bin/rm -f  $(LIBDIR)/lib/$(TARGET) 
	/bin/rm -f  *.f90 
	/bin/rm -f  *.lst 
	/bin/rm -f  *.o 
	/bin/rm -f  *.mod 
	/bin/rm -f  lm 
	/bin/rm -f  map
	/bin/rm -f depend

#
# DEPENDENCIES :
#
include depend
