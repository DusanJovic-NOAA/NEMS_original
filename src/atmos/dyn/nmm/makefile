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
#         TARGET   - name of the executable
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
# CPP, Compiler, and Linker Options
#
include ../../../../build/compile_settings
MAKEFILE=makefile
TARGET= nmm_dyn.a
#

MODULES=	\
		module_DYNAMICS_OUTPUT.o \
		module_DYN_PHY_CPL_DATA.o \
		module_DYN_PHY_CPL_COMP.o \
		module_GET_CONFIG_DYN.o \
		module_DYNAMICS_INIT_READ_BIN.o \
		module_DYNAMICS_INIT_READ_NEMSIO.o \
		module_DYNAMICS_GRID_COMP.o 

MODULE_SMP=     \
		module_DYNAMICS_INTERNAL_STATE.o \
		module_DYNAMICS_ROUTINES.o 


MODULES_f90 = ${MODULES:.o=.f90}
MODULES_SMP_f90 = ${MODULES_SMP:.o=.f90}

.SUFFIXES:	.F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)


$(LIBDIR)/lib/$(TARGET): $(MODULE_SMP) $(MODULES)  
	$(AR) $(ARFLAGS) $@ $?

$(MODULES): %.o: %.f90 $(DEPS)
	$(F90) $(FFLAG_NMM_DYN) -c $*.f90

$(MODULE_SMP): %.o: %.f90  $(DEPS)
	$(F90) $(FFLAG_NMM_DYN) $(SMP) $(MAF) -c $*.f90

depend: $(MODULE_SMP_f90) $(MODULES_f90)
	$(MAKEDEPEXEC) 

clean:	
	/bin/rm -f $(LIBDIR)/lib/$(TARGET) 
	/bin/rm -f *.f90 
	/bin/rm -f *.lst 
	/bin/rm -f *.o 
	/bin/rm -f *.mod 
	/bin/rm -f lm 
	/bin/rm -f map
	/bin/rm -f depend
#
# DEPENDENCIES :
#
include depend 
