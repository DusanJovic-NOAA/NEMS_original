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
TARGET = nmmb_esmf.x
#
# CPP, Compiler, and Linker Options
#
EXEC     = nmm_util.a
UTIL     = nmm_util
FC       = mpxlf95_r
PROF     = -g -pg
PROF     =
ARCH     = AIX
CPP      = /lib/cpp -P
CPPFLAGS = -DENABLE_SMP -DIBM
OPTS     = -O2 -qmaxmem=-1 -qarch=auto -qtune=auto $(PROF) -NS2000
#TRAPS    = -qinitauto=FF911299 -qflttrap=ov:zero:inv:en -qsigtrap -C
TRAPS    = 
#
# Assemble Options
#
FINCO    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
FINCS	 = $(FINCO) 
FFLAGS   = $(OPTS) $(FINCS) $(TRAPS)
#

PRINTMEM=       print_memory.o 
 
MODULEIO=       nemsio_module.o

MODULES=	\
		module_ATM_INTERNAL_STATE.o  \
		module_CLOCKTIMES.o \
		module_CONSTANTS.o \
		module_CONTROL.o \
		module_DIAGNOSE.o \
		module_DM_PARALLEL.o \
		module_EXCHANGE.o \
		module_INCLUDE.o \
		module_GET_CONFIG_WRITE.o \
		module_WRITE_GRID_COMP.o \
		module_WRITE_ROUTINES.o \
		module_WRITE_INTERNAL_STATE.o \
		module_DIGITAL_FILTER_NMM.o \
		module_ERR_MSG.o \
		looplimits.o

INCLUDES=	kind.inc

DEPS= $(INCLUDES)

AR = ar
ARFLAGS = -r


.SUFFIXES:	.F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

$(EXEC): $(MODULEIO) $(MODULES) 
	$(AR) $(ARFLAGS) $(EXEC) $(MODULEIO) $(MODULES) 

$(MODULEIO):     $(DEPS)
	$(FC) $(FFLAGS) -c $<

$(MODULES):	$(DEPS)
	$(FC) $(FFLAGS) -c $<

clean:	
	/bin/rm -f  $(TARGET) *.f90 *.lst *.o *.mod *.a lm map

#
# DEPENDENCIES :
#
module_WRITE_GRID_COMP.o:       module_INCLUDE.o \
                                module_CONTROL.o \
                                module_WRITE_ROUTINES.o \
                                module_WRITE_INTERNAL_STATE.o \
                                module_ERR_MSG.o

module_GET_CONFIG_WRITE.o:     	\
                                module_WRITE_INTERNAL_STATE.o \
                                module_ERR_MSG.o

module_WRITE_ROUTINES.o :       module_INCLUDE.o \
                                module_CONTROL.o \
                                module_WRITE_INTERNAL_STATE.o \
                                module_ERR_MSG.o \
				nemsio_module.o

module_CONTROL.o :      module_INCLUDE.o \
                        module_CONSTANTS.o \
                        module_DM_PARALLEL.o \
                        module_EXCHANGE.o

module_DIGITAL_FILTER_NMM.o: module_INCLUDE.o \
			     module_DM_PARALLEL.o

module_EXCHANGE.o :     module_INCLUDE.o \
                        module_DM_PARALLEL.o

module_DM_PARALLEL.o :  module_INCLUDE.o

module_CONSTANTS.o :    module_INCLUDE.o
