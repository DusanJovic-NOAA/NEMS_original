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
EXEC     = nmm_dyn.a
UTIL = nmm_util
UTILINCS = -I../$(UTIL)
GCINCS  = $(UTILINCS)
FC       = mpxlf95_r
PROF     = -g -pg
PROF     =
ARCH     = AIX
CPP      = /lib/cpp -P
CPPFLAGS = -DENABLE_SMP -DIBM
SMP      = -qsmp=omp
MAF      = -qfloat=nomaf
OPTS     = -O2 -qmaxmem=-1 -qarch=auto -qtune=auto $(PROF) -NS2000
#
# Assemble Options
#
FINCO    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
FINCS	 = $(FINCO) $(GCINCS)
FFLAGS   = $(OPTS) $(FINCS)  
#

MODULES=	\
		module_DYNAMICS_FIELDS.o \
		module_DYNAMICS_OUTPUT.o \
		module_DYN_PHY_CPL_COMP.o \
		module_DYN_PHY_CPL_DATA.o \
		module_GET_CONFIG_DYN.o \

MODULE_DGC=	module_DYNAMICS_GRID_COMP.o

MODULE_SMP=     module_DYNAMICS_INTERNAL_STATE.o \
                module_DYNAMICS_ROUTINES.o \
                module_FLTBNDS.o

INCLUDES=	kind.inc


DEPS= $(INCLUDES)

AR = ar
ARFLAGS = -r


.SUFFIXES:	.F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

$(EXEC): $(MODULES) $(MODULE_DGC) $(MODULE_SMP) $(MODULE_SMP1)
	$(AR) $(ARFLAGS) $(EXEC) $(MODULES) $(MODULE_DGC) $(MODULE_SMP) $(MODULE_SMP1)

$(MODULES):	$(DEPS)
	$(FC) $(FFLAGS) -c $<

$(MODULE_DGC):	$(DEPS)
	$(FC) $(FFLAGS) -c $<

$(MODULE_SMP):	$(DEPS)
	$(FC) $(FFLAGS) $(SMP) $(MAF) -qsource -qlist -c $<

clean:	
	/bin/rm -f  $(TARGET) *.f90 *.lst *.o *.mod *.a lm map

#
# DEPENDENCIES :
#
 
module_DYN_PHY_CPL_COMP.o:	\
			  	module_DYNAMICS_GRID_COMP.o \
			  	module_DYN_PHY_CPL_DATA.o 

module_DYNAMICS_GRID_COMP.o:	\
			   	module_FLTBNDS.o \
			   	module_GET_CONFIG_DYN.o \
			   	module_DYNAMICS_INTERNAL_STATE.o \
			   	module_DYNAMICS_FIELDS.o \
			   	module_DYNAMICS_OUTPUT.o \
			   	module_DYNAMICS_ROUTINES.o 

module_GET_CONFIG_DYN.o:    	module_DYNAMICS_INTERNAL_STATE.o 

module_DYNAMICS_FIELDS.o:	\
				module_DYNAMICS_INTERNAL_STATE.o


module_DYNAMICS_OUTPUT.o:	module_DYNAMICS_INTERNAL_STATE.o 

module_DYNAMICS_INTERNAL_STATE.o:	\
				 	module_FLTBNDS.o 
