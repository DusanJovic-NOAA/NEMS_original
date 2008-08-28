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
UTIL_LIB   = ../$(UTIL)/nmm_util.a
GCLIBS  = $(UTIL_LIB)
GCINCS  = $(UTILINCS)
FC       = mpxlf95_r
PROF     = -g -pg
PROF     =
ARCH     = AIX
CPP      = /lib/cpp -P
CPPFLAGS = -DENABLE_SMP -DIBM
#CPPFLAGS = -DENABLE_SMP
SMP      = -qsmp=omp
MAF      = -qfloat=nomaf
OPTS     = -O0 -qmaxmem=-1 -qarch=auto -qtune=auto $(PROF) -NS2000
#OPTS     = -O3 -qnostrict -qmaxmem=-1 -qarch=auto -qtune=auto $(PROF)
########     O3 failed in producing bit-identical results with diff. decomposition
TRAPS    = -NS2000 -qinitauto=FF911299 -qflttrap=ov:zero:inv:en -qsigtrap -C 
TRACE    = -L /usr/lib -lmpitrace
PROFILE  = 
DEBUG    =
MEM      =
MAP      = -bmap:map -bloadmap:lm
MASS     = -L/usrx/local/mass -lmass -lmassvp5
ESSL     = -L/usr/lpp/essl.rte.rs1/lib -lessl_r
#
# Assemble Options
#
FINCO    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
FINCS	 = $(FINCO) $(GCINCS)
FFLAGS   = $(OPTS) $(FINCS) $(LIST) $(PROFILE) $(DEBUG)
FFLAGSDB   = $(OPTS) $(FINCS) $(LIST) $(PROFILE) $(DEBUG) $(TRAPS)
FFLAGSNT = $(OPTS) $(FINCS) $(LIST) $(PROFILE) $(DEBUG)
LDFLAGS  = $(MEM) $(MAP) $(SMP) $(PROF) $(TRACE)
ESMFLIB  = -L/meso/save/wx20rv/ESMF_libs/esmf_test_310r/lib/libO/AIX.default.64.mpi.default                  \
             /meso/save/wx20rv/ESMF_libs/esmf_test_310r/lib/libO/AIX.default.64.mpi.default/libesmf.a
LIBS     = -lC ${ESSL} ${MASS} ${ESMFLIB} $(GCLIBS)
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
	$(FC) $(FFLAGSNT) -c $<

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
