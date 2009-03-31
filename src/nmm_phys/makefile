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
EXEC     = nmm_phys.a
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
#CPPFLAGS = -DIBM
CPPFLAGS = -DENABLE_SMP -DIBM
SMP      = -qsmp=omp
#SMP      =
MAF      = -qfloat=nomaf
OPTS     = -O2 -qmaxmem=-1 -qarch=auto -qtune=auto $(PROF) -NS2000
#TRAPS    = -qinitauto=FF911299 -qflttrap=ov:zero:inv:en -qsigtrap -C
TRAPS    = 
MAP      = -bmap:map -bloadmap:lm
#
# Assemble Options
#
FINCO    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
FINCS	 = $(FINCO) $(GCINCS)
FFLAGS   = $(OPTS) $(FINCS) $(TRAPS)
#

MODULES=	\
		module_GET_CONFIG_PHY.o \
		module_H_TO_V.o \
		module_SF_URBAN.o \
		module_LANDSURFACE.o \
		module_PHYSICS_FIELDS.o \
		module_PHYSICS_GRID_COMP.o \
		module_PHYSICS_OUTPUT.o \
		module_GWD.o \

MODULE_SMP=     \
                module_FLTBNDS.o

MODULE_SMP1=    module_CONVECTION.o \
                module_MICROPHYSICS.o \
                module_PHYSICS_INTERNAL_STATE.o \
                module_RADIATION.o \
		module_SURFACE_LAYER.o \
		module_TURBULENCE.o

INCLUDES=	kind.inc


DEPS= $(INCLUDES)

AR = ar
ARFLAGS = -r


.SUFFIXES:	.F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

$(EXEC): $(MODULES) $(MODULE_SMP) $(MODULE_SMP1)
	$(AR) $(ARFLAGS) $(EXEC) $(MODULES) $(MODULE_SMP) $(MODULE_SMP1)

$(MODULES):	$(DEPS)
	$(FC) $(FFLAGS) -c $<

$(MODULE_SMP):	$(DEPS)
	$(FC) $(FFLAGS) $(SMP) $(MAF) -qsource -qlist -c $<

$(MODULE_SMP1):	$(DEPS)
	$(FC) $(FFLAGS) $(SMP) -qsource -qlist -c $<

clean:	
	/bin/rm -f  $(TARGET) *.f90 *.lst *.o *.mod *.a lm map

#
# DEPENDENCIES :
#
 
module_PHYSICS_GRID_COMP.o:	\
			   	module_FLTBNDS.o \
			   	module_GET_CONFIG_PHY.o \
			   	module_PHYSICS_INTERNAL_STATE.o \
			   	module_PHYSICS_FIELDS.o \
			   	module_PHYSICS_OUTPUT.o \
			   	module_H_TO_V.o \
			   	module_CONVECTION.o \
				module_SF_URBAN.o \
			   	module_LANDSURFACE.o \
			   	module_MICROPHYSICS.o \
			   	module_RADIATION.o \
			   	module_SURFACE_LAYER.o \
			   	module_TURBULENCE.o \
				module_GWD.o 

module_GET_CONFIG_PHY.o:        module_PHYSICS_INTERNAL_STATE.o 

module_PHYSICS_FIELDS.o:	\
				module_PHYSICS_INTERNAL_STATE.o

module_PHYSICS_OUTPUT.o:	module_PHYSICS_INTERNAL_STATE.o 

module_PHYSICS_INTERNAL_STATE.o:	\
					module_SF_URBAN.o \
					module_LANDSURFACE.o \
					module_MICROPHYSICS.o 

module_TURBULENCE.o:	\
		    	module_H_TO_V.o \
			module_SF_URBAN.o \
		    	module_LANDSURFACE.o \
			module_SURFACE_LAYER.o \
		    	module_GWD.o 

module_RADIATION.o:	\
		    	module_MICROPHYSICS.o

module_SURFACE_LAYER.o:	\
			module_SF_URBAN.o \
		    	module_LANDSURFACE.o 

module_LAND_SURFACE.o:	\
			module_SF_URBAN.o 
