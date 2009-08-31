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
include ../../ush/compile_settings
MAKEFILE = makefile
TARGET= nmm_phys.a
#

MODULES=	\
		module_GET_CONFIG_PHY.o \
		module_H_TO_V.o \
		module_SF_URBAN.o \
		module_surface.o \
		module_LANDSURFACE.o \
		module_PHYSICS_FIELDS.o \
		module_PRECIP_ADJUST.o \
		module_PHYSICS_INIT_READ.o \
		module_PHYSICS_GRID_COMP.o \
		module_PHYSICS_OUTPUT.o \
		module_GWD.o 

MODULE_SMP=     \
		module_FLTBNDS.o

MODULE_SMP1=    \
		module_CONVECTION.o \
		module_MICROPHYSICS.o \
		module_PHYSICS_INTERNAL_STATE.o \
		module_RADIATION.o \
		module_SURFACE_LAYER.o \
		module_TURBULENCE.o

.SUFFIXES:	.F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)

$(LIBDIR)/lib/$(TARGET): $(MODULES) $(MODULE_SMP) $(MODULE_SMP1)
	$(AR) $(ARFLAGS) $@ $?

$(MODULES): %.o: %.f90 $(DEPS)
	$(FC) $(FFLAG_NMM_PHY) -c $*.f90

$(MODULE_SMP): %.o: %.f90 $(DEPS)
	$(FC) $(FFLAG_NMM_PHY) $(SMP) $(MAF) -c $*.f90

$(MODULE_SMP1): %.o: %.f90 $(DEPS)
	$(FC) $(FFLAG_NMM_PHY) $(SMP) -c $*.f90

clean:	
	/bin/rm -f  $(LIBDIR)/lib/$(TARGET) 
	/bin/rm -f  *.f90 
	/bin/rm -f  *.lst 
	/bin/rm -f  *.o 
	/bin/rm -f  *.mod 
	/bin/rm -f  lm 
	/bin/rm -f  map

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
				module_surface.o \
			   	module_LANDSURFACE.o \
			   	module_MICROPHYSICS.o \
			   	module_RADIATION.o \
			   	module_SURFACE_LAYER.o \
			   	module_TURBULENCE.o \
				module_PRECIP_ADJUST.o \
				module_GWD.o \
				module_PHYSICS_INIT_READ.o 

module_GET_CONFIG_PHY.o:        \
				module_PHYSICS_INTERNAL_STATE.o 

module_PHYSICS_FIELDS.o:	\
				module_PHYSICS_INTERNAL_STATE.o

module_PHYSICS_OUTPUT.o:	\
				module_PHYSICS_INTERNAL_STATE.o 

module_PHYSICS_INTERNAL_STATE.o:	\
					module_SF_URBAN.o \
					module_surface.o \
					module_LANDSURFACE.o \
					module_MICROPHYSICS.o 

module_TURBULENCE.o:	\
		    	module_H_TO_V.o \
			module_SF_URBAN.o \
			module_surface.o \
		    	module_LANDSURFACE.o \
			module_SURFACE_LAYER.o \
		    	module_GWD.o 

module_RADIATION.o:	\
		    	module_MICROPHYSICS.o

module_SURFACE_LAYER.o:	\
			module_SF_URBAN.o \
			module_surface.o \
		    	module_LANDSURFACE.o 

module_LAND_SURFACE.o:	\
			module_surface.o \
			module_SF_URBAN.o 

module_PHYSICS_INIT_READ.o: \
			module_MICROPHYSICS.o
