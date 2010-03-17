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
MAKEFILE = makefile
TARGET= nmm_phys.a
#

MODULES=	\
		module_GET_CONFIG_PHY.o \
		module_H_TO_V.o \
		module_SF_URBAN.o \
		module_LS_NOAHLSM.o \
		module_LS_LISS.o \
		module_PRECIP_ADJUST.o \
		module_PHYSICS_INIT_READ_BIN.o module_PHYSICS_INIT_READ_NEMSIO.o \
		module_PHYSICS_GRID_COMP.o \
		module_PHYSICS_OUTPUT.o \
		module_GWD.o

MODULE_SMP=    \
		module_CU_BMJ.o module_CONVECTION.o \
		module_CU_BMJ_DEV.o module_CONVECTION_DEV.o \
		module_MP_ETANEW.o module_MP_WSM6.o module_MICROPHYSICS.o \
		module_PHYSICS_INTERNAL_STATE.o \
		module_RADIATION.o \
		module_SF_JSFC.o module_BL_MYJPBL.o module_TURBULENCE.o

MODULE_RRTM=    \
		module_RADIATION_RRTM.o

MODULE_RRTM_77= \
		machine.o \
		resol_def.o

MODULE_RRTM_90= \
		radiation_astronomy.o \
		radiation_gases.o \
		radiation_surface.o \
		module_bfmicrophysics.o \
		radiation_clouds.o \
		radiation_aerosols.o \
		radlw_param.o \
		radlw_main.o \
		radlw_datatb.o \
		grrad.o \
		radsw_param.o \
		radsw_main.o \
		radsw_datatb.o \
		physcons_v.o \
		funcphys_v.o \
		iounitdef.o

MODULES_f90 = ${MODULES:.o=.f90}
MODULE_SMP_f90 = ${MODULE_SMP:.o=.f90}
MODULE_RRTM_f90 = ${MODULE_RRTM:.o=.f90}
MODULE_RRTM_77_f77 = ${MODULE_RRTM_77:.o=.f}
MODULE_RRTM_90_f90=${MODULE_RRTM_90:.o=.f}

.SUFFIXES:	.F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)

$(LIBDIR)/lib/$(TARGET): $(MODULES) $(MODULE_SMP) $(MODULE_RRTM_77) $(MODULE_RRTM_90) $(MODULE_RRTM)
	$(AR) $(ARFLAGS) $@ $?

$(MODULES): %.o: %.f90 $(DEPS)
	$(F90) $(FFLAG_NMM_PHY) -c $*.f90

$(MODULE_SMP): %.o: %.f90 $(DEPS)
	$(F90) $(FFLAG_NMM_PHY) $(SMP) -c $*.f90

$(MODULE_RRTM): %.o: %.f90 $(DEPS)
	$(F90) $(OPT_NMM) -I$(ESMF_MOD) $(TRAPS) -I$(NMM_UTIL) -I$(NMM_IO) $(SMP) -c $*.f90

$(MODULE_RRTM_77): %.o: %.f $(DEPS)
	$(F77) $(OPT_NMM) -I$(ESMF_MOD) $(TRAPS) -qxlf77=leadzero -qrealsize=4 -qnolm -qsmp=noauto -qnosave -c $*.f

$(MODULE_RRTM_90): %.o: %.f $(DEPS)
	$(F90) $(OPT_NMM) -I$(ESMF_MOD) $(TRAPS) -qrealsize=8 -c $*.f

depend: $(MODULES_f90) $(MODULE_SMP_f90) $(MODULE_RRTM_f90) $(MODULE_RRTM_77_f77) $(MODULE_RRTM_90_f90)
	$(MAKEDEPEXEC)

clean:	
	/bin/rm -f  $(LIBDIR)/lib/$(TARGET) 
	/bin/rm -f  *.f90 
	/bin/rm -f  *.lst 
	/bin/rm -f  *.o 
	/bin/rm -f  *.mod 
	/bin/rm -f  lm 
	/bin/rm -f  map
	/bin/rm -f  depend

#
# DEPENDENCIES :
#
include depend 
