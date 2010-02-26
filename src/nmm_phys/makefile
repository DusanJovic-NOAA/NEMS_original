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
		module_LS_NOAHLSM.o \
		module_LS_LISS.o \
		module_PRECIP_ADJUST.o \
		module_PHYSICS_INIT_READ_BIN.o module_PHYSICS_INIT_READ_NEMSIO.o \
		module_PHYSICS_GRID_COMP.o \
		module_PHYSICS_OUTPUT.o \
		module_GWD.o 

MODULE_SMP=     \
		module_FLTBNDS.o

MODULE_SMP1=    \
		module_CU_BMJ.o module_CONVECTION.o \
		module_CU_BMJ_DEV.o module_CONVECTION_DEV.o \
		module_MP_ETANEW.o module_MICROPHYSICS.o \
		module_PHYSICS_INTERNAL_STATE.o \
		module_RADIATION.o \
		module_SF_JSFC.o module_BL_MYJPBL.o module_TURBULENCE.o

MODULE_RRTM=    \
		module_RADIATION_RRTM.o

MODULE_RRTM_77=	\
		machine.o \
		resol_def.o

MODULE_RRTM_90=	\
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

.SUFFIXES:	.F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)

$(LIBDIR)/lib/$(TARGET): $(MODULES) $(MODULE_SMP) $(MODULE_SMP1) $(MODULE_RRTM_77) $(MODULE_RRTM_90) $(MODULE_RRTM)
	$(AR) $(ARFLAGS) $@ $?

$(MODULES): %.o: %.f90 $(DEPS)
	$(FC) $(FFLAG_NMM_PHY) -c $*.f90

$(MODULE_SMP): %.o: %.f90 $(DEPS)
	$(FC) $(FFLAG_NMM_PHY) $(SMP) $(MAF) -c $*.f90

$(MODULE_SMP1): %.o: %.f90 $(DEPS)
	$(FC) $(FFLAG_NMM_PHY) $(SMP) -c $*.f90

$(MODULE_RRTM): %.o: %.f90 $(DEPS)
	$(FC) $(OPT_NMM) -I$(ESMF_MOD) $(TRAPS) -I$(NMM_UTIL) -I$(NMM_IO) $(SMP) -c $*.f90

$(MODULE_RRTM_77): %.o: %.f $(DEPS)
	$(F77) $(OPT_NMM) -I$(ESMF_MOD) $(TRAPS) -qxlf77=leadzero -qrealsize=4 -qnolm -qsmp=noauto -qnosave -c $*.f

$(MODULE_RRTM_90): %.o: %.f $(DEPS)
	$(F90) $(OPT_NMM) -I$(ESMF_MOD) $(TRAPS) -qrealsize=8 -c $*.f

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
			   	module_PHYSICS_OUTPUT.o \
			   	module_H_TO_V.o \
			   	module_CONVECTION.o \
			   	module_CONVECTION_DEV.o \
				module_SF_URBAN.o \
			   	module_LS_NOAHLSM.o \
			   	module_LS_LISS.o \
			   	module_MICROPHYSICS.o \
			   	module_RADIATION.o \
			   	module_RADIATION_RRTM.o \
			   	module_TURBULENCE.o \
				module_PRECIP_ADJUST.o \
				module_GWD.o \
				module_PHYSICS_INIT_READ_BIN.o \
				module_PHYSICS_INIT_READ_NEMSIO.o \
                                        module_MP_ETANEW.o

module_GET_CONFIG_PHY.o:        \
				module_PHYSICS_INTERNAL_STATE.o 

module_PHYSICS_FIELDS.o:	\
				module_PHYSICS_INTERNAL_STATE.o

module_PHYSICS_OUTPUT.o:	\
				module_PHYSICS_INTERNAL_STATE.o 

module_PHYSICS_INTERNAL_STATE.o:	\
					module_SF_URBAN.o \
					module_LS_NOAHLSM.o \
					module_LS_LISS.o \
					module_MICROPHYSICS.o 

module_CONVECTION.o:	\
		    	module_CU_BMJ.o

module_CONVECTION_DEV.o:	\
		    	module_CU_BMJ_DEV.o

module_TURBULENCE.o:	\
		    	module_H_TO_V.o \
			module_SF_URBAN.o \
                        module_BL_MYJPBL.o \
                        module_SF_JSFC.o \
		    	module_LS_NOAHLSM.o \
		    	module_LS_LISS.o \
		    	module_GWD.o 

module_RADIATION.o:	\
		    	module_MICROPHYSICS.o \
                              module_MP_ETANEW.o

module_RADIATION_RRTM.o:	\
		    	module_MICROPHYSICS.o \
			module_MP_ETANEW.o \
		    	grrad.o \
		    	radiation_astronomy.o

module_PHYSICS_INIT_READ.o: \
			module_MICROPHYSICS.o

module_MICROPHYSICS.o: \
			module_MP_ETANEW.o

grrad.o:		\
			machine.o \
			physcons_v.o \
			funcphys_v.o \
			resol_def.o \
			radiation_astronomy.o \
			radiation_gases.o \
			radiation_aerosols.o \
			radiation_surface.o \
			radiation_clouds.o \
			radsw_param.o \
			radsw_main.o \
			radlw_param.o \
			radlw_main.o

physcons_v.o:		\
			machine.o

funcphys_v.o:		\
			machine.o \
			physcons_v.o

radiation_astronomy.o:	\
			machine.o \
			physcons_v.o \
			iounitdef.o

radiation_gases.o:	\
			machine.o \
			physcons_v.o \
			iounitdef.o \
			resol_def.o

radiation_aerosols.o:	\
			machine.o \
			physcons_v.o \
			iounitdef.o \
			radsw_param.o \
			radlw_param.o

radiation_surface.o:	\
			machine.o \
			physcons_v.o \
			iounitdef.o

radiation_clouds.o:	\
			machine.o \
			physcons_v.o \
			iounitdef.o \
			module_bfmicrophysics.o

module_bfmicrophysics.o:	\
			machine.o \
			physcons_v.o \
			funcphys_v.o

radsw_param.o:		\
			machine.o

radsw_main.o:		\
			machine.o \
			physcons_v.o \
			radsw_param.o \
			radsw_datatb.o

radsw_datatb.o:		\
			machine.o \
			radsw_param.o

radlw_param.o:		\
			machine.o

radlw_main.o:		\
			machine.o \
			physcons_v.o \
			radlw_param.o \
			radlw_datatb.o

radlw_datatb.o:		\
			machine.o \
			radlw_param.o
