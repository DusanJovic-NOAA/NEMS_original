include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../share -I../phys

LIBRARY  = libnmm.a

MODULES  = \
                module_BGRID_INTERP.o \
                module_CLOCKTIMES.o \
                module_CONTROL.o \
                module_DIAGNOSE.o \
                module_DIGITAL_FILTER_NMM.o \
                module_DM_PARALLEL.o \
                module_DOMAIN_GRID_COMP.o \
                module_DOMAIN_INTERNAL_STATE.o \
                module_DYNAMICS_GRID_COMP.o \
                module_DYNAMICS_INIT_READ_BIN.o \
                module_DYNAMICS_INIT_READ_NEMSIO.o \
                module_DYNAMICS_OUTPUT.o \
                module_DYN_PHY_CPL_COMP.o \
                module_DYN_PHY_CPL_DATA.o \
                module_EXCHANGE.o \
                module_GET_CONFIG_DYN.o \
                module_GET_CONFIG_PHY.o \
                module_GET_CONFIG_WRITE.o \
                module_GWD.o \
                module_H_TO_V.o \
                module_NESTING.o \
                module_NMM_CORE_SETUP.o \
                module_NMM_GRID_COMP.o \
                module_NMM_INTEGRATE.o \
                module_NMM_INTERNAL_STATE.o \
                module_PARENT_CHILD_CPL_COMP.o \
                module_PHYSICS_GRID_COMP.o \
                module_PHYSICS_INIT_READ_BIN.o \
                module_PHYSICS_INIT_READ_NEMSIO.o \
                module_PHYSICS_OUTPUT.o \
                module_PRECIP_ADJUST.o \
                module_TIMESERIES.o \
                module_VARS.o \
                module_VARS_STATE.o \
                module_WRITE_GRID_COMP.o \
                module_WRITE_INTERNAL_STATE.o \
                module_WRITE_ROUTINES.o


MODULE_SMP=    \
                module_CONVECTION.o \
                module_CONVECTION_DEV.o \
                module_MICROPHYSICS.o \
                module_PHYSICS_INTERNAL_STATE.o \
                module_RADIATION.o \
                module_TURBULENCE.o

MODULE_SMP_MAF= \
                module_DYNAMICS_INTERNAL_STATE.o \
                module_DYNAMICS_ROUTINES.o \
                module_FLTBNDS.o

MODULES_STUB  = module_NMM_GRID_COMP_stub.o

MODULES_f90 = ${MODULES:.o=.f90}
MODULE_SMP_f90 = ${MODULE_SMP:.o=.f90}
MODULE_SMP_MAF_f90 = ${MODULE_SMP_MAF:.o=.f90}


.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

all: depend $(MODULES) $(MODULE_SMP) $(MODULE_SMP_MAF)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES) $(MODULE_SMP) $(MODULE_SMP_MAF)
	
stub: $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_STUB)
	
$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS_NMM) $(UTILINCS) -c $*.f90

$(MODULE_SMP): %.o: %.f90
	$(FC) $(FFLAGS_NMM) $(SMP) $(UTILINCS) -c $*.f90

$(MODULE_SMP_MAF): %.o: %.f90
	$(FC) $(FFLAGS_NMM) $(SMP) $(MAF) $(UTILINCS) -c $*.f90

$(MODULES_STUB): %.o: %.f90
	$(FC) $(FFLAGS_NMM) $(UTILINCS) -c $*.f90

depend: $(MODULES_f90) $(MODULE_SMP_f90) $(MODULE_SMP_MAF_f90)
	$(MAKEDEPF90) -nosrc $^ > depend

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst lm map

#
# DEPENDENCIES :
#
include depend
