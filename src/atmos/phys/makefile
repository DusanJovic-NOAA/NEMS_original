include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../share

LIBRARY  = libphys.a

MODULES_NMM      = module_LS_LISS.o \
                   module_LS_NOAHLSM.o \
                   module_SF_URBAN.o

MODULES_NMM_SMP  = module_BL_MYJPBL.o \
                   module_CU_BMJ.o \
                   module_CU_BMJ_DEV.o \
                   module_MP_ETANEW.o \
                   module_MP_WSM6.o \
                   module_RA_GFDL.o \
                   module_RA_RRTM.o \
                   module_SF_JSFC.o

MODULES_NMM_RRTM = machine.o \
                   resol_def.o \
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


MODULES_GFS =

MODULES_ALL = $(MODULES) $(MODULES_NMM) $(MODULES_NMM_SMP) $(MODULES_NMM_RRTM) $(MODULES_GFS)

DEPEND_FILES = ${MODULES_NMM:.o=.F90} ${MODULES_NMM_SMP:.o=.F90} ${MODULES_NMM_RRTM:.o=.f}


all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES_ALL)
	$(AR) $(ARFLAGS) $@ $(MODULES_ALL)

$(MODULES): %.o: %.f90
	$(FC) $(FFLAGS) $(UTILINCS) -c $*.f90

$(MODULES_NMM): %.o: %.f90
	$(FC) $(FFLAGS_NMM) $(UTILINCS) -c $*.f90

$(MODULES_NMM_SMP): %.o: %.f90
	$(FC) $(FFLAGS_NMM) $(SMP) $(UTILINCS) -c $*.f90

$(MODULES_NMM_RRTM): %.o: %.f
	$(FC) $(FREE) $(FFLAGS_NMM) $(UTILINCS) $(R8) -c $*.f

$(MODULES_GFS): %.o: %.f90
	$(FC) $(FFLAGS_GFS) $(UTILINCS) -c $*.f90

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst lm map depend

include ../../conf/make.rules

include depend
