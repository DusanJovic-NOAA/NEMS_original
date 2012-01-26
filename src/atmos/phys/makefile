include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../share
GOCARTINCS = -I../../chem/gocart/AIX/include/Chem_Base

LIBRARY  = libphys.a

MODULES_NMM      = module_LS_LISS.o \
                   module_LS_NOAHLSM.o \
                   module_SF_URBAN.o

MODULES_NMM_SMP  = module_BL_MYJPBL.o \
                   module_BL_GFSPBL.o \
                   module_CU_SAS.o \
                   module_CU_BMJ.o \
                   module_CU_BMJ_DEV.o \
                   module_MP_ETANEW.o \
                   module_MP_FER_HIRES.o \
                   module_MP_WSM6.o \
                   module_MP_GFS.o \
                   module_RA_GFDL.o \
                   module_RA_RRTM.o \
                   module_SF_JSFC.o

MODULES_NMM_RRTM = machine.o \
                   grrad_nmmb.o              \
                   radiation_astronomy_nmmb.o\
                   radiation_aerosols.o \
                   ozphys.o \
                   ozinterp.o \
                   radiation_gases.o \
                   radiation_surface.o \
                   module_bfmicrophysics.o \
                   radiation_clouds.o \
                   radlw_param.o \
                   radlw_main.o \
                   radlw_datatb.o \
                   radsw_param.o \
                   radsw_main.o \
                   radsw_datatb.o \
                   physcons.o \
                   funcphys.o \
                   iounitdef.o

MODULES_GFS =                                \
	           machine.o                 \
	           num_parthds.o             \
	           physcons.o                \
	           funcphys.o                \
	           module_nst_parameters.o   \
	           module_nst_water_prop.o   \
	           module_nst_model.o        \
	           calpreciptype.o

MODULES_GFSF=                                \
	           iounitdef.o               \
	           ozne_def.o                \
	           progtm_module.o           \
	           date_def.o                \
	           rascnvv2.o                \
	           namelist_soilveg.o        \
	           module_bfmicrophysics.o   \
	           mersenne_twister.o        \
	           tracer_const_h.o          \
	           gfs_phy_tracer_config.o   \
	           radlw_param.o             \
 	           radlw_datatb.o            \
	           radlw_main.o              \
	           radsw_param.o             \
	           radsw_datatb.o            \
	           radsw_main.o              \
	           radiation_astronomy.o     \
	           radiation_aerosols.o      \
	           radiation_gases.o         \
	           radiation_clouds.o        \
	           radiation_surface.o       \
	           grrad.o                   \
	           ozphys.o                  \
	           ozinterp.o                \
	           sfcsub.o                  \
	           gbphys.o                  \
	           dcyc2.o                   \
	           dcyc2.pre.rad.o           \
	           set_soilveg.o             \
	           sfc_drv.o                 \
	           sfc_land.o                \
	           progt2.o                  \
	           sfc_sice.o                \
	           sfc_ocean.o               \
	           sfc_nst.o                 \
	           sfc_diff.o                \
	           sfc_diag.o                \
	           sflx.o                    \
	           moninp.o                  \
	           moninp1.o                 \
	           moninq.o                  \
	           gwdps.o                   \
	           gwdc.o                    \
	           sascnv.o                  \
	           sascnvn.o                 \
	           cnvc90.o                  \
		   tridi2t3.o                \
	           shalcv.o                  \
	           shalcv_opr.o              \
	           shalcnv.o                 \
	           lrgsclr.o                 \
	           gscond.o                  \
	           precpd.o                  \
	           mstadb.o                  \
	           mstadbtn.o                \
	           mstcnv.o                  \
	           get_prs.o                 \
	           gsmddrive.o

MODULES_TRAC=      gocart_tracer_config.o

MODULES_TRAC_stub= gocart_tracer_config_stub.o

MODULES_ALL = $(MODULES_NMM) $(MODULES_GFS) $(MODULES_TRAC_stub) $(MODULES_GFSF) $(MODULES_NMM_SMP) $(MODULES_NMM_RRTM)
MODULES_GOCART = $(MODULES_NMM) $(MODULES_GFS) $(MODULES_TRAC) $(MODULES_GFSF) $(MODULES_NMM_SMP) $(MODULES_NMM_RRTM)

ifeq ($(GOCART_MODE),full)
	DEPEND_FILES = ${MODULES_NMM:.o=.F90} ${MODULES_GFS:.o=.f} ${MODULES_TRAC:.o=.f} ${MODULES_NMM_SMP:.o=.F90} ${MODULES_NMM_RRTM:.o=.f} ${MODULES_GFSF:.o=.f}
else
	DEPEND_FILES = ${MODULES_NMM:.o=.F90} ${MODULES_GFS:.o=.f} ${MODULES_TRAC_stub:.o=.f} ${MODULES_NMM_SMP:.o=.F90} ${MODULES_NMM_RRTM:.o=.f} ${MODULES_GFSF:.o=.f}
endif

all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULES_ALL)
	$(AR) $(ARFLAGS) $@ $(MODULES_ALL)

gocart_stub: $(MODULES_ALL)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_ALL)

gocart_full: $(MODULES_GOCART)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_GOCART)

$(MODULES_NMM): %.o: %.f90
	$(FC) $(FFLAGS_NMM) $(UTILINCS) -c $*.f90

$(MODULES_NMM_SMP): %.o: %.f90
	$(FC) $(FFLAGS_NMM) $(SMP) $(UTILINCS) -c $*.f90

$(MODULES_NMM_RRTM): %.o: %.f
	$(FC) $(FREE) $(FFLAGS_NMM) $(UTILINCS) $(R8) -c $*.f

$(MODULES_GFS): %.o: %.f
	$(FC) $(FFLAGS_GFS) $(UTILINCS) $(R8) -c $*.f

$(MODULES_TRAC): %.o: %.f
	$(FC) $(FFLAGS_GFSF) $(UTILINCS) $(GOCARTINCS) $(R8) -c $*.f

$(MODULES_TRAC_stub): %.o: %.f
	$(FC) $(FFLAGS_GFSF) $(UTILINCS) $(R8) -c $*.f

$(MODULES_GFSF): %.o: %.f
	$(FC) $(FFLAGS_GFSF) $(UTILINCS) $(R8) -c $*.f

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst lm map depend

include ../../conf/make.rules

include depend
