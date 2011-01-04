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
                   module_MP_WSM6.o \
                   module_MP_GFS.o \
                   module_RA_GFDL.o \
                   module_RA_RRTM.o \
                   module_SF_JSFC.o

MODULES_NMM_RRTM = n_machine.o \
                   n_radiation_astronomy.o \
                   n_radiation_gases.o \
                   n_radiation_surface.o \
                   n_module_bfmicrophysics.o \
                   n_radiation_clouds.o \
                   n_radiation_aerosols.o \
                   n_radlw_param.o \
                   n_radlw_main.o \
                   n_radlw_datatb.o \
                   n_grrad.o \
                   n_radsw_param.o \
                   n_radsw_main.o \
                   n_radsw_datatb.o \
                   physcons_v.o \
                   funcphys_v.o \
                   n_iounitdef.o

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

MODULES_ALL = $(MODULES) $(MODULES_NMM) $(MODULES_GFS) $(MODULES_GFSF) $(MODULES_NMM_SMP) $(MODULES_NMM_RRTM)

DEPEND_FILES = ${MODULES_NMM:.o=.F90} ${MODULES_GFS:.o=.f} ${MODULES_GFSF:.o=.f} ${MODULES_NMM_SMP:.o=.F90} ${MODULES_NMM_RRTM:.o=.f}


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

$(MODULES_GFS): %.o: %.f
	$(FC) $(FFLAGS_GFS) $(UTILINCS) $(R8) -c $*.f

$(MODULES_GFSF): %.o: %.f
	$(FC) $(FFLAGS_GFSF) $(UTILINCS) $(GOCARTINCS) $(R8) -c $*.f

clean:
	$(RM) -f $(LIBRARY) *.f90 *.o *.mod *.lst lm map depend

include ../../conf/make.rules

include depend
