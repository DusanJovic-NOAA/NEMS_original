####   requires necessary grid components
 ROOT    = MAIN_NEMS
 EXEC    = $(ROOT).x
 OBJS    = $(ROOT).o
 F77     = mpxlf95_r
 LDR     = $(F77) -qsmp=noauto
#
 NMM_UTIL = nmm_util
 NMM_PHY = nmm_phys
 NMM_DYN = nmm_dyn
 ATM = atmos
 GFS_DYN = gfs_dyn
 GFS_PHY = gfs_phys
 GFS_UTIL = gfs_util
 NMM_UTIL_LIB = $(NMM_UTIL)/nmm_util.a
 NMM_PHY_LIB = $(NMM_PHY)/nmm_phys.a
 NMM_DYN_LIB = $(NMM_DYN)/nmm_dyn.a
 ATM_LIB = $(ATM)/atmos.a
 GFS_DYN_LIB = $(GFS_DYN)/gfs_dyn.a
 GFS_PHY_LIB = $(GFS_PHY)/gfs_phys.a
 GFS_UTIL_LIB = $(GFS_UTIL)/gfs_util.a
 GCINCS  = -I$(NMM_UTIL) -I$(NMM_DYN) -I$(NMM_PHY) -I$(ATM) -I$(GFS_DYN) -I$(GFS_PHY) -I$(GFS_UTIL)
 GCLIBS = $(NMM_UTIL_LIB) $(NMM_DYN_LIB) $(NMM_PHY_LIB) $(ATM_LIB) $(GFS_DYN_LIB) $(GFS_PHY_LIB) $(GFS_UTIL_LIB)
 PROF     =
 SMP      = -qsmp=omp
 MAF      = -qfloat=nomaf
 TRACE    = -L /usr/lib -lmpitrace
 PROFILE  =
 MEM      =
 MAP      = -bmap:map -bloadmap:lm
 MASS     = -L/usrx/local/mass -lmass -lmassvp5
 ESSL     = -L/usr/lpp/essl.rte.rs1/lib -lessl_r
 FINC0    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
 FINCS    = $(FINC0) $(GCINCS)
 OPT_NMM   = -O0 -qmaxmem=-1 -qarch=auto -qtune=auto $(PROF) -NS2000
 OPT_GFS   = -qsuffix=cpp=f -O0 -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto -qnosave
 FFLAG_NMM   = $(OPT_NMM) $(FINCS)
 FFLAG_GFS   = $(OPT_GFS) $(FINCS) -qfree -NS2048
 ESMFLIB  = /meso/save/wx20rv/ESMF_libs/esmf_test_310r/lib/libO/AIX.default.64.mpi.default
 ESMFLIBS = ${ESMFLIB}/libesmf.a 
 LDFLAG_GFS = -lessl_r -lmass -qsmp=noauto ${ESMFLIBS}
 LDFLAG_NMM  = $(MEM) $(MAP) $(SMP) $(PROF) $(TRACE)
 LIB_COM    = ${GCLIBS} -L /nwprod/lib/ -l w3_d -l bacio_4 -lsp_d
 LIB_GFS    = -lC $(GCLIBS) $(LIB_COM)
 LIB_NMM     = -lC ${ESSL} ${MASS} ${ESMFLIBS} $(LIB_COM)
 LDR     = mpxlf95_r -qsmp=noauto
 .SUFFIXES:      .F90 .f90 .o

 .F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90
#
#
# *****************************************************************
#
NMM: 
	cd $(GFS_DYN) && make -f makefile && cd ..
	cd $(GFS_PHY) && make -f makefile && cd ..
	cd $(GFS_UTIL) && make -f makefile && cd ..
	cd $(NMM_UTIL) && make -f makefile && cd ..
	cd $(NMM_DYN) && make -f makefile && cd ..
	cd $(NMM_PHY) && make -f makefile && cd ..
	cd $(ATM) && make -f makefile && cd ..
	$(F77) $(FFLAG_NMM) -c $(ROOT).F90
	$(LDR) $(LDFLAG_NMM) -o $(EXEC) $(OBJS) $(LIB_NMM)
	echo $(EXEC)" is created for NMM core."

GFS: 
	cd $(GFS_DYN) && make -f makefile && cd ..
	cd $(GFS_PHY) && make -f makefile && cd ..
	cd $(GFS_UTIL) && make -f makefile && cd ..
	cd $(NMM_UTIL) && make -f makefile && cd ..
	cd $(NMM_DYN) && make -f makefile && cd ..
	cd $(NMM_PHY) && make -f makefile && cd ..
	cd $(ATM) && make -f makefile && cd ..
	$(F77) $(FFLAG_GFS) -c $(ROOT).F90
	$(LDR) $(LDFLAG_GFS) -o $(EXEC) $(OBJS) $(LIB_GFS)
	echo $(EXEC)" is created for GFS core."

clean:
	cd $(GFS_DYN) && make -f makefile clean && cd ..
	cd $(GFS_PHY) && make -f makefile clean && cd ..
	cd $(GFS_UTIL) && make -f makefile clean && cd ..
	cd $(NMM_UTIL) && make -f makefile clean && cd ..
	cd $(NMM_DYN) && make -f makefile clean && cd ..
	cd $(NMM_PHY) && make -f makefile clean && cd ..
	cd $(ATM) && make -f makefile clean && cd ..
	rm -f *.o *.mod *.f *.lst $(EXEC)
	echo "directories cleaned."
