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
 FINC0    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
 FINCS    = $(FINC0) $(GCINCS)
 OPT   = -O2 -qmaxmem=-1 -qarch=auto -qtune=auto -NS2000
 FFLAG  = $(OPT) $(FINCS) -qfree -NS2048
 ESMFLIB  = /meso/save/wx20rv/ESMF_libs/esmf_test_310r/lib/libO/AIX.default.64.mpi.default
 ESMFLIBS = ${ESMFLIB}/libesmf.a 
 LDFLAG = -lessl_r -lmass -qsmp=noauto ${ESMFLIBS}
 LIB_COM    = ${GCLIBS} -L /nwprod/lib/ -l w3_d -l bacio_4 -lsp_d
 LIB    = -lC $(LIB_COM)
 LDR     = mpxlf95_r -qsmp=noauto
 .SUFFIXES:      .F90 .f90 .o

 .F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90
#
#
# *****************************************************************
#
all: 
	cd $(GFS_DYN) && make -f makefile && cd ..
	cd $(GFS_PHY) && make -f makefile && cd ..
	cd $(GFS_UTIL) && make -f makefile && cd ..
	cd $(NMM_UTIL) && make -f makefile && cd ..
	cd $(NMM_DYN) && make -f makefile && cd ..
	cd $(NMM_PHY) && make -f makefile && cd ..
	cd $(ATM) && make -f makefile && cd ..
	$(F77) $(FFLAG) -c $(ROOT).F90
	$(LDR) $(LDFLAG) -o $(EXEC) $(OBJS) $(LIB)
	echo $(EXEC)" is created for NMM/GFS core."


clean:
	cd $(GFS_DYN) && make -f makefile clean && cd ..
	cd $(GFS_PHY) && make -f makefile clean && cd ..
	cd $(GFS_UTIL) && make -f makefile clean && cd ..
	cd $(NMM_UTIL) && make -f makefile clean && cd ..
	cd $(NMM_DYN) && make -f makefile clean && cd ..
	cd $(NMM_PHY) && make -f makefile clean && cd ..
	cd $(ATM) && make -f makefile clean && cd ..
	rm -f *.o *.mod *.f *.lst *.a lm map $(EXEC)
	echo "directories cleaned."
