#!/usr/bin/ksh

COMPILE_TEMPLATE_DYN=$1
COMPILE_TEMPLATE_PHY=$2
COMPILE_TEMPLATE_CHEM=$3
OS=$4
cp build_setup/compile_settings_${OS} build_setup/scratch/compile_settings_orig
case $COMPILE_TEMPLATE_DYN in
        "1")
 	if [[ $COMPILE_TEMPLATE_PHY -eq 1 ]] ; then
        	cat build_setup/scratch/compile_settings_orig | sed s:_NMMUTILLIB_:'$(LIBDIR)/lib/nmm_libutil.a':g \
        	| sed s:_NMMDYNLIB_:'$(LIBDIR)/lib/nmm_dyn.a':g \
        	| sed s:_NMMIOLIB_:'$(LIBDIR)/lib/nmm_io.a':g | sed s:_GFSUTILLIB_:'$(LIBDIR)/libstub/gfs_libutil_stub.a':g \
        	| sed s:_GFSDYNLIB_:'$(LIBDIR)/libstub/gfs_dyn_stub.a':g \
        	| sed s:_GFSIOLIB_:'$(LIBDIR)/libstub/gfs_io_stub.a':g | sed s:_GEFSCPL_:'$(LIBDIR)/libstub/GEFS_Cpl_stub.a':g \
		| sed s:_ATMOSSTUBLIB_:'$(LIBDIR)/libstub/atmos_gfs_stub.a':g \
        	| sed s:_W3LIB_:'$(W3_LIB)':g | sed s:_OPTS1_:'$(OPT_NMM)':g | sed s:_FINCS1_:'$(FINCS_NMM)':g | sed s:_OPTS2_:'$(OPT_GFS)':g \
        	| sed s:_FINCS2_:'$(FINCS_GFS)':g | sed s:_NEMSEXE_:'NEMS.x':g \
		> build_setup/scratch/compile_settings_dyn
	else
                cat build_setup/scratch/compile_settings_orig | sed s:_NMMUTILLIB_:'$(LIBDIR)/lib/nmm_libutil.a':g \
                | sed s:_NMMDYNLIB_:'$(LIBDIR)/lib/nmm_dyn.a':g \
                | sed s:_NMMIOLIB_:'$(LIBDIR)/lib/nmm_io.a':g | sed s:_GFSUTILLIB_:'$(LIBDIR)/lib/gfs_libutil.a':g \
                | sed s:_GFSDYNLIB_:'$(LIBDIR)/libstub/gfs_dyn_stub.a':g \
                | sed s:_GFSIOLIB_:'$(LIBDIR)/libstub/gfs_io_stub.a':g | sed s:_GEFSCPL_:'$(LIBDIR)/libstub/GEFS_Cpl_stub.a':g \
		| sed s:_ATMOSSTUBLIB_:'$(LIBDIR)/libstub/atmos_gfs_stub.a':g \
                | sed s:_W3LIB_:'$(W3_LIB)':g | sed s:_OPTS1_:'$(OPT_NMM)':g | sed s:_FINCS1_:'$(FINCS_NMM)':g | sed s:_OPTS2_:'$(OPT_GFS)':g \
                | sed s:_FINCS2_:'$(FINCS_GFS)':g | sed s:_NEMSEXE_:'NEMS.x':g \
                > build_setup/scratch/compile_settings_dyn
        fi
	;;
	"2")
        cat build_setup/scratch/compile_settings_orig | sed s:_NMMUTILLIB_:'$(LIBDIR)/libstub/nmm_libutil_stub.a':g \
        | sed s:_NMMDYNLIB_:'$(LIBDIR)/libstub/nmm_dyn_stub.a':g \
        | sed s:_NMMIOLIB_:'$(LIBDIR)/libstub/nmm_io_stub.a':g | sed s:_GFSUTILLIB_:'$(LIBDIR)/lib/gfs_libutil.a':g \
        | sed s:_GFSDYNLIB_:'$(LIBDIR)/lib/gfs_dyn.a':g \
        | sed s:_GFSIOLIB_:'$(LIBDIR)/lib/gfs_io.a':g | sed s:_GEFSCPL_:'$(LIBDIR)/lib/GEFS_Cpl.a':g \
	| sed s:_ATMOSSTUBLIB_:'$(LIBDIR)/libstub/atmos_nmm_stub.a':g \
        | sed s:_W3LIB_:' ':g | sed s:_OPTS1_:'$(OPT_GFS)':g | sed s:_FINCS1_:'$(FINCS_GFS)':g | sed s:_OPTS2_:'$(OPT_NMM)':g \
        | sed s:_FINCS2_:'$(FINCS_NMM)':g | sed s:_NEMSEXE_:'NEMS.x':g \
	>  build_setup/scratch/compile_settings_dyn
	;;
	"3")
        cat build_setup/scratch/compile_settings_orig | sed s:_NMMUTILLIB_:'$(LIBDIR)/lib/nmm_libutil.a':g \
        | sed s:_NMMDYNLIB_:'$(LIBDIR)/lib/nmm_dyn.a':g \
        | sed s:_NMMIOLIB_:'$(LIBDIR)/lib/nmm_io.a':g | sed s:_GFSUTILLIB_:'$(LIBDIR)/lib/gfs_libutil.a':g \
        | sed s:_GFSDYNLIB_:'$(LIBDIR)/lib/gfs_dyn.a':g \
        | sed s:_GFSIOLIB_:'$(LIBDIR)/lib/gfs_io.a':g | sed s:_GEFSCPL_:'$(LIBDIR)/lib/GEFS_Cpl.a':g \
	| sed s:_ATMOSSTUBLIB_:'$(LIBDIR)/lib/atmos.a':g \
        | sed s:_W3LIB_:'$(W3_LIB)':g | sed s:_OPTS1_:'$(OPT_NMM)':g | sed s:_FINCS1_:'$(FINCS)':g | sed s:_OPTS2_:'':g \
        | sed s:_FINCS2_:'':g | sed s:_NEMSEXE_:'NEMS.x':g \
	> build_setup/scratch/compile_settings_dyn
	;;
esac
case $COMPILE_TEMPLATE_PHY in
        "1")
        cat build_setup/scratch/compile_settings_dyn | sed s:_NMMPHYLIB_:'$(LIBDIR)/lib/nmm_phys.a':g \
        | sed s:_GFSPHYLIB_:'$(LIBDIR)/libstub/gfs_phys_stub.a':g \
	> build_setup/scratch/compile_settings_phys
        ;;
        "2")
	if [[ $COMPILE_TEMPLATE_DYN -eq 1 ]] ; then
        	cat build_setup/scratch/compile_settings_dyn | sed s:_NMMPHYLIB_:'$(LIBDIR)/lib/nmm_phys.a':g \
        	| sed s:_GFSPHYLIB_:'$(LIBDIR)/lib/gfs_phys.a':g \
        	> build_setup/scratch/compile_settings_phys
	else
        	cat build_setup/scratch/compile_settings_dyn | sed s:_NMMPHYLIB_:'$(LIBDIR)/libstub/nmm_phys_stub.a':g \
        	| sed s:_GFSPHYLIB_:'$(LIBDIR)/lib/gfs_phys.a':g \
		> build_setup/scratch/compile_settings_phys
	fi
        ;;
        "3")
        cat build_setup/scratch/compile_settings_dyn | sed s:_NMMPHYLIB_:'$(LIBDIR)/lib/nmm_phys.a':g \
        | sed s:_GFSPHYLIB_:'$(LIBDIR)/lib/gfs_phys.a':g \
	> build_setup/scratch/compile_settings_phys
        ;;
esac
case $COMPILE_TEMPLATE_CHEM in
	"0")
	cat build_setup/scratch/compile_settings_phys | sed s:_NEMSTOPDIR_:"${NEMSTOPDIR}":g \
	| sed s:_GOCARTMODE_:'stub':g > build_setup/scratch/compile_settings_chem
	;;
	"1")
        cat build_setup/scratch/compile_settings_phys | sed s:_GOCARTMODE_:'full':g > build_setup/scratch/compile_settings_chem
	cat build_setup/scratch/compile_settings_phys | sed s:_NEMSTOPDIR_:"${NEMSTOPDIR}":g \
	| sed s:_GOCARTMODE_:'full':g > build_setup/scratch/compile_settings_chem
	;;
esac

cp build_setup/scratch/compile_settings_chem compile_settings
rm build_setup/scratch/compile_settings_orig
rm build_setup/scratch/compile_settings_dyn
rm build_setup/scratch/compile_settings_phys
rm build_setup/scratch/compile_settings_chem
