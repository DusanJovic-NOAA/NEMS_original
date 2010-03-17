#!/usr/bin/ksh

COMPILE_TEMPLATE_DYN=$1
COMPILE_TEMPLATE_PHY=$2

if [[ $COMPILE_TEMPLATE_DYN -eq 1 ]]; then
	cat ../src/atmos/drivers/makefile_stub.IN | sed s:_STUBMODULES_:'module_GFS_WRITE_STUB.o module_GFS_CORE_SETUP_STUB.o module_GOCART_ROUTINES_STUB.o module_GFS_INTEGRATE_STUB.o':g \
        | sed s:_NEMSCORE_:'gfs':g \
        > ../src/atmos/drivers/makefile_stub
        cat ../src/atmos/drivers/makefile.IN | sed s:_MODULES_:'module_NMM_CORE_SETUP.o module_NESTING.o module_NMM_INTEGRATE.o module_PARENT_CHILD_CPL_COMP.o module_ATM_GRID_COMP.o module_ATM_DRIVER_INTERNAL_STATE.o module_ATM_DRIVER_COMP.o':g \
        > ../src/atmos/drivers/makefile
	case $COMPILE_TEMPLATE_PHY in
	"1")
	cat ../src/main/Makefile_stub.IN | sed s:_MODULES_:'cd $(GFS_DYN) \&\& gmake -f makefile_stub ; cp $(GFS_DYN)/gfs_dynamics_grid_comp_mod.mod $(MODDIR)/modstub/gfs ; cd $(GFS_PHY) \&\& gmake -f makefile_stub ; cp $(GFS_PHY)/gfs_physics_grid_comp_mod.mod $(MODDIR)/modstub/gfs ; cd $(GFS_IO) \&\& gmake -f makefile_stub ; cp $(GFS_IO)/module_write_routines_gfs.mod $(MODDIR)/modstub/gfs ; cp $(GFS_IO)/module_write_grid_comp_gfs.mod $(MODDIR)/modstub/gfs ; cd $(GFS_UTIL) \&\& gmake -f makefile_stub ; cp $(GFS_UTIL)/atmos_dyn_phy_cpl_comp_mod.mod $(MODDIR)/modstub/gfs ; cp $(GFS_UTIL)/module_gfs_mpi_def.mod $(MODDIR)/modstub/gfs ; cp $(GFS_UTIL)/atmos_phy_chem_cpl_comp_mod.mod $(MODDIR)/modstub/gfs ; cp $(GFS_UTIL)/atmos_dyn_chem_cpl_comp_mod.mod $(MODDIR)/modstub/gfs ; cd $(GEFS_CPL) \&\& gmake -f makefile_stub ; cp $(GEFS_CPL)/gefs_cplcomp_esmfmod.mod $(MODDIR)/modstub/gfs ; cd $(ATM) \&\& gmake -f makefile_stub ; cp module_gfs_core_setup.mod $(MODDIR)/modstub/gfs ; cp module_gocart_integrate.mod $(MODDIR)/modstub/gfs ; cp module_gfs_integrate.mod $(MODDIR)/modstub/gfs':g \
        | sed s:_NEMSCORE_:'gfs':g \
        | sed s:_CLEANMODULES_:'cd $(GFS_DYN) \&\& gmake -f makefile_stub clean ; cd $(GFS_PHY) \&\& gmake -f makefile_stub clean ; cd $(GFS_IO) \&\& gmake -f makefile_stub clean ; cd $(GFS_UTIL) \&\& gmake -f makefile_stub clean ; cd $(GEFS_CPL) \&\& gmake -f makefile_stub clean ; cd $(ATM) \&\& gmake -f makefile_stub clean':g \
        > ../src/main/Makefile_stub
        cat ../src/main/Makefile_main.IN | sed s:_MODULES_:'cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src install GOCART_MODE=$(GOCART_MODE) ; cd $(NMM_UTIL) \&\& gmake depend \&\& gmake ; cd $(NMM_IO) \&\& gmake depend \&\& gmake ; cd $(NMM_DYN) \&\& gmake depend \&\& gmake ; cd $(NMM_PHY) \&\& gmake depend \&\& gmake ; cd $(ATM) \&\& gmake depend \&\& gmake':g \
        | sed s:_NEMSCORE_:'nmm':g \
        | sed s:_FLAGS_:'$(FFLAG_NMM_MAIN)':g \
        | sed s:_CLEANMODULES_:'cd $(NMM_UTIL) \&\& gmake clean ; cd $(NMM_IO) \&\& gmake clean ; cd $(NMM_DYN) \&\& gmake clean ; cd $(NMM_PHY) \&\& gmake clean ; cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src realclean GOCART_MODE=$(GOCART_MODE) ; cd $(ATM) \&\& gmake clean':g \
        > ../src/main/Makefile_main
	;;
	"2")
        cat ../src/main/Makefile_stub.IN | sed s:_MODULES_:'cd $(GFS_DYN) \&\& gmake -f makefile_stub ; cp $(GFS_DYN)/gfs_dynamics_grid_comp_mod.mod $(MODDIR)/modstub/gfs ; cd $(GFS_IO) \&\& gmake -f makefile_stub ; cp $(GFS_IO)/module_write_routines_gfs.mod $(MODDIR)/modstub/gfs ; cp $(GFS_IO)/module_write_grid_comp_gfs.mod $(MODDIR)/modstub/gfs ; cd $(GEFS_CPL) \&\& gmake -f makefile_stub ; cp $(GEFS_CPL)/gefs_cplcomp_esmfmod.mod $(MODDIR)/modstub/gfs ; cd $(ATM) \&\& gmake -f makefile_stub ; cp module_gfs_core_setup.mod $(MODDIR)/modstub/gfs ; cp module_gocart_integrate.mod $(MODDIR)/modstub/gfs ; cp module_gfs_integrate.mod $(MODDIR)/modstub/gfs':g \
        | sed s:_NEMSCORE_:'gfs':g \
        | sed s:_CLEANMODULES_:'cd $(GFS_DYN) \&\& gmake -f makefile_stub clean ; cd $(GFS_IO) \&\& gmake -f makefile_stub clean ; cd $(GEFS_CPL) \&\& gmake -f makefile_stub clean ; cd $(ATM) \&\& gmake -f makefile_stub clean':g \
        > ../src/main/Makefile_stub
        cat ../src/main/Makefile_main.IN | sed s:_MODULES_:'cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src install GOCART_MODE=$(GOCART_MODE) ; cd $(GFS_UTIL) \&\& gmake depend \&\& gmake ; cd $(GFS_PHY) \&\& gmake depend \&\& gmake ; cd $(NMM_UTIL) \&\& gmake depend \&\& gmake ; cd $(NMM_IO) \&\& gmake depend \&\& gmake ; cd $(NMM_DYN) \&\& gmake depend \&\& gmake ; cd $(NMM_PHY) \&\& gmake depend \&\& gmake ; cd $(ATM) \&\& gmake depend \&\& gmake':g \
        | sed s:_NEMSCORE_:'nmm':g \
        | sed s:_FLAGS_:'$(FFLAG_NMM_MAIN)':g \
        | sed s:_CLEANMODULES_:'cd $(GFS_UTIL) \&\& gmake clean ; cd $(GFS_PHY) \&\& gmake clean ; cd $(NMM_UTIL) \&\& gmake clean ; cd $(NMM_IO) \&\& gmake clean  ; cd $(NMM_DYN) \&\& gmake clean ; cd $(NMM_PHY) \&\& gmake clean ; cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src realclean GOCART_MODE=$(GOCART_MODE) ; cd $(ATM) \&\& gmake clean':g \
        > ../src/main/Makefile_main
	;;
 	esac
fi
if [[ $COMPILE_TEMPLATE_DYN -eq 2 ]]; then       
        cat ../src/atmos/drivers/makefile_stub.IN | sed s:_STUBMODULES_:'module_NMM_CORE_SETUP_STUB.o module_NESTING_STUB.o module_NMM_INTEGRATE_STUB.o module_PARENT_CHILD_CPL_COMP_STUB.o':g \
        | sed s:_NEMSCORE_:'nmm':g \
        > ../src/atmos/drivers/makefile_stub
        cat ../src/atmos/drivers/makefile.IN | sed s:_MODULES_:'module_GFS_WRITE.o module_GFS_CORE_SETUP.o module_GOCART_ROUTINES.o module_GFS_INTEGRATE.o module_ATM_GRID_COMP.o module_ATM_DRIVER_INTERNAL_STATE.o module_ATM_DRIVER_COMP.o':g \
        > ../src/atmos/drivers/makefile
	case $COMPILE_TEMPLATE_PHY in
	"1")      
        cat ../src/main/Makefile_stub.IN | sed s:_MODULES_:'cd $(NMM_UTIL) \&\& gmake -f makefile_stub ; cp $(NMM_UTIL)/module_atm_internal_state.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_clocktimes.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_control.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_diagnose.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_dm_parallel.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_include.mod $(MODDIR)/modstub/nmm ; cd $(NMM_IO) \&\& gmake -f makefile_stub ; cp $(NMM_IO)/module_write_grid_comp.mod $(MODDIR)/modstub/nmm ; cp $(NMM_IO)/module_write_routines.mod $(MODDIR)/modstub/nmm ; cp $(NMM_IO)/module_get_config_write.mod $(MODDIR)/modstub/nmm ; cd $(NMM_DYN) \&\& gmake -f makefile_stub ; cp $(NMM_DYN)/module_dynamics_grid_comp.mod $(MODDIR)/modstub/nmm ; cp $(NMM_DYN)/module_dyn_phy_cpl_comp.mod $(MODDIR)/modstub/nmm ; cp $(NMM_DYN)/module_get_config_dyn.mod $(MODDIR)/modstub/nmm ; cd $(NMM_PHY) \&\& gmake -f makefile_stub ; cp $(NMM_PHY)/module_get_config_phy.mod $(MODDIR)/modstub/nmm ; cp $(NMM_PHY)/module_physics_grid_comp.mod $(MODDIR)/modstub/nmm ; cd $(ATM) \&\& gmake -f makefile_stub ; cp module_nmm_core_setup.mod $(MODDIR)/modstub/nmm ; cp module_nesting.mod $(MODDIR)/modstub/nmm ; cp module_nmm_integrate.mod $(MODDIR)/modstub/nmm ; cp module_parent_child_cpl_comp.mod $(MODDIR)/modstub/nmm':g \
        | sed s:_NEMSCORE_:'nmm':g \
        | sed s:_CLEANMODULES_:'cd $(NMM_UTIL) \&\& gmake -f makefile_stub clean ; cd $(NMM_IO) \&\& gmake -f makefile_stub clean ; cd $(NMM_DYN) \&\& gmake -f makefile_stub clean ; cd $(NMM_PHY) \&\& gmake -f makefile_stub clean ; cd $(ATM) \&\& gmake -f makefile_stub clean':g \
        > ../src/main/Makefile_stub
	cat ../src/main/Makefile_main.IN | sed s:_MODULES_:'cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src install GOCART_MODE=$(GOCART_MODE) ; cd $(NMM_PHY) \&\& gmake depend \&\& gmake ; cd $(GFS_UTIL) \&\& gmake depend \&\& gmake ; cd $(GFS_DYN) \&\& gmake depend \&\& gmake ; cd $(GFS_PHY) \&\& gmake depend \&\& gmake ; cd $(GFS_IO) \&\& gmake depend \&\& gmake ; cd $(GEFS_CPL) \&\& gmake ; cd $(ATM) \&\& gmake depend \&\& gmake':g \
        | sed s:_NEMSCORE_:'gfs':g \
        | sed s:_FLAGS_:'$(FFLAG_GFS_MAIN)':g \
        | sed s:_CLEANMODULES_:'cd $(NMM_PHY) \&\& gmake clean ; cd $(GFS_DYN) \&\& gmake clean ; cd $(GFS_PHY) \&\& gmake clean ; cd $(GFS_IO) \&\& gmake clean ; cd $(GFS_UTIL) \&\& gmake clean ; cd $(GEFS_CPL) \&\& gmake clean ;  cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src realclean GOCART_MODE=$(GOCART_MODE) ; cd $(ATM) \&\& gmake clean':g \
        > ../src/main/Makefile_main
	;;
	"2")
	cat ../src/main/Makefile_stub.IN | sed s:_MODULES_:'cd $(NMM_UTIL) \&\& gmake -f makefile_stub ; cp $(NMM_UTIL)/module_atm_internal_state.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_clocktimes.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_control.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_diagnose.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_dm_parallel.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_include.mod $(MODDIR)/modstub/nmm ; cp $(NMM_UTIL)/module_digital_filter_nmm.mod $(MODDIR)/modstub/nmm ; cd $(NMM_IO) \&\& gmake -f makefile_stub ; cp $(NMM_IO)/module_write_grid_comp.mod $(MODDIR)/modstub/nmm ; cp $(NMM_IO)/module_write_routines.mod $(MODDIR)/modstub/nmm ; cp $(NMM_IO)/module_get_config_write.mod $(MODDIR)/modstub/nmm ; cd $(NMM_DYN) \&\& gmake -f makefile_stub ; cp $(NMM_DYN)/module_dynamics_grid_comp.mod $(MODDIR)/modstub/nmm ; cp $(NMM_DYN)/module_dyn_phy_cpl_comp.mod $(MODDIR)/modstub/nmm ; cp $(NMM_DYN)/module_get_config_dyn.mod $(MODDIR)/modstub/nmm ; cd $(NMM_PHY) \&\& gmake -f makefile_stub ; cp $(NMM_PHY)/module_get_config_phy.mod $(MODDIR)/modstub/nmm ; cp $(NMM_PHY)/module_physics_grid_comp.mod $(MODDIR)/modstub/nmm ; cd $(ATM) \&\& gmake -f makefile_stub ; cp module_nmm_core_setup.mod $(MODDIR)/modstub/nmm ; cp module_nesting.mod $(MODDIR)/modstub/nmm ; cp module_nmm_integrate.mod $(MODDIR)/modstub/nmm ; cp module_parent_child_cpl_comp.mod $(MODDIR)/modstub/nmm':g \
        | sed s:_NEMSCORE_:'nmm':g \
        | sed s:_CLEANMODULES_:'cd $(NMM_UTIL) \&\& gmake -f makefile_stub clean ; cd $(NMM_IO) \&\& gmake -f makefile_stub clean ; cd $(NMM_DYN) \&\& gmake -f makefile_stub clean ; cd $(NMM_PHY) \&\& gmake -f makefile_stub clean ; cd $(ATM) \&\& gmake -f makefile_stub clean':g \
        > ../src/main/Makefile_stub
	cat ../src/main/Makefile_main.IN | sed s:_MODULES_:'cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src install GOCART_MODE=$(GOCART_MODE) ; cd $(GFS_UTIL) \&\& gmake depend \&\& gmake ; cd $(GFS_DYN) \&\& gmake depend \&\& gmake ; cd $(GFS_PHY) \&\& gmake depend \&\& gmake ; cd $(GFS_IO) \&\& gmake depend \&\& gmake ; cd $(GEFS_CPL) \&\& gmake ; cd $(ATM) \&\& gmake depend \&\& gmake':g \
        | sed s:_NEMSCORE_:'gfs':g \
        | sed s:_FLAGS_:'$(FFLAG_GFS_MAIN)':g \
        | sed s:_CLEANMODULES_:'cd $(GFS_DYN) \&\& gmake clean ; cd $(GFS_PHY) \&\& gmake clean ; cd $(GFS_IO) \&\& gmake clean ; cd $(GFS_UTIL) \&\& gmake clean ; cd $(GEFS_CPL) \&\& gmake clean ; cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src realclean GOCART_MODE=$(GOCART_MODE) ; cd $(ATM) \&\& gmake clean':g \
        > ../src/main/Makefile_main
	;;
	esac
fi
if [[ $COMPILE_TEMPLATE_DYN -eq 3 && $COMPILE_TEMPLATE_PHY -eq 3 ]]; then
        cat ../src/atmos/drivers/makefile.IN | sed s:_MODULES_:'module_GFS_WRITE.o module_GFS_CORE_SETUP.o module_GOCART_ROUTINES.o module_GFS_INTEGRATE.o module_NMM_CORE_SETUP.o module_NESTING.o module_NMM_INTEGRATE.o module_PARENT_CHILD_CPL_COMP.o module_ATM_GRID_COMP.o  module_ATM_DRIVER_INTERNAL_STATE.o module_ATM_DRIVER_COMP.o':g \
	> ../src/atmos/drivers/makefile
        cat ../src/main/Makefile_main.IN | sed s:_MODULES_:'cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src install GOCART_MODE=$(GOCART_MODE) ; cd $(GFS_UTIL) \&\& gmake depend \&\& gmake ; cd $(GFS_DYN) \&\& gmake depend \&\& gmake ; cd $(GFS_PHY) \&\& gmake depend \&\& gmake ; cd $(GFS_IO) \&\& gmake depend \&\& gmake ; cd $(GEFS_CPL) \&\& gmake ; cd $(NMM_UTIL) \&\& gmake depend \&\& gmake ; cd $(NMM_IO) \&\& gmake depend \&\& gmake ; cd $(NMM_DYN) \&\& gmake depend \&\& gmake ; cd $(NMM_PHY) \&\& gmake depend \&\& gmake ; cd $(ATM) \&\& gmake depend \&\& gmake':g \
        | sed s:_NEMSCORE_:'nmm_gfs':g \
        | sed s:_FLAGS_:'$(FFLAG_MAIN)':g \
        | sed s:_CLEANMODULES_:'cd $(GFS_UTIL) \&\& gmake clean ; cd $(GFS_DYN) \&\& gmake clean ; cd $(GFS_PHY) \&\& gmake clean ; cd $(GFS_IO) \&\& gmake clean ; cd $(GEFS_CPL) \&\& gmake clean ; cd $(NMM_UTIL) \&\& gmake clean ; cd $(NMM_IO) \&\& gmake clean ; cd $(NMM_DYN) \&\& gmake clean ; cd $(NMM_PHY) \&\& gmake clean ; cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src realclean GOCART_MODE=$(GOCART_MODE) ; cd $(ATM) \&\& gmake clean':g \
        > ../src/main/Makefile_main
fi
if [[ $COMPILE_TEMPLATE_DYN -eq 1 && $COMPILE_TEMPLATE_PHY -eq 3 ]]; then
        cat ../src/atmos/drivers/makefile_stub.IN | sed s:_STUBMODULES_:'module_GFS_WRITE_STUB.o module_GFS_CORE_SETUP_STUB.o module_GOCART_ROUTINES_STUB.o module_GFS_INTEGRATE_STUB.o':g \
        | sed s:_NEMSCORE_:'gfs':g \
        > ../src/atmos/drivers/makefile_stub
        cat ../src/atmos/drivers/makefile.IN | sed s:_MODULES_:'module_NMM_CORE_SETUP.o module_NESTING.o module_NMM_INTEGRATE.o module_PARENT_CHILD_CPL_COMP.o module_ATM_GRID_COMP.o module_ATM_DRIVER_INTERNAL_STATE.o module_ATM_DRIVER_COMP.o':g \
        > ../src/atmos/drivers/makefile
       cat ../src/main/Makefile_stub.IN | sed s:_MODULES_:'cd $(GFS_DYN) \&\& gmake -f makefile_stub ; cp $(GFS_DYN)/gfs_dynamics_grid_comp_mod.mod $(MODDIR)/modstub/gfs ; cp $(GFS_PHY)/gfs_physics_grid_comp_mod.mod $(MODDIR)/modstub/gfs ; cd $(GFS_IO) \&\& gmake -f makefile_stub ; cp $(GFS_IO)/module_write_routines_gfs.mod $(MODDIR)/modstub/gfs ; cp $(GFS_IO)/module_write_grid_comp_gfs.mod $(MODDIR)/modstub/gfs ; cd $(GEFS_CPL) \&\& gmake -f makefile_stub ; cp $(GEFS_CPL)/gefs_cplcomp_esmfmod.mod $(MODDIR)/modstub/gfs ; cd $(ATM) \&\& gmake -f makefile_stub ; cp module_gfs_core_setup.mod $(MODDIR)/modstub/gfs ; cp module_gocart_integrate.mod $(MODDIR)/modstub/gfs ; cp module_gfs_integrate.mod $(MODDIR)/modstub/gfs':g \
        | sed s:_NEMSCORE_:'gfs':g \
        | sed s:_CLEANMODULES_:'cd $(GFS_DYN) \&\& gmake -f makefile_stub clean ; cd $(GFS_PHY) \&\& gmake -f makefile_stub clean ; cd $(GFS_IO) \&\& gmake -f makefile_stub clean ; cd $(GFS_UTIL) \&\& gmake -f makefile_stub clean ; cd $(GEFS_CPL) \&\& gmake -f makefile_stub clean ; cd $(ATM) \&\& gmake -f makefile_stub clean':g \
        > ../src/main/Makefile_stub
        cat ../src/main/Makefile_main.IN | sed s:_MODULES_:'cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src install GOCART_MODE=$(GOCART_MODE) ; cd $(GFS_UTIL) \&\& gmake depend \&\& gmake ; cd $(GFS_PHY) \&\& gmake depend \&\& gmake ; cd $(GEFS_CPL) \&\& gmake ; cd $(NMM_UTIL) \&\& gmake depend \&\& gmake ; cd $(NMM_IO) \&\& gmake depend \&\& gmake ; cd $(NMM_DYN) \&\& gmake depend \&\& gmake ; cd $(NMM_PHY) \&\& gmake depend \&\& gmake ; cd $(ATM) \&\& gmake depend \&\& gmake':g \
        | sed s:_NEMSCORE_:'nmm':g \
        | sed s:_FLAGS_:'$(FFLAG_NMM_MAIN)':g \
        | sed s:_CLEANMODULES_:'cd $(GFS_PHY) \&\& gmake ; cd $(NMM_UTIL) \&\& gmake clean ; cd $(NMM_IO) \&\& gmake clean  ; cd $(NMM_DYN) \&\& gmake clean ; cd $(NMM_PHY) \&\& gmake clean ; cd $(GOCART) \&\& export ESMADIR=`pwd` \&\& gmake -C src realclean GOCART_MODE=$(GOCART_MODE) ; cd $(ATM) \&\& gmake clean':g \
        > ../src/main/Makefile_main
fi
