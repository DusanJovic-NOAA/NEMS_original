#!/bin/ksh

################################################################################
# Script name:         prepare
# Script description:  Sets up Load Leveler commands for serial batch jobs
#                      in /scripts needed to run the WRF-NMM launcher
#
# Authors:    Ed Colon
# CODE REVISION:
# 06Nov2009   Sarah Lu, change make to gmake
# 06Nov2009   Sarah Lu, NEMSTOPDIR passed onto compile_settings
# 06Nov2009   Sarah Lu,   GOCART_MODE passed onto compile_settings
################################################################################
#
#  Used to set the following Load Leveler commands in serial batch jobs in the
#  /scripts subdirectory
#
#
#-- *** IMPORTANT *** IMPORTANT *** IMPORTANT *** IMPORTANT ***
#
#   Make sure that the settings used here are consistent with those set in prelaunch!!!
#
###############################################################################################################
set -e

OS=`uname`
echo "platform is "${OS}
settings_file=compile_settings
if [[ -a "$settings_file" ]]
then
echo "modifying compile settings file"
rm $settings_file
fi

# ------ NEMSTOPDIR -----------------
NEMSTOPDIR=`cd ..; pwd`
echo "NEMSTOPDIR = ${NEMSTOPDIR}" > temp.mk

# ------ NCEP IBM SP ----------------
    if [[ $1 = 'nmm' ]]
    then
    echo "Setup for NMM core"

### the GOCART option is not yet available for NMMb core
    echo "GOCART option is turned off"
    echo "GOCART_MODE=stub" >> temp.mk

    if [[ $2 = 'gfs_physics' ]]
    then
    cat compile_options/nmm_compile_options_${OS}.IN | sed s:_GFS_PHY_LIB_:'$(LIBDIR)/lib/gfs_phys.a':g > compile_options/nmm_compile_options_${OS} 
####cat compile_options/gen_compile_options_${OS} compile_options/nmm_compile_options_${OS} > compile_settings
    cat temp.mk compile_options/gen_compile_options_${OS} compile_options/nmm_compile_options_${OS} > compile_settings
    cat makefile_templates/main_makefile_nmm_stub.IN | sed s:_GFS_PHYSICS_:#:g \
	| sed s:_GFS_PHY_LIB_LOC_:#:g >  makefile_templates/main_makefile_nmm_stub
    cp -f makefile_templates/main_makefile_nmm_stub ../src/main/Makefile_stub
###cat makefile_templates/main_makefile_nmm.IN | sed s:_GFS_PHYSICS_:'cd $(GFS_PHY) \&\& make -f makefile_gfs_full_physics':g > makefile_templates/main_makefile_nmm 
    cat makefile_templates/main_makefile_nmm.IN | sed s:_GFS_PHYSICS_:'cd $(GFS_PHY) \&\& gmake -f makefile_gfs_full_physics':g > makefile_templates/main_makefile_nmm 
    else 
    echo 'ok'
    cat compile_options/nmm_compile_options_${OS}.IN | sed s:_GFS_PHY_LIB_:'$(LIBDIR)/libstub/gfs_phys_stub.a':g > compile_options/nmm_compile_options_${OS}
####cat compile_options/gen_compile_options_${OS} compile_options/nmm_compile_options_${OS} > compile_settings
    cat temp.mk compile_options/gen_compile_options_${OS} compile_options/nmm_compile_options_${OS} > compile_settings
####cat makefile_templates/main_makefile_nmm_stub.IN | sed s:_GFS_PHYSICS_:'cd $(GFS_PHY) \&\& make -f makefile_stub':g \
    cat makefile_templates/main_makefile_nmm_stub.IN | sed s:_GFS_PHYSICS_:'cd $(GFS_PHY) \&\& gmake -f makefile_stub':g \
        | sed s:_GFS_PHY_LIB_LOC_:'cp $(GFS_PHY)/gfs_physics_grid_comp_mod.mod $(MODDIR)/modstub/gfs':g >  makefile_templates/main_makefile_nmm_stub
    cat makefile_templates/main_makefile_nmm.IN | sed s:_GFS_PHYSICS_:#:g > makefile_templates/main_makefile_nmm
    fi 
    cp -f makefile_templates/main_makefile_nmm_stub ../src/main/Makefile_stub
    cp -f makefile_templates/main_makefile_nmm ../src/main/Makefile_main
    cp -f makefile_templates/nmm_atmos_makefile ../src/main/atmos/makefile
    cp -f makefile_templates/makefile_nmm_stub ../src/main/atmos/makefile_stub 
    elif [[ $1 = 'gfs' ]]
    then
    echo "Setup for GFS core"

    if [[ $2 = 'gocart' ]]
    then
      echo "GOCART option is turned on"
      echo "GOCART_MODE=full" >> temp.mk
    else
      echo "GOCART option is turned off"
      echo "GOCART_MODE=stub" >> temp.mk
    fi

####cat compile_options/gen_compile_options_${OS} compile_options/gfs_compile_options_${OS} > compile_settings
    cat temp.mk compile_options/gen_compile_options_${OS} compile_options/gfs_compile_options_${OS} > compile_settings
    cp -f makefile_templates/main_makefile_gfs ../src/main/Makefile_main
    cp -f makefile_templates/main_makefile_gfs_stub ../src/main/Makefile_stub
    cp -f makefile_templates/gfs_atmos_makefile ../src/main/atmos/makefile
    cp -f makefile_templates/makefile_gfs_stub ../src/main/atmos/makefile_stub
    else
    echo "Invalid core selection"
    fi
    ./clean_stub.sh
    ./compile_stub.sh	
    rm -rf temp.mk
#=======================================================================
exit
