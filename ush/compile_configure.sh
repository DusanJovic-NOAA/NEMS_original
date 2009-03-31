#!/bin/ksh

################################################################################
# Script name:         prepare
# Script description:  Sets up Load Leveler commands for serial batch jobs
#                      in /scripts needed to run the WRF-NMM launcher
#
# Authors:    Ed Colon
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
# ------ NCEP IBM SP ----------------
    if [[ $1 = 'nmm' ]]
    then
    echo "Setup for NMM core"
    rm compile_settings
    cat compile_options/gen_compile_options compile_options/nmm_compile_options > compile_settings
    cp -f makefile_templates/main_makefile_nmm ../src/main/Makefile_main
    cp -f makefile_templates/nmm_atmos_makefile ../src/main/atmos/atmos_makefile
    elif [[ $1 = 'gfs' ]]
    then
    echo "Setup for GFS core"
    rm compile_settings
    cat compile_options/gen_compile_options compile_options/gfs_compile_options > compile_settings
    cp -f makefile_templates/main_makefile_gfs ../src/main/Makefile_main
    cp -f makefile_templates/gfs_atmos_makefile ../src/main/atmos/atmos_makefile
    else
    echo "Invalid core selection"
    fi

#=======================================================================
exit
