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
# ------ NCEP IBM SP ----------------
#=======================================================================
if [ ! -a ../src/main/Makefile_main ] ; then
   echo Makefile_main not found. Make sure clean.sh has not previously been invoked. Exiting.
   exit
fi
    rm -f configure_flags
    cd ../src/main && gmake -f Makefile_main clean
    rm -f Makefile_main
    cd ../../build 
    rm -f compile_settings
echo 'NEMS cleaned'
exit
