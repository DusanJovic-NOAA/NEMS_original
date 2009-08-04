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
    rm -f ../lib/libstub/*.a
    rm -f ../mod/modstub/gfs/*.mod
    rm -f ../mod/modstub/nmm/*.mod
    rm -f ../mod/modstub/gfs/*.o
    rm -f ../mod/modstub/nmm/*.o
    cd ../src/main && make -f Makefile_stub clean
echo 'NEMS cleaned'

exit
