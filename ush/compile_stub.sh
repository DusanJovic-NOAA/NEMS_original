#!/bin/ksh

################################################################################
# Authors:    Ed Colon
################################################################################
#
#
#
set +x
#####
# ------ NCEP IBM SP ----------------
#=======================================================================
    cd ../src/main/ && make -f Makefile_stub

echo 'NEMS built'

exit
