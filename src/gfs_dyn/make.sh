#!/bin/ksh
set -x 
sorc_dir=$(pwd)
exec_dir=$(pwd)
#
export EXEC="$exec_dir/gfs_dynamics.a"
#
 node=`hostname | cut -d. -f1,1`
 mac=`echo $node | cut -c1-1`
 if [ $mac = b -o $mac = w ] ; then
  suf=ncep
 fi
 export F77=${suf}mpxlf
 export F90=${suf}mpxlf90
#
make -f Makefile
