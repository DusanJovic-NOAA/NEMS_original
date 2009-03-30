#!/bin/ksh
set -x 
sorc_dir=$(pwd)
exec_dir=$(pwd)
#exec_dir=/ptmp/wx23sm
mkdir -p $exec_dir
make_dir=/ptmp/$(logname)/sorc/$(basename $sorc_dir)
mkdir -p $make_dir
cd $make_dir
cd $make_dir || exit 99
[ $? -ne 0 ] && exit 8
#
# rm $make_dir/*.o $make_dir/*.mod
#
#
tar -cf- -C$sorc_dir .|tar -xf-
#
export EXEC="$exec_dir/global_fcst"
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
