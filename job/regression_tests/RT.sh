#!/bin/ksh
##########################################
# temporary fix for AIX
# ratko, May 6th, 2013
##########################################

argn=$#

if [ `uname` = AIX ]; then
 echo "#!/bin/ksh" > tmp.sh
 cat RT.sh_IN >> tmp.sh
 chmod 755 tmp.sh
 cp gfs_fcst_run.IN_IBM gfs_fcst_run.IN
else
 echo "#!/bin/ksh -l" > tmp.sh
 cat RT.sh_IN >> tmp.sh
 chmod 755 tmp.sh
 cp gfs_fcst_run.IN_Linux gfs_fcst_run.IN
fi

  if [ $argn = 0 ]; then
    tmp.sh
  elif [ $argn = 1 ]; then
    tmp.sh $1
  elif [ $argn = 2 ]; then
    tmp.sh $1 $2
  fi

  rm -f tmp.sh gfs_fcst_run.IN

exit
