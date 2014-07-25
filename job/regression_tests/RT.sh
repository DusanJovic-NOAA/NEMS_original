#!/bin/ksh
##########################################
# temporary fix for WCOSS
# ratko, May 6th, 2013
##########################################

#export ACCNR=null
#rv export ACCNR=nems
#rv . ./detect_machine.sh
#rv export MACHINE_ID

#rv
export ACCNR=rm
export MACHINE_ID=zeus
#rv
if [ $MACHINE_ID = wcoss ]; then
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

./tmp.sh $*

rm -f tmp.sh gfs_fcst_run.IN

exit
