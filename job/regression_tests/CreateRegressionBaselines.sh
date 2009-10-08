#!/bin/ksh

date > RegressionBaselines.log
echo "Start Regression Baseline Creation" >> RegressionBaselines.log

# RTPWD - path with previous stored data
export ME=`whoami`
#
export RTPWD=/meso/noscrub/wx20rv/REGRESSION_TEST
export PATHRT=`pwd`
cd ../../
export PATHTR=`pwd`
rm -rf /ptmp/${ME}/RT

clear;echo;echo

####################################################################################################
# Clean and compile NMMB core
####################################################################################################
echo "Preparing NMMB core for baselines"
printf %s "Compiling NMMB core (this will take ~10 minutes)......."
cd ush

./clean_stub.sh > /dev/null 2>&1
./clean.sh > /dev/null 2>&1
./compile_configure.sh nmm > 1.1 2>&1
./compile.sh > /dev/null 2>&1
echo "   NMMB core Compiled";echo;echo

cd $PATHRT

####################################################################################################
# Submit baselines 1
####################################################################################################

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_glob_ll.IN         | sed s:_TPN_:32:g         \
                         | sed s:_THRD_:1:g         \
                         | sed s:_GS_:#:g           \
                         | sed s:_RUND_:/ptmp/${ME}/RT/NMM_CNTRL:g           \
			 | sed s:_RTPWD_:${RTPWD}:g \
                         | sed s:_SRCD_:${PATHTR}:g \
                         | sed s:_DIR_:NMM_CNTRL:g  >  nmm_glob_ll

cat nmm_glob_conf.IN     | sed s:_INPES_:06:g       \
                         | sed s:_WTPG_:2:g         \
                         | sed s:_FCSTL_:48:g       \
                         | sed s:_NEMSI_:false:g    \
                         | sed s:_RSTRT_:false:g    \
                         | sed s:_gfsP_:false:g     \
                         | sed s:_RGS_:false:g      \
                         | sed s:_WGS_:true:g       >  configure_file

job_id=`llsubmit nmm_glob_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Global Baseline" >> RegressionBaselines.log
echo "Global Baseline"

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. Baseline 1 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. Baseline 1 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. Baseline 1 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. Baseline 1 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. Baseline 1 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do

if [ -f /ptmp/${ME}/RT/NMM_CNTRL/$i ] ; then

  echo "All history files present" >> RegressionBaselines.log
  echo "All history files present"

else

  echo "Missing " /ptmp/${ME}/RT/NMM_CNTRL/$i " output file" >> RegressionBaselines.log
  echo "Missing " /ptmp/${ME}/RT/NMM_CNTRL/$i " output file"

  exit

fi

done


sleep 4
clear;echo;echo

####################################################################################################
# Submit baseline 2
####################################################################################################

rm -f err out configure_file nmm_glob_ll nmm_reg_ll  gfs_fcst_run  gfs_ll

cat nmm_reg_ll.IN         | sed s:_TPN_:32:g           \
                          | sed s:_THRD_:1:g         \
                           | sed s:_GS_:#:g             \
                           | sed s:_RUND_:/ptmp/${ME}/RT/NMM_REG_CTL:g           \
                           | sed s:_RTPWD_:${RTPWD}:g   \
                           | sed s:_SRCD_:${PATHTR}:g   \
                           | sed s:_DIR_:NMM_REG_CTL:g  > nmm_reg_ll

cat nmm_reg_conf.IN | sed s:_INPES_:06:g         \
                           | sed s:_WTPG_:2:g           \
                           | sed s:_FCSTL_:48:g         \
                           | sed s:_NEMSI_:false:g      \
                           | sed s:_RSTRT_:false:g      \
                           | sed s:_gfsP_:false:g     \
                           | sed s:_RGS_:false:g        \
                           | sed s:_WGS_:true:g         >  configure_file

job_id=`llsubmit nmm_reg_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "REGIONAL BASELINE" >> RegressionBaselines.log
echo "REGIONAL BASELINE"

job_running=1

n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. BASELINE 2 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. BASELINE 2 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. BASELINE 2 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. BASELINE 2 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. BASELINE 2 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done


for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do

if [ -f /ptmp/${ME}/RT/NMM_REG_CTL/$i ] ; then

  echo "All history files present" >> RegressionBaselines.log
  echo "All history files present"
else
  echo "Missing " /ptmp/${ME}/RT/NMM_REG_CTL/$i " output file" >> RegressionBaselines.log
  echo "Missing " /ptmp/${ME}/RT/NMM_REG_CTL/$i " output file"
  exit
fi

done

sleep 4
clear;echo;echo

####################################################################################################
# Clean and compile GFS core
####################################################################################################
echo "Preparing GFS core for baselines" >> RegressionBaselines.log
echo "Preparing GFS core for baselines"
printf %s "Compiling GFS core (this will take ~10 minutes)......."
cd ${PATHTR}/ush

./clean_stub.sh > /dev/null 2>&1
./clean.sh > /dev/null 2>&1
./compile_configure.sh gfs > 1.1 2>&1
./compile.sh > /dev/null 2>&1
echo "   GFS core Compiled";echo;echo

cd $PATHRT

####################################################################################################
####################################################################################################
# Submit baselines 3
####################################################################################################

export RUNDIR=/ptmp/${ME}/RT/GFS_32
mkdir -p ${RUNDIR}
cp gfs_configfile ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configure_file nmm_glob_ll nmm_reg_ll gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:32:g         \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:32:g         \
                    | sed s:_PE1_:30:g           \
                    | sed s:_WPG_:2:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_NSOUT_:1:g          \
                    | sed s:_QUILTING_:.true.:g          \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "GFS baseline" >> RegressionBaselines.log
echo "GFS baseline"
(echo "GFS, 32 proc, 1 thread";echo;echo)>> RegressionBaselines.log
 echo "GFS, 32 proc, 1 thread";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. BASELINE 3 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. BASELINE 3 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. BASELINE 3 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. BASELINE 3 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. BASELINE 3 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done
for i in sigf03 sigf06 sigf09 sigf12 sigf24 sigf48 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 sfcf48 \
         flxf03 flxf06 sigf09 flxf12 flxf24 flxf48

do

if [ -f /ptmp/${ME}/RT/GFS_32/$i ] ; then
  echo "All history files present" >> RegressionBaselines.log
  echo "All history files present"
else
  echo "Missing " /ptmp/${ME}/RT/GFS_32/$i " output file" >> RegressionBaselines.log
  echo "Missing " /ptmp/${ME}/RT/GFS_32/$i " output file"
  exit
fi

done

sleep 4
clear;echo;echo

echo "Finished creating new NEMS regression test baselines located at " /ptmp/${ME}/RT/ "."

exit
