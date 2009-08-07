#!/bin/ksh

date > RegressionTests.log
echo "Start Regression test" >> RegressionTests.log

# RTPWD - path with previous stored data
export RTPWD=/meso/noscrub/wx20rv/REGRESSION_TEST
export ME=`whoami`
#
export PATHRT=`pwd`
cd ../../
export PATHTR=`pwd`
rm -rf /ptmp/${ME}/RT

clear;echo;echo

####################################################################################################
# Clean and compile NMMB core
####################################################################################################
echo "Preparing NMMB core for regression tests"
printf %s "Compiling NMMB core (this will take ~10 minutes)......."
cd ush

./clean.sh > /dev/null 2>&1
./compile_configure.sh nmm > 1.1 2>&1
./compile.sh > /dev/null 2>&1
echo "   NMMB core Compiled";echo;echo

cd $PATHRT

####################################################################################################
# Submit test 1
####################################################################################################

rm -f err out configfile_global configfile_regional runglobal runregional  gfs_fcst_run  gfs_ll

cat runglobal.IN         | sed s:_TPN_:32:g         \
                         | sed s:_GS_:#:g           \
                         | sed s:_ME_:${ME}:g       \
                         | sed s:_SRCD_:${PATHTR}:g \
                         | sed s:_DIR_:NMM_CNTRL:g  >  runglobal

cat configfile_global.IN | sed s:_INPES_:06:g       \
                         | sed s:_WTPG_:2:g         \
                         | sed s:_FCSTL_:48:g       \
                         | sed s:_NEMSI_:false:g    \
                         | sed s:_RSTRT_:false:g    \
                         | sed s:_RGS_:false:g      \
                         | sed s:_WGS_:true:g       >  configfile_global

job_id=`llsubmit runglobal 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 1" >> RegressionTests.log
echo "Test 1"
echo "Compare results with previous trunk version" >> RegressionTests.log
echo "Compare results with previous trunk version"
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 1 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 1 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 1 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 1 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 1 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 1 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 1 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f /ptmp/${ME}/RT/NMM_CNTRL/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i /ptmp/${ME}/RT/NMM_CNTRL/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " /ptmp/${ME}/RT/NMM_CNTRL/$i " output file" >> RegressionTests.log
  echo "Missing " /ptmp/${ME}/RT/NMM_CNTRL/$i " output file"
 (echo;echo " Test 1 failed ")>> RegressionTests.log
  echo;echo " Test 1 failed "
  exit

fi

done

echo " Test 1 passed " >> RegressionTests.log
echo " Test 1 passed "
export timing1=`grep total_tim $PATHRT/err | tail -1 | awk '{ print $5 }'`
export timingc=`cat /meso/noscrub/wx20rv/REGRESSION_TEST/NMMB_glob/timing.txt`
(echo;echo " Original timing: " $timingc " , test1 timing: " $timing1)>> RegressionTests.log
 echo;echo " Original timing: " $timingc " , test1 timing: " $timing1

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 2
####################################################################################################

rm -f err out configfile_global configfile_regional runglobal runregional  gfs_fcst_run  gfs_ll

cat runglobal.IN         | sed s:_TPN_:32:g         \
                         | sed s:_GS_:#:g           \
                         | sed s:_ME_:${ME}:g       \
                         | sed s:_SRCD_:${PATHTR}:g \
                         | sed s:_DIR_:NMM_NEMSIO:g >  runglobal

cat configfile_global.IN | sed s:_INPES_:06:g       \
                         | sed s:_WTPG_:2:g         \
                         | sed s:_FCSTL_:48:g       \
                         | sed s:_NEMSI_:true:g     \
                         | sed s:_RSTRT_:false:g    \
                         | sed s:_RGS_:false:g      \
                         | sed s:_WGS_:false:g      >  configfile_global

job_id=`llsubmit runglobal 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 2" >> RegressionTests.log
echo "Test 2"
echo "Test NEMSIO as input file" >> RegressionTests.log
echo "Test NEMSIO as input file"
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 2 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 2 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 2 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 2 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 2 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 2 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 2 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f /ptmp/${ME}/RT/NMM_NEMSIO/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i /ptmp/${ME}/RT/NMM_NEMSIO/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " /ptmp/${ME}/RT/NMM_NEMSIO/$i " output file" >> RegressionTests.log
  echo "Missing " /ptmp/${ME}/RT/NMM_NEMSIO/$i " output file"
 (echo;echo " Test 2 failed ")>> RegressionTests.log
  echo;echo " Test 2 failed "
  exit

fi

done

echo " Test 2 passed " >> RegressionTests.log
echo " Test 2 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 3
####################################################################################################

rm -f err out configfile_global configfile_regional runglobal runregional  gfs_fcst_run  gfs_ll

cat runglobal.IN         | sed s:_TPN_:32:g         \
                         | sed s:_GS_:#:g           \
                         | sed s:_ME_:${ME}:g       \
                         | sed s:_SRCD_:${PATHTR}:g \
                         | sed s:_DIR_:NMM_REST:g   >  runglobal

cat configfile_global.IN | sed s:_INPES_:06:g       \
                         | sed s:_WTPG_:2:g         \
                         | sed s:_FCSTL_:48:g       \
                         | sed s:_NEMSI_:false:g    \
                         | sed s:_RSTRT_:true:g     \
                         | sed s:_RGS_:false:g      \
                         | sed s:_WGS_:false:g      >  configfile_global

job_id=`llsubmit runglobal 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 3" >> RegressionTests.log
echo "Test 3"
echo "Test restart run" >> RegressionTests.log
echo "Test restart run"
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 3 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 3 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 3 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 3 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 3 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 3 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 3 results ...."

#
# NOTE! first hitory files differ because of PSFC, others should be identical:
#
for i in nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f /ptmp/${ME}/RT/NMM_REST/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i /ptmp/${ME}/RT/NMM_REST/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " /ptmp/${ME}/RT/NMM_REST/$i " output file" >> RegressionTests.log
  echo "Missing " /ptmp/${ME}/RT/NMM_REST/$i " output file"
 (echo;echo " Test 3 failed ")>> RegressionTests.log
  echo;echo " Test 3 failed "
  exit

fi

done

echo " Test 3 passed " >> RegressionTests.log
echo " Test 3 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 4
####################################################################################################

rm -f err out configfile_global configfile_regional runglobal runregional  gfs_fcst_run  gfs_ll

cat runglobal.IN         | sed s:_TPN_:16:g         \
                         | sed s:_GS_::g            \
                         | sed s:_ME_:${ME}:g       \
                         | sed s:_SRCD_:${PATHTR}:g \
                         | sed s:_DIR_:NMM_DECOMP:g >  runglobal

cat configfile_global.IN | sed s:_INPES_:03:g       \
                         | sed s:_WTPG_:1:g         \
                         | sed s:_FCSTL_:24:g       \
                         | sed s:_NEMSI_:false:g    \
                         | sed s:_RSTRT_:false:g    \
                         | sed s:_RGS_:true:g       \
                         | sed s:_WGS_:false:g      >  configfile_global

job_id=`llsubmit runglobal 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 4" >> RegressionTests.log
echo "Test 4"
echo "Test different decomposition" >> RegressionTests.log
echo "Test different decomposition"
(echo "NMMB, global case";echo;echo)>> RegressionTests.log
 echo "NMMB, global case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 4 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 4 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 4 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 4 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 4 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 4 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 4 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 \
         nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f /ptmp/${ME}/RT/NMM_DECOMP/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_glob/$i /ptmp/${ME}/RT/NMM_DECOMP/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " /ptmp/${ME}/RT/NMM_DECOMP/$i " output file" >> RegressionTests.log
  echo "Missing " /ptmp/${ME}/RT/NMM_DECOMP/$i " output file"
 (echo;echo " Test 4 failed ")>> RegressionTests.log
  echo;echo " Test 4 failed "
  exit

fi

done

echo " Test 4 passed " >> RegressionTests.log
echo " Test 4 passed "

sleep 4
clear;echo;echo


####################################################################################################
# Submit test 5
####################################################################################################

rm -f err out configfile_global configfile_regional runglobal runregional  gfs_fcst_run  gfs_ll

cat runregional.IN         | sed s:_TPN_:32:g           \
                           | sed s:_GS_:#:g             \
                           | sed s:_ME_:${ME}:g         \
                           | sed s:_SRCD_:${PATHTR}:g   \
                           | sed s:_DIR_:NMM_REG_CTL:g  >  runregional

cat configfile_regional.IN | sed s:_INPES_:06:g         \
                           | sed s:_WTPG_:2:g           \
                           | sed s:_FCSTL_:48:g         \
                           | sed s:_NEMSI_:false:g      \
                           | sed s:_RSTRT_:false:g      \
                           | sed s:_RGS_:false:g        \
                           | sed s:_WGS_:true:g         >  configfile_regional

job_id=`llsubmit runregional 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 5" >> RegressionTests.log
echo "Test 5"
echo "Compare results with previous trunk version" >> RegressionTests.log
echo "Compare results with previous trunk version"
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 5 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 5 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 5 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 5 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 5 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 5 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 5 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
         nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
         nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f /ptmp/${ME}/RT/NMM_REG_CTL/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i /ptmp/${ME}/RT/NMM_REG_CTL/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " /ptmp/${ME}/RT/NMM_REG_CTL/$i " output file" >> RegressionTests.log
  echo "Missing " /ptmp/${ME}/RT/NMM_REG_CTL/$i " output file"
 (echo;echo " Test 5 failed ")>> RegressionTests.log
  echo;echo " Test 5 failed "
  exit

fi

done

echo " Test 5 passed " >> RegressionTests.log
echo " Test 5 passed "
export timing1=`grep total_tim $PATHRT/err | tail -1 | awk '{ print $5 }'`
export timingc=`cat /meso/noscrub/wx20rv/REGRESSION_TEST/NMMB_reg/timing.txt`
(echo;echo " Original timing: " $timingc " , test1 timing: " $timing1)>> RegressionTests.log
 echo;echo " Original timing: " $timingc " , test1 timing: " $timing1

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 6
####################################################################################################

rm -f err out configfile_global configfile_regional runglobal runregional  gfs_fcst_run  gfs_ll

cat runregional.IN         | sed s:_TPN_:32:g              \
                           | sed s:_GS_:#:g                \
                           | sed s:_ME_:${ME}:g            \
                           | sed s:_SRCD_:${PATHTR}:g      \
                           | sed s:_DIR_:NMM_REG_NEMSIO:g  >  runregional

cat configfile_regional.IN | sed s:_INPES_:06:g            \
                           | sed s:_WTPG_:2:g              \
                           | sed s:_FCSTL_:12:g            \
                           | sed s:_NEMSI_:true:g          \
                           | sed s:_RSTRT_:false:g         \
                           | sed s:_RGS_:false:g           \
                           | sed s:_WGS_:false:g           >  configfile_regional

job_id=`llsubmit runregional 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 6" >> RegressionTests.log
echo "Test 6"
echo "Test NEMSIO as input file" >> RegressionTests.log
echo "Test NEMSIO as input file"
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 6 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 6 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 6 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 6 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 6 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 6 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 6 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.009 nmm_b_history.012    \
         nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                   \
         nmm_b_history_nemsio.009 nmm_b_history_nemsio.012
 
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f /ptmp/${ME}/RT/NMM_REG_NEMSIO/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i /ptmp/${ME}/RT/NMM_REG_NEMSIO/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " /ptmp/${ME}/RT/NMM_REG_NEMSIO/$i " output file" >> RegressionTests.log
  echo "Missing " /ptmp/${ME}/RT/NMM_REG_NEMSIO/$i " output file"
 (echo;echo " Test 6 failed ")>> RegressionTests.log
  echo;echo " Test 6 failed "
  exit

fi

done

echo " Test 6 passed " >> RegressionTests.log
echo " Test 6 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 7
####################################################################################################

rm -f err out configfile_global configfile_regional runglobal runregional  gfs_fcst_run  gfs_ll

cat runregional.IN         | sed s:_TPN_:32:g           \
                           | sed s:_GS_:#:g             \
                           | sed s:_ME_:${ME}:g         \
                           | sed s:_SRCD_:${PATHTR}:g   \
                           | sed s:_DIR_:NMM_REG_RST:g  >  runregional

cat configfile_regional.IN | sed s:_INPES_:06:g         \
                           | sed s:_WTPG_:2:g           \
                           | sed s:_FCSTL_:48:g         \
                           | sed s:_NEMSI_:false:g      \
                           | sed s:_RSTRT_:true:g       \
                           | sed s:_RGS_:false:g        \
                           | sed s:_WGS_:false:g        >  configfile_regional

job_id=`llsubmit runregional 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 7" >> RegressionTests.log
echo "Test 7"
echo "Test restart run" >> RegressionTests.log
echo "Test restart run"
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 7 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 7 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 7 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 7 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 7 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 7 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 7 results ...."

#
# NOTE! first hitory files differ because of PSFC, others should be identical:
#
for i in nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f /ptmp/${ME}/RT/NMM_REG_RST/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i /ptmp/${ME}/RT/NMM_REG_RST/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " /ptmp/${ME}/RT/NMM_REG_RST/$i " output file" >> RegressionTests.log
  echo "Missing " /ptmp/${ME}/RT/NMM_REG_RST/$i " output file"
 (echo;echo " Test 7 failed ")>> RegressionTests.log
  echo;echo " Test 7 failed "
  exit

fi

done

echo " Test 7 passed " >> RegressionTests.log
echo " Test 7 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 8
####################################################################################################

rm -f err out configfile_global configfile_regional runglobal runregional  gfs_fcst_run  gfs_ll

cat runregional.IN         | sed s:_TPN_:16:g              \
                           | sed s:_GS_::g                 \
                           | sed s:_ME_:${ME}:g            \
                           | sed s:_SRCD_:${PATHTR}:g      \
                           | sed s:_DIR_:NMM_REG_DECOMP:g  >  runregional

cat configfile_regional.IN | sed s:_INPES_:03:g            \
                           | sed s:_WTPG_:1:g              \
                           | sed s:_FCSTL_:12:g            \
                           | sed s:_NEMSI_:false:g         \
                           | sed s:_RSTRT_:false:g         \
                           | sed s:_RGS_:true:g            \
                           | sed s:_WGS_:false:g           >  configfile_regional

job_id=`llsubmit runregional 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 8" >> RegressionTests.log
echo "Test 8"
echo "Test different decomposition" >> RegressionTests.log
echo "Test different decomposition"
(echo "NMMB, regional case";echo;echo)>> RegressionTests.log
 echo "NMMB, regional case";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $8 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST 8 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 8 is running,            ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST 8 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 8 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 8 is finished,           ID: " $job_id " Status: " $status  ", Finished " $FnshHrs "hours"
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 8 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 8 results ...."

for i in nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.009 nmm_b_history.012

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f /ptmp/${ME}/RT/NMM_REG_DECOMP/$i ] ; then

  d=`cmp ${RTPWD}/NMMB_reg/$i /ptmp/${ME}/RT/NMM_REG_DECOMP/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " /ptmp/${ME}/RT/NMM_REG_DECOMP/$i " output file" >> RegressionTests.log
  echo "Missing " /ptmp/${ME}/RT/NMM_REG_DECOMP/$i " output file"
 (echo;echo " Test 8 failed ")>> RegressionTests.log
  echo;echo " Test 8 failed "
  exit

fi

done

echo " Test 8 passed " >> RegressionTests.log
echo " Test 8 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Clean and compile GFS core
####################################################################################################
echo "Preparing GFS core for regression tests" >> RegressionTests.log
echo "Preparing GFS core for regression tests"
printf %s "Compiling GFS core (this will take ~10 minutes)......."
cd ${PATHTR}/ush

./clean.sh > /dev/null 2>&1
./compile_configure.sh gfs > 1.1 2>&1
./compile.sh > /dev/null 2>&1
echo "   GFS core Compiled";echo;echo

cd $PATHRT

####################################################################################################
# Submit test 9
####################################################################################################

export RUNDIR=/ptmp/${ME}/RT/GFS_32
mkdir -p ${RUNDIR}
cp configfile_global_gfs ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configfile_global configfile_regional runglobal runregional gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:32:g         \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:32:g         \
                    | sed s:_PE1_:30:g           \
                    | sed s:_WPG_:2:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 9" >> RegressionTests.log
echo "Test 9"
echo "Compare results with previous trunk version" >> RegressionTests.log
echo "Compare results with previous trunk version"
(echo "GFS, 32 proc, 1 thread";echo;echo)>> RegressionTests.log
 echo "GFS, 32 proc, 1 thread";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 9 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 9 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 9 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 9 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 9 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 9 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 9 results ...."

for i in sigf003_nemsio sigf006_nemsio sigf012_nemsio sigf024_nemsio sigf048_nemsio \
         sfcf003_nemsio sfcf006_nemsio sfcf012_nemsio sfcf024_nemsio sfcf048_nemsio \
         flxf003_nemsio flxf006_nemsio flxf012_nemsio flxf024_nemsio flxf048_nemsio

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 9 failed ")>> RegressionTests.log
  echo;echo " Test 9 failed "
  exit

fi

done

echo " Test 9 passed " >> RegressionTests.log
echo " Test 9 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 10
####################################################################################################

export RUNDIR=/ptmp/${ME}/RT/GFS_60
mkdir -p ${RUNDIR}
cp configfile_global_gfs ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configfile_global configfile_regional runglobal runregional gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:60:g         \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:60:g         \
                    | sed s:_PE1_:58:g           \
                    | sed s:_WPG_:2:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 10" >> RegressionTests.log
echo "Test 10"
echo "Test different decomposition" >> RegressionTests.log
echo "Test different decomposition"
(echo "GFS, 60 proc, 1 thread";echo;echo)>> RegressionTests.log
 echo "GFS, 60 proc, 1 thread";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 10 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 10 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 10 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 10 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 10 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 10 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 10 results ...."

for i in sigf003_nemsio sigf006_nemsio sigf012_nemsio sigf024_nemsio sigf048_nemsio \
         sfcf003_nemsio sfcf006_nemsio sfcf012_nemsio sfcf024_nemsio sfcf048_nemsio \
         flxf003_nemsio flxf006_nemsio flxf012_nemsio flxf024_nemsio flxf048_nemsio

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 10 failed ")>> RegressionTests.log
  echo;echo " Test 10 failed "
  exit

fi

done

echo " Test 10 passed " >> RegressionTests.log
echo " Test 10 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 11
####################################################################################################

export RUNDIR=/ptmp/${ME}/RT/GFS_16
mkdir -p ${RUNDIR}
cp configfile_global_gfs ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configfile_global configfile_regional runglobal runregional gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:16:g         \
                    | sed s:_THRDS_:2:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:16:g         \
                    | sed s:_PE1_:14:g           \
                    | sed s:_WPG_:2:g            \
                    | sed s:_THRDS_:2:g          \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:2:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 11" >> RegressionTests.log
echo "Test 11"
echo "Test threads" >> RegressionTests.log
echo "Test threads"
(echo "GFS, 16 proc, 2 threads";echo;echo)>> RegressionTests.log
 echo "GFS, 16 proc, 2 threads";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 11 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 11 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 11 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 11 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 11 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 11 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 11 results ...."

for i in sigf003_nemsio sigf006_nemsio sigf012_nemsio sigf024_nemsio sigf048_nemsio \
         sfcf003_nemsio sfcf006_nemsio sfcf012_nemsio sfcf024_nemsio sfcf048_nemsio \
         flxf003_nemsio flxf006_nemsio flxf012_nemsio flxf024_nemsio flxf048_nemsio

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 11 failed ")>> RegressionTests.log
  echo;echo " Test 11 failed "
  exit

fi

done

echo " Test 11 passed " >> RegressionTests.log
echo " Test 11 passed "

sleep 4
clear;echo;echo

####################################################################################################
# Submit test 12
####################################################################################################

export RUNDIR=/ptmp/${ME}/RT/GFS_01
mkdir -p ${RUNDIR}
cp configfile_global_gfs ${RUNDIR}/configure_file
cp ${RTPWD}/GFS/gfsanl.2009072400 ${RUNDIR}/.
cp ${RTPWD}/GFS/sfcanl.2009072400 ${RUNDIR}/.
rm -f err out configfile_global configfile_regional runglobal runregional gfs_fcst_run  gfs_ll

cat gfs_ll.IN       | sed s:_TASKS_:1:g          \
                    | sed s:_THRDS_:1:g          >  gfs_ll

cat gfs_fcst_run.IN | sed s:_TASKS_:1:g          \
                    | sed s:_PE1_:1:g            \
                    | sed s:_WPG_:0:g            \
                    | sed s:_THRDS_:1:g          \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_NDAYS_:1:g          >  gfs_fcst_run

job_id=`llsubmit gfs_ll 2>&1 | grep submitted`
job_id=`echo $job_id | cut -d\" -f2 | cut -d. -f1,5 `

echo "Test 12" >> RegressionTests.log
echo "Test 12"
echo "Test single processor" >> RegressionTests.log
echo "Test single processor"
(echo "GFS, 1 proc, 1 threads";echo;echo)>> RegressionTests.log
 echo "GFS, 1 proc, 1 threads";echo;echo

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

sleep 60
export status=`llq| grep $job_id | awk '{ print$5}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST 12 is waiting in a queue, ID: " $job_id " Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST 12 is running,            ID: " $job_id " Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST 12 is ready to run,       ID: " $job_id " Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST 12 is finished,           ID: " $job_id " Status: " $status
else                          echo $n "min. TEST 12 is finished,           ID: " $job_id " Status: " $status
fi

job_running=`llq | grep $job_id | wc -l`
  (( n=n+1 ))
done

(echo;echo;echo "Checking test 12 results ....")>> RegressionTests.log
 echo;echo;echo "Checking test 12 results ...."

for i in sigf003_nemsio sigf006_nemsio sigf012_nemsio sigf024_nemsio \
         sfcf003_nemsio sfcf006_nemsio sfcf012_nemsio sfcf024_nemsio \
         flxf003_nemsio flxf006_nemsio flxf012_nemsio flxf024_nemsio

do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/GFS/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test 12 failed ")>> RegressionTests.log
  echo;echo " Test 12 failed "
  exit

fi

done

echo " Test 12 passed " >> RegressionTests.log
echo " Test 12 passed "


####################################################################################################
####################################################################################################
####################################################################################################


rm -f err out configfile_global configfile_regional runglobal runregional gfs_fcst_run  gfs_ll configure_file
cd ${PATHTR}/ush
./clean.sh > /dev/null 2>&1

date

exit
