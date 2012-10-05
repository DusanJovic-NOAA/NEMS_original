#!/bin/ksh
set -ua

mkdir -p ${RUNDIR}

## specify the regression test setting

export NEMSDIR=${PATHTR}
export WORKDIR=${RUNDIR}
export REGSDIR=${RTPWD}/GFS_GOCART_POST
export PARA_CONFIG=${REGSDIR}/ngac_para_config

####################################################################################################
# Submit test
####################################################################################################

JBNME=RT_${TEST_NR}_$$

cat ngac_ll.IN      | sed s:_JBNME_:${JBNME}:g   \
                    | sed s:_CLASS_:${CLASS}:g   \
                    | sed s:_GROUP_:${GROUP}:g   \
                    | sed s:_NEMSDIR_:${NEMSDIR}:g   \
                    | sed s:_WORKDIR_:${WORKDIR}:g   \
                    | sed s:_REGSDIR_:${REGSDIR}:g   \
                    | sed s:_CONFIG_:${PARA_CONFIG}:g   \
                    | sed s:_ACCNR_:${ACCNR}:g   \
                    | sed s:_WLCLK_:${WLCLK}:g   \
                    | sed s:_TASKS_:${TASKS}:g   \
                    | sed s:_THRDS_:${THRD}:g    >  ngac_ll

llsubmit ngac_ll 2>&1 | grep submitted > /dev/null

echo "Test ${TEST_NR}" >> ${REGRESSIONTEST_LOG}
echo "Test ${TEST_NR}"
echo ${TEST_DESCR} >> ${REGRESSIONTEST_LOG}
echo ${TEST_DESCR}
(echo "GFS, ${TASKS} proc, ${THRD} thread";echo;echo)>> ${REGRESSIONTEST_LOG}
 echo "GFS, ${TASKS} proc, ${THRD} thread";echo;echo

# wait for the job to enter the queue
job_running=0
until [ $job_running -eq 1 ]
do
echo "TEST is waiting to enter the queue"
job_running=`llq -u ${LOGIN} -f %st %jn | grep ${JBNME} | wc -l`;sleep 5
done

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

export status=`llq -u ${LOGIN} -f %st %jn | grep ${JBNME} | awk '{ print $1}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST ${TEST_NR} is waiting in a queue, Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST ${TEST_NR} is running,            Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST ${TEST_NR} is ready to run,       Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
else                          echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
fi

sleep 60

job_running=`llq -u ${LOGIN} -f %st %jn | grep ${JBNME} | wc -l`
  (( n=n+1 ))
done

####################################################################################################
# Check results
####################################################################################################

(echo;echo;echo "Checking test ${TEST_NR} results ....")>> ${REGRESSIONTEST_LOG}
 echo;echo;echo "Checking test ${TEST_NR} results ...."

#
     if [ ${CREATE_BASELINE} = false ]; then
#
# --- regression test comparison ----
#

for i in ${LIST_FILES}

do
printf %s " Comparing " $i "....." >> ${REGRESSIONTEST_LOG}
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/${CNTL_DIR}/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> ${REGRESSIONTEST_LOG}
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit 2
  fi

  echo "....OK" >> ${REGRESSIONTEST_LOG}
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> ${REGRESSIONTEST_LOG}
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test ${TEST_NR} failed ")>> ${REGRESSIONTEST_LOG}
  echo;echo " Test ${TEST_NR} failed "
  exit 2

fi

done

#
     else
#
# --- create baselines
#

 echo;echo;echo "Moving set ${TEST_NR} files ...."

for i in ${LIST_FILES}
do
  printf %s " Moving " $i "....."
  if [ -f ${RUNDIR}/$i ] ; then
    cp ${RUNDIR}/${i} /stmp/${LOGIN}/REGRESSION_TEST/${CNTL_DIR}/${i}
  else
    echo "Missing " ${RUNDIR}/$i " output file"
    echo;echo " Set ${TEST_NR} failed "
    exit 2
  fi
done

# ---
     fi
# ---

echo " Test ${TEST_NR} passed " >> ${REGRESSIONTEST_LOG}
echo " Test ${TEST_NR} passed "

sleep 4
clear;echo;echo

####################################################################################################
# End test
####################################################################################################

exit 0
