#!/bin/ksh
set -ua

export GEN_ENSEMBLE=${GEN_ENSEMBLE:-0}
echo "GEN_ENSEMBLE=" $GEN_ENSEMBLE

mkdir -p ${RUNDIR}

if [ $GEN_ENSEMBLE = 0 ] ; then

####################################################################################################
# For the stand alone gen regression tests.
####################################################################################################

####################################################################################################
# Make configure and run files
####################################################################################################


echo 'RUNDIR=' $RUNDIR

cat gen_fcst_run_GEN_m1.IN \
                    | sed s:_SRCDIR_:${PATHTR}:g \
                    | sed s:_RUNDIR_:${RUNDIR}:g > gen_fcst_run

cp atmos.configure_gen ${RUNDIR}/atmos.configure

else

####################################################################################################
# For the concurrency ensemble GEN regression test.
####################################################################################################

cd $PATHRT

echo 'RUNDIR=' $RUNDIR

cat gen_fcst_run_GEN_m4.IN \
                    | sed s:_SRCDIR_:${PATHTR}:g \
                    | sed s:_RUNDIR_:${RUNDIR}:g > gen_fcst_run

cp atmos.configure_gen ${RUNDIR}/atmos.configure

fi

####################################################################################################
# Submit test
####################################################################################################


JBNME=RT_${TEST_NR}_$$

cat gen_ll.IN       | sed s:_JBNME_:${JBNME}:g   \
                    | sed s:_CLASS_:${CLASS}:g   \
                    | sed s:_GROUP_:${GROUP}:g   \
                    | sed s:_ACCNR_:${ACCNR}:g   \
                    | sed s:_WLCLK_:${WLCLK}:g   \
                    | sed s:_TASKS_:${TASKS}:g   \
                    | sed s:_THRDS_:${THRD}:g    >  gen_ll

cat gen_fcst_run >> gen_ll

llsubmit gen_ll 2>&1 | grep submitted > /dev/null

echo "Test ${TEST_NR}" >> RegressionTests.log
echo "Test ${TEST_NR}"
echo ${TEST_DESCR} >> RegressionTests.log
echo ${TEST_DESCR}
(echo "GEN, ${TASKS} proc, ${THRD} thread";echo;echo)>> RegressionTests.log
 echo "GEN, ${TASKS} proc, ${THRD} thread";echo;echo

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

(echo;echo;echo "Checking test ${TEST_NR} results ....")>> RegressionTests.log
 echo;echo;echo "Checking test ${TEST_NR} results ...."

d=`grep 'Test run in the gen finalize routine' ${RUNDIR}/NEMS.out | wc -l`

  if [[ $d -eq 0 ]] ; then
   (echo " ......NOT OK" ; echo " Failed!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo " Failed!   " ; exit 2
  fi

echo " Test ${TEST_NR} passed " >> RegressionTests.log
echo " Test ${TEST_NR} passed "


(echo;echo;echo;)>> RegressionTests.log
 echo;echo;echo;


sleep 4
clear;echo;echo

####################################################################################################
# End test
####################################################################################################

exit 0
