#!/bin/ksh
set -ua

export GEN_ENSEMBLE=${GEN_ENSEMBLE:-0}
echo "GEN_ENSEMBLE=" $GEN_ENSEMBLE

export RT_BASEDIR=`pwd`

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
                    | sed s:_SCHEDULER_:${SCHEDULER}:g \
                    | sed s:_RUNDIR_:${RUNDIR}:g > gen_fcst_run

chmod 755 gen_fcst_run
cp gen_fcst_run ${RUNDIR}

cp atmos.configure_gen ${RUNDIR}/atmos.configure
cp ocean.configure ${RUNDIR}/ocean.configure

else

####################################################################################################
# For the concurrency ensemble GEN regression test.
####################################################################################################

cd $PATHRT

echo 'RUNDIR=' $RUNDIR

cat gen_fcst_run_GEN_m4.IN \
                    | sed s:_SRCDIR_:${PATHTR}:g \
                    | sed s:_SCHEDULER_:${SCHEDULER}:g \
                    | sed s:_RUNDIR_:${RUNDIR}:g > gen_fcst_run

chmod 755 gen_fcst_run
cp gen_fcst_run ${RUNDIR}

cp atmos.configure_gen ${RUNDIR}/atmos.configure
cp ocean.configure ${RUNDIR}/ocean.configure

fi

####################################################################################################
# Submit test
####################################################################################################


JBNME=RT_${TEST_NR}_$$

if [ $SCHEDULER = 'loadleveler' ]; then

cat gen_ll.IN       | sed s:_JBNME_:${JBNME}:g   \
                    | sed s:_CLASS_:${CLASS}:g   \
                    | sed s:_GROUP_:${GROUP}:g   \
                    | sed s:_ACCNR_:${ACCNR}:g   \
                    | sed s:_WLCLK_:${WLCLK}:g   \
                    | sed s:_TASKS_:${TASKS}:g   \
                    | sed s:_THRDS_:${THRD}:g    >  gen_ll

elif [ $SCHEDULER = 'moab' ]; then

cat gen_msub.IN         | sed s:_JBNME_:${JBNME}:g   \
                        | sed s:_WLCLK_:${WLCLK}:g   \
                        | sed s:_TPN_:${TPN}:g       \
                        | sed s:_THRD_:${THRD}:g     >  gen_msub

elif [ $SCHEDULER = 'pbs' ]; then

cat gen_qsub.IN         | sed s:_JBNME_:${JBNME}:g   \
                        | sed s:_ACCNR_:${ACCNR}:g   \
                        | sed s:_WLCLK_:${WLCLK}:g   \
                        | sed s:_TASKS_:${TASKS}:g       \
                        | sed s:_THRD_:${THRD}:g     \
                        | sed s:_RUND_:${RUNDIR}:g   \
                        | sed s:_SCHED_:${SCHEDULER}:g   >  gen_qsub

fi

if [ $SCHEDULER = 'loadleveler' ]; then
  llsubmit gen_ll 2>&1 | grep submitted > /dev/null
elif [ $SCHEDULER = 'moab' ]; then
  msub gen_msub > /dev/null
elif [ $SCHEDULER = 'pbs' ]; then
  rm -f $PATHRT/err $PATHRT/out
  qsub $PATHRT/gen_qsub > /dev/null
fi

echo "Test ${TEST_NR}" >> ${REGRESSIONTEST_LOG}
echo "Test ${TEST_NR}"
echo ${TEST_DESCR} >> ${REGRESSIONTEST_LOG}
echo ${TEST_DESCR}
(echo "GEN, ${TASKS} proc, ${THRD} thread";echo;echo)>> ${REGRESSIONTEST_LOG}
 echo "GEN, ${TASKS} proc, ${THRD} thread";echo;echo

# wait for the job to enter the queue
job_running=0
until [ $job_running -eq 1 ]
do
echo "TEST is waiting to enter the queue"
if [ $SCHEDULER = 'loadleveler' ]; then
job_running=`llq -u ${USER} -f %st %jn | grep ${JBNME} | wc -l`;sleep 5
elif [ $SCHEDULER = 'moab' -o $SCHEDULER = 'pbs' ]; then
job_running=`showq -u ${USER} -n | grep ${JBNME} | wc -l`;sleep 5
fi

done

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ]
do

if [ $SCHEDULER = 'loadleveler' ]; then

export status=`llq -u ${LOGIN} -f %st %jn | grep ${JBNME} | awk '{ print $1}'` ; export status=${status:--}

if   [ $status = 'I' ];  then echo $n "min. TEST ${TEST_NR} is waiting in a queue, Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST ${TEST_NR} is running,            Status: " $status
elif [ $status = 'ST' ]; then echo $n "min. TEST ${TEST_NR} is ready to run,       Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
else                          echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
fi

elif [ $SCHEDULER = 'moab' -o $SCHEDULER = 'pbs' ]; then

  status=`showq -u ${USER} -n | grep ${JBNME} | awk '{print $3}'` ; status=${status:--}
  if [ -f ${RUNDIR}/err ] ; then FnshHrs=`grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $6 }'` ; fi
  if   [ $status = 'Idle' ];       then echo $n "min. TEST ${TEST_NR} is waiting in a queue, Status: " $status
  elif [ $status = 'Running' ];    then echo $n "min. TEST ${TEST_NR} is running,            Status: " $status  
  elif [ $status = 'Starting' ];   then echo $n "min. TEST ${TEST_NR} is ready to run,       Status: " $status
  elif [ $status = 'Completed' ];  then echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
  else                                  echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
  fi
fi

sleep 60
if [ $SCHEDULER = 'loadleveler' ]; then
  job_running=`llq -u ${USER} -f %st %jn | grep ${JBNME} | wc -l`
elif [ $SCHEDULER = 'moab' -o $SCHEDULER = 'pbs' ]; then
  job_running=`showq -u ${USER} -n | grep ${JBNME} | wc -l`
fi
  (( n=n+1 ))
done

####################################################################################################
# Check results
####################################################################################################

(echo;echo;echo "Checking test ${TEST_NR} results ....")>> ${REGRESSIONTEST_LOG}
 echo;echo;echo "Checking test ${TEST_NR} results ...."

d=`grep 'Test run in the gen finalize routine' ${RUNDIR}/NEMS.out | wc -l`

  if [[ $d -eq 0 ]] ; then
   (echo " ......NOT OK" ; echo " Failed!   ")>> ${REGRESSIONTEST_LOG}
    echo " ......NOT OK" ; echo " Failed!   " ; exit 2
  fi

echo " Test ${TEST_NR} passed " >> ${REGRESSIONTEST_LOG}
echo " Test ${TEST_NR} passed "


(echo;echo;echo;)>> ${REGRESSIONTEST_LOG}
 echo;echo;echo;


sleep 4
clear;echo;echo

####################################################################################################
# End test
####################################################################################################

exit 0
