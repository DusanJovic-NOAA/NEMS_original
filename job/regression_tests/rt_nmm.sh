#!/bin/ksh
set -ua

####################################################################################################
# Make configure and run files
####################################################################################################

JBNME=NEMS_RT_${TEST_NR}_$$

cat nmm_${GBRG}_ll.IN   | sed s:_JBNME_:${JBNME}:g   \
                        | sed s:_CLASS_:${CLASS}:g   \
                        | sed s:_ACCNR_:${ACCNR}:g   \
                        | sed s:_TPN_:${TPN}:g       \
                        | sed s:_THRD_:${THRD}:g     \
                        | sed s:_GS_:${GS}:g         \
                        | sed s:_CPPCP_:${CPPCP}:g   \
                        | sed s:_RTPWD_:${RTPWD}:g   \
                        | sed s:_SRCD_:${PATHTR}:g   \
                        | sed s:_RUND_:${RUNDIR}:g   >  nmm_ll

cat nmm_${GBRG}_conf.IN | sed s:_INPES_:${INPES}:g   \
                        | sed s:_JNPES_:${JNPES}:g   \
                        | sed s:_WTPG_:${WTPG}:g     \
                        | sed s:_FCSTL_:${FCSTL}:g   \
                        | sed s:_NEMSI_:${NEMSI}:g   \
                        | sed s:_RSTRT_:${RSTRT}:g   \
                        | sed s:_gfsP_:${gfsP}:g     \
                        | sed s:_PCPFLG_:${PCPFLG}:g \
                        | sed s:_WPREC_:${WPREC}:g   \
                        | sed s:_RGS_:${RGS}:g       \
                        | sed s:_WGS_:${WGS}:g       \
                        | sed s:_NCHILD_:${NCHILD}:g >  configure_file

####################################################################################################
# Submit test
####################################################################################################

llsubmit nmm_ll 2>&1 | grep submitted > /dev/null

echo "Test ${TEST_NR}" >> RegressionTests.log
echo "Test ${TEST_NR}"
echo ${TEST_DESCR} >> RegressionTests.log
echo ${TEST_DESCR}


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

if [ -f $PATHRT/err ] ; then FnshHrs=`grep Finished $PATHRT/err | tail -1 | awk '{ print $7 }'` ; fi
export FnshHrs=${FnshHrs:-0}

if   [ $status = 'I' ];  then echo $n "min. TEST ${TEST_NR} is waiting in a queue, Status: " $status
elif [ $status = 'R' ];  then echo $n "min. TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
elif [ $status = 'ST' ]; then echo $n "min. TEST ${TEST_NR} is ready to run,       Status: " $status
elif [ $status = 'C' ];  then echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
else                          echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
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

#
     if [ ${CREATE_BASELINE} = false ]; then
#
# --- regression test comparison ----
#

for i in ${LIST_FILES}
do
printf %s " Comparing " $i "....." >> RegressionTests.log
printf %s " Comparing " $i "....."

if [ -f ${RUNDIR}/$i ] ; then

  d=`cmp ${RTPWD}/${CNTL_DIR}/$i ${RUNDIR}/$i | wc -l`

  if [[ $d -ne 0 ]] ; then
   (echo " ......NOT OK" ; echo ; echo "   $i differ!   ")>> RegressionTests.log
    echo " ......NOT OK" ; echo ; echo "   $i differ!   " ; exit 2
  fi

  echo "....OK" >> RegressionTests.log
  echo "....OK"

else

  echo "Missing " ${RUNDIR}/$i " output file" >> RegressionTests.log
  echo "Missing " ${RUNDIR}/$i " output file"
 (echo;echo " Test ${TEST_NR} failed ")>> RegressionTests.log
  echo;echo " Test ${TEST_NR} failed "
  exit 2

fi

done

#
     else
#
# --- create baselines
#

#mkdir -p /stmp/${LOGIN}/REGRESSION_TEST/${CNTL_DIR}

 echo;echo;echo "Moving set ${TEST_NR} files ...."

for i in ${LIST_FILES}
do
  printf %s " Moving " $i "....."
  if [ -f ${RUNDIR}/$i ] ; then
    mv ${RUNDIR}/${i} /stmp/${LOGIN}/REGRESSION_TEST/${CNTL_DIR}/${i}
  else
    echo "Missing " ${RUNDIR}/$i " output file"
    echo;echo " Set ${TEST_NR} failed "
    exit 2
  fi
done

# ---
     fi
# ---

echo " Test ${TEST_NR} passed " >> RegressionTests.log
(echo;echo;echo)                >> RegressionTests.log
echo " Test ${TEST_NR} passed "

sleep 4
clear;echo;echo

####################################################################################################
# End test
####################################################################################################

exit 0
