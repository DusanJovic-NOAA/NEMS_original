#!/bin/sh
#set -uaex

####################################################################################################
# Make configure and run files
####################################################################################################

JBNME=NEMS_RT_${TEST_NR}_$$

cat nmm_${GBRG}_run.IN  | sed s:_JBNME_:${JBNME}:g   \
                        | sed s:_TS_:${TS}:g         \
                        | sed s:_CPPCP_:${CPPCP}:g   \
                        | sed s:_POST_:${POST}:g     \
                        | sed s:_RTPWD_:${RTPWD}:g   \
                        | sed s:_SRCD_:${PATHTR}:g   \
                        | sed s:_RUND_:${RUNDIR}:g   >  nmm_run

cat nmm_${GBRG}_conf.IN | sed s:_INPES_:${INPES}:g   \
                        | sed s:_JNPES_:${JNPES}:g   \
                        | sed s:_WTPG_:${WTPG}:g     \
                        | sed s:_FCSTL_:${FCSTL}:g   \
                        | sed s:_NEMSI_:${NEMSI}:g   \
                        | sed s:_RSTRT_:${RSTRT}:g   \
                        | sed s:_gfsP_:${gfsP}:g     \
                        | sed s:_CONVC_:${CONVC}:g   \
                        | sed s:_MICRO_:${MICRO}:g   \
                        | sed s:_TURBL_:${TURBL}:g   \
                        | sed s:_PCPFLG_:${PCPFLG}:g \
                        | sed s:_WPREC_:${WPREC}:g   \
                        | sed s:_NCHILD_:${NCHILD}:g >  configure_file_01

if [ $SCHEDULER = 'loadleveler' ]; then

cat nmm_ll.IN           | sed s:_JBNME_:${JBNME}:g   \
                        | sed s:_CLASS_:${CLASS}:g   \
                        | sed s:_GROUP_:${GROUP}:g   \
                        | sed s:_ACCNR_:${ACCNR}:g   \
                        | sed s:_WLCLK_:${WLCLK}:g   \
                        | sed s:_TPN_:${TPN}:g       \
                        | sed s:_THRD_:${THRD}:g     \
                        | sed s:_AFFNP_:${AFFN}:g    \
                        | sed s:_NODE_:${NODE}:g     >  nmm_ll

elif [ $SCHEDULER = 'moab' ]; then

cat nmm_msub.IN         | sed s:_JBNME_:${JBNME}:g   \
                        | sed s:_WLCLK_:${WLCLK}:g   \
                        | sed s:_TPN_:${TPN}:g       \
                        | sed s:_THRD_:${THRD}:g     >  nmm_msub

elif [ $SCHEDULER = 'pbs' ]; then

cat nmm_qsub.IN         | sed s:_JBNME_:${JBNME}:g   \
                        | sed s:_WLCLK_:${WLCLK}:g   \
                        | sed s:_TPN_:${TPN}:g       \
                        | sed s:_THRD_:${THRD}:g     \
                        | sed s:_RUND_:${RUNDIR}:g   >  nmm_qsub

fi

if [ ${GBRG} = nests ]; then
  cat ${RTPWD}/NMMB_nests/configure_file_02.IN | sed s:_RSTRT_:${RSTRT}:g > configure_file_02
  cat ${RTPWD}/NMMB_nests/configure_file_03.IN | sed s:_RSTRT_:${RSTRT}:g > configure_file_03
  cat ${RTPWD}/NMMB_nests/configure_file_04.IN | sed s:_RSTRT_:${RSTRT}:g > configure_file_04
fi

if [ ${GBRG} = mnests ]; then
  rm -f configure_file_02 configure_file_03 configure_file_04
  cat ${RTPWD}/NMMB_mvg_nests/configure_file_02.IN | sed s:_RSTRT_:${RSTRT}:g > configure_file_02
  cat ${RTPWD}/NMMB_mvg_nests/configure_file_03.IN | sed s:_RSTRT_:${RSTRT}:g > configure_file_03
  cat ${RTPWD}/NMMB_mvg_nests/configure_file_04.IN | sed s:_RSTRT_:${RSTRT}:g > configure_file_04
fi

if [ ${GBRG} = fltr ]; then
  rm -f configure_file_02 configure_file_03 configure_file_04
  cat ${RTPWD}/NMMB_reg_filt/configure_file_02.IN | sed s:_RSTRT_:${RSTRT}:g > configure_file_02
  cat ${RTPWD}/NMMB_reg_filt/configure_file_03.IN | sed s:_RSTRT_:${RSTRT}:g > configure_file_03
fi

####################################################################################################
# Submit test
####################################################################################################

sh ./nmm_run

echo "Test ${TEST_NR}" >> RegressionTests.log
echo "Test ${TEST_NR}"
echo ${TEST_DESCR} >> RegressionTests.log
echo ${TEST_DESCR}


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

  status=`llq -u ${USER} -f %st %jn | grep ${JBNME} | awk '{ print $1}'` ; status=${status:--}
  if [ -f ${RUNDIR}/err ] ; then FnshHrs=`grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $7 }'` ; fi
  FnshHrs=${FnshHrs:-0}
  if   [ $status = 'I' ];  then echo $n "min. TEST ${TEST_NR} is waiting in a queue, Status: " $status
  elif [ $status = 'R' ];  then echo $n "min. TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
  elif [ $status = 'ST' ]; then echo $n "min. TEST ${TEST_NR} is ready to run,       Status: " $status
  elif [ $status = 'C' ];  then echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
  else                          echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
  fi

elif [ $SCHEDULER = 'moab' -o $SCHEDULER = 'pbs' ]; then

  status=`showq -u ${USER} -n | grep ${JBNME} | awk '{print $3}'` ; status=${status:--}
  if [ -f ${RUNDIR}/err ] ; then FnshHrs=`grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $6 }'` ; fi
  FnshHrs=${FnshHrs:-0}
  if   [ $status = 'Idle' ];       then echo $n "min. TEST ${TEST_NR} is waiting in a queue, Status: " $status
  elif [ $status = 'Running' ];    then echo $n "min. TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
  elif [ $status = 'Starting' ];   then echo $n "min. TEST ${TEST_NR} is ready to run,       Status: " $status  ", Finished " $FnshHrs "hours"
  elif [ $status = 'Completed' ];  then echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status
  else                                  echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
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

 echo;echo;echo "Moving set ${TEST_NR} files ...."

for i in ${LIST_FILES}
do
  printf %s " Moving " $i "....."
  if [ -f ${RUNDIR}/$i ] ; then
    cp ${RUNDIR}/${i} /stmp/${USER}/REGRESSION_TEST/${CNTL_DIR}/${i}
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
