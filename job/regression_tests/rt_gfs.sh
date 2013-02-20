#!/bin/ksh
set -ua

export GEFS_ENSEMBLE=${GEFS_ENSEMBLE:-0}
echo "GEFS_ENSEMBLE=" $GEFS_ENSEMBLE

mkdir -p ${RUNDIR}

if [ $GEFS_ENSEMBLE = 0 ] ; then

####################################################################################################
# For the stand alone GFS regression tests.
####################################################################################################

####################################################################################################
# Make configure and run files
####################################################################################################

## determine GOCART and TRACER from gocart_aerosol and passive_tracer 
export gocart_aerosol=${gocart_aerosol:-NO}
export passive_tracer=${passive_tracer:-NO}
if [ $gocart_aerosol = 'YES' ]; then
 export GOCART=1 
else
 export GOCART=0 
fi
if  [ $passive_tracer = 'YES' ]; then
 export TRACER=.true.
else
 export TRACER=.false.
fi
##

JBNME=NEMS_RT_${TEST_NR}_$$

cd $PATHRT

cat gfs_fcst_run.IN | sed s:_TASKS_:${TASKS}:g   \
                    | sed s:_PE1_:${PE1}:g       \
                    | sed s:_NEMSIOIN_:${NEMSIOIN}:g      \
                    | sed s:_NEMSIOOUT_:${NEMSIOOUT}:g      \
                    | sed s:_SIGIOOUT_:${SIGIOOUT}:g      \
                    | sed s:_SFCIOOUT_:${SFCIOOUT}:g      \
                    | sed s:_WPG_:${WTPG}:g      \
                    | sed s:_WRTGP_:${WRTGP}:g   \
                    | sed s:_wrtdopost_:${WRITE_DOPOST}:g   \
                    | sed s:_postgrbvs_:${POST_GRIBVERSION}:g   \
                    | sed s:_aer2post_:${GOCART_AER2POST}:g   \
                    | sed s:_WRTGP_:${WRTGP}:g   \
                    | sed s:_THRDS_:${THRD}:g    \
                    | sed s:_NSOUT_:${NSOUT}:g   \
                    | sed s:_QUILT_:${QUILT}:g   \
                    | sed s:_IAER_:${IAER}:g       \
                    | sed s:_wavecoef_:${wavecoef}:g     \
                    | sed s:_wavegrid_:${wavegrid}:g     \
                    | sed s:_lm_:${lm}:g         \
                    | sed s:_lsoil_:${lsoil}:g   \
                    | sed s:_MEMBER_NAMES_:${MEMBER_NAMES}:g   \
                    | sed s:_CP2_:${CP2}:g       \
                    | sed s:_RUNDIR_:${RUNDIR}:g \
                    | sed s:_PATHTR_:${PATHTR}:g \
                    | sed s:_FDFI_:${FDFI}:g \
                    | sed s:_FHRES_:${FHRES}:g \
                    | sed s:_REDUCEDGRID_:${REDUCEDGRID}:g \
                    | sed s:_ADIAB_:${ADIAB}:g \
                    | sed s:_NSTFCST_:${NST_FCST}:g \
                    | sed s:_GOCART_:${GOCART}:g \
                    | sed s:_TRACER_:${TRACER}:g \
                    | sed s:_SFCPRESSID_:${SFCPRESS_ID}:g \
                    | sed s:_THERMODYNID_:${THERMODYN_ID}:g \
                    | sed s:_IDVC_:${IDVC}:g \
                    | sed s:_NDSLFV_:${NDSLFV}:g \
                    | sed s:_SPECTRALLOOP_:${SPECTRALLOOP}:g \
                    | sed s:_IDEA_:${IDEA}:g \
                    | sed s:_CDATE_:${CDATE}:g \
                    | sed s:_NDAYS_:${NDAYS}:g   >  gfs_fcst_run

chmod 755 gfs_fcst_run

cp gfs_fcst_run ${RUNDIR}

####################################################################################################
# Copy init files
####################################################################################################

cp atmos.configure_gfs ${RUNDIR}/atmos.configure
cp ocean.configure ${RUNDIR}/ocean.configure
cp MAPL.rc ${RUNDIR}/MAPL.rc
cp Chem_Registry.rc ${RUNDIR}/Chem_Registry.rc

if [ $GOCART = 1 ] ; then
 export EXTDIR=/global/save/wx23lu/NEMS/fix
 export RCSDIR=/global/save/wx23lu/NEMS/Chem_Registry
 cp ${RCSDIR}/*.rc ${RUNDIR}/.
 cp -r  ${EXTDIR}/ExtData ${RUNDIR}/.
fi

if [ "$NEMSIOIN" = ".true." ]; then
  if [ $IDVC = 2 ] ; then
    cp ${RTPWD}/GFS_DFI_REDUCEDGRID_HYB/gfsanl.2010010100 ${RUNDIR}/.
    cp ${RTPWD}/GFS_DFI_REDUCEDGRID_HYB/sfcanl.2010010100 ${RUNDIR}/.
#to run gfs test
    if [ "$rungfstest" = ".true." ]; then
      cp /climate/save/wx20wa/esmf/nems/20120913/data/nemsio/gfsanl.2012010100 ${RUNDIR}/.
      cp /climate/save/wx20wa/esmf/nems/20120913/data/nemsio/sfcanl.2012010100 ${RUNDIR}/.
      export CDATE=2012010100
    fi

#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62hyb/gfsanl.2010010100 ${RUNDIR}/.
#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62hyb/sfnanl.2010010100 ${RUNDIR}/sfcanl.2010010100


  elif [ $IDVC = 3 ] ; then
    cp ${RTPWD}/GFS_NODFI/gfsanl.2010010100 ${RUNDIR}/.
    cp ${RTPWD}/GFS_NODFI/sfcanl.2010010100 ${RUNDIR}/.

#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62/gfsanl.2010010100 ${RUNDIR}/.
#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62/sfnanl.2010010100 ${RUNDIR}/sfcanl.2010010100
  fi

#no nemsio input
else
   if [ "$rungfstest" = ".true." ]; then
     cp /climate/save/wx20wa/esmf/nems/20120913/data/siganl.2012010100 ${RUNDIR}/.
     cp /climate/save/wx20wa/esmf/nems/20120913/data/sfcanl.2012010100 ${RUNDIR}/.
     export CDATE=2012010100
   fi
fi

if [ "$IDEA" = ".true." ]; then
  cp /climate/save/wx20wa/wam/Fei_Wu/data/*anl*${CDATE} ${RUNDIR}/.
fi

#ensembl
else


####################################################################################################
# For the concurrency ensemble GEFS regression test.
####################################################################################################

cd $PATHRT

 cp ${RTPWD}/GEFS_data_2008082500/* $RUNDIR
# cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/GEFS_data_2008082500/gfsanl* $RUNDIR
# cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/GEFS_data_2008082500/sfcanl* $RUNDIR
cp $PATHRT/gfs_configfile_190 $RUNDIR/configure_file

cat gfs_fcst_run_GEFS.IN \
                    | sed s:_SRCDIR_:${PATHTR}:g \
                    | sed s:_NDSLFV_:${NDSLFV}:g \
                    | sed s:_NEMSIOIN_:${NEMSIOIN}:g \
                    | sed s:_IDEA_:${IDEA}:g \
                    | sed s:_RUNDIR_:${RUNDIR}:g > gfs_fcst_run


cp gfs_fcst_run ${RUNDIR}
chmod +x ${RUNDIR}/gfs_fcst_run
cp Chem_Registry.rc ${RUNDIR}/Chem_Registry.rc
cp atmos.configure_gfs ${RUNDIR}/atmos.configure
cp ocean.configure ${RUNDIR}/ocean.configure

fi

####################################################################################################
# Submit test
####################################################################################################

JBNME=RT_${TEST_NR}_$$

if [ $SCHEDULER = 'loadleveler' ]; then

cat gfs_ll.IN       | sed s:_JBNME_:${JBNME}:g   \
                    | sed s:_CLASS_:${CLASS}:g   \
                    | sed s:_GROUP_:${GROUP}:g   \
                    | sed s:_ACCNR_:${ACCNR}:g   \
                    | sed s:_WLCLK_:${WLCLK}:g   \
                    | sed s:_TASKS_:${TASKS}:g   \
                    | sed s:_THRDS_:${THRD}:g    >  gfs_ll

elif [ $SCHEDULER = 'moab' ]; then

cat gfs_msub.IN         | sed s:_JBNME_:${JBNME}:g   \
                        | sed s:_WLCLK_:${WLCLK}:g   \
                        | sed s:_TPN_:${TPN}:g       \
                        | sed s:_THRD_:${THRD}:g     >  gfs_msub

elif [ $SCHEDULER = 'pbs' ]; then

cat gfs_qsub.IN         | sed s:_JBNME_:${JBNME}:g   \
                        | sed s:_ACCNR_:${ACCNR}:g   \
                        | sed s:_WLCLK_:${WLCLK}:g   \
                        | sed s:_TASKS_:${TASKS}:g       \
                        | sed s:_THRD_:${THRD}:g     \
                        | sed s:_RUND_:${RUNDIR}:g   \
                        | sed s:_SCHED_:${SCHEDULER}:g   >  gfs_qsub

fi

  cp ../exglobal_fcst.sh.sms_nems ${RUNDIR}

export RUNDIR=${RUNDIR}

cd $PATHRT

if [ $SCHEDULER = 'loadleveler' ]; then
  llsubmit gfs_ll 2>&1 | grep submitted > /dev/null
elif [ $SCHEDULER = 'moab' ]; then
  msub gfs_msub > /dev/null
elif [ $SCHEDULER = 'pbs' ]; then
  rm -f $PATHRT/err $PATHRT/out
  qsub $PATHRT/gfs_qsub > /dev/null
fi

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
if [ $SCHEDULER = 'loadleveler' ]; then
  job_running=`llq -u ${USER} -f %st %jn | grep ${JBNME} | wc -l`;sleep 5
elif [ $SCHEDULER = 'moab' ]; then
  job_running=`showq -u ${USER} -n | grep ${JBNME} | wc -l`;sleep 5
elif [ $SCHEDULER = 'pbs' ]; then
  job_running=`qstat -u ${USER} -n | grep ${JBNME} | wc -l`;sleep 5
elif [ $SCHEDULER = 'lsf' ]; then
  job_running=`bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep " dev " | wc -l`;sleep 5
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

echo " Test ${TEST_NR} passed " >> ${REGRESSIONTEST_LOG}
echo " Test ${TEST_NR} passed "

sleep 4
clear;echo;echo

####################################################################################################
# End test
####################################################################################################

exit 0
