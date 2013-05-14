#!/bin/ksh
#set -xua

export GEFS_ENSEMBLE=${GEFS_ENSEMBLE:-0}
echo "GEFS_ENSEMBLE=" $GEFS_ENSEMBLE

mkdir -p ${RUNDIR}

export CDATE=${CDATE:-2012010100}
export NEMSIOIN=${NEMSIOIN:-.false.}
export SIGIOIN=${SIGIOIN:-.true.}
export NEMSIOOUT=${NEMSIOOUT:-.false.}
export SIGIOOUT=${SIGIOOUT:-.true.}
export MACHINE_ID=${MACHINE_ID:-ccs}
export SIGHDR=${SIGHDR:-/nwprod/exec/global_sighdr}
export MACHINE_ID=${MACHINE_ID:-WCOSS}
export SCHEDULER=${SCHEDULER:-lsf}

export IEMS=0
export ISOL=1
export ICO2=2
export IAER=111
if [ $GEFS_ENSEMBLE = 0 ] ; then

################################################################################
# For the stand alone GFS regression tests.
################################################################################

################################################################################
# Make configure and run files
################################################################################

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

 cat gfs_fcst_run.IN | sed s:_TASKS_:${TASKS}:g                   \
                     | sed s:_PE1_:${PE1}:g                       \
                     | sed s:_NEMSIOIN_:${NEMSIOIN}:g             \
                     | sed s:_NEMSIOOUT_:${NEMSIOOUT}:g           \
                     | sed s:_SIGIOIN_:${SIGIOIN}:g               \
                     | sed s:_SIGIOOUT_:${SIGIOOUT}:g             \
                     | sed s:_SFCIOOUT_:${SFCIOOUT}:g             \
                     | sed s:_WTPG_:${WTPG}:g                     \
                     | sed s:_WRTGP_:${WRTGP}:g                   \
                     | sed s:_wrtdopost_:${WRITE_DOPOST}:g        \
                     | sed s:_postgrbvs_:${POST_GRIBVERSION}:g    \
                     | sed s:_aer2post_:${GOCART_AER2POST}:g      \
                     | sed s:_WRTGP_:${WRTGP}:g                   \
                     | sed s:_THRDS_:${THRD}:g                    \
                     | sed s:_NSOUT_:${NSOUT}:g                   \
                     | sed s:_QUILT_:${QUILT}:g                   \
                     | sed s:_IAER_:${IAER}:g                     \
                     | sed s:_wave_:${wave}:g                     \
                     | sed s:_lm_:${lm}:g                         \
                     | sed s:_lsoil_:${lsoil}:g                   \
                     | sed s:_MEMBER_NAMES_:${MEMBER_NAMES}:g     \
                     | sed s:_CP2_:${CP2}:g                       \
                     | sed s:_RUNDIR_:${RUNDIR}:g                 \
                     | sed s:_PATHTR_:${PATHTR}:g                 \
                     | sed s:_FDFI_:${FDFI}:g                     \
                     | sed s:_FHRES_:${FHRES}:g                   \
                     | sed s:_REDUCEDGRID_:${REDUCEDGRID}:g       \
                     | sed s:_ADIAB_:${ADIAB}:g                   \
                     | sed s:_NSTFCST_:${NST_FCST}:g              \
                     | sed s:_GOCART_:${GOCART}:g                 \
                     | sed s:_TRACER_:${TRACER}:g                 \
                     | sed s:_SFCPRESSID_:${SFCPRESS_ID}:g        \
                     | sed s:_THERMODYNID_:${THERMODYN_ID}:g      \
                     | sed s:_IDVC_:${IDVC}:g                     \
                     | sed s:_NDSLFV_:${NDSLFV}:g                 \
                     | sed s:_SPECTRALLOOP_:${SPECTRALLOOP}:g     \
                     | sed s:_IDEA_:${IDEA}:g                     \
                     | sed s:_CDATE_:${CDATE}:g                   \
                     | sed s:_IEMS_:${IEMS}:g                     \
                     | sed s:_ISOL_:${ISOL}:g                     \
                     | sed s:_ICO2_:${ICO2}:g                     \
                     | sed s:_IAER_:${IAER}:g                     \
                     | sed s:_SIGHDR_:${SIGHDR}:g                 \
                     | sed s:_MACHINE_ID_:${MACHINE_ID}:g          \
                     | sed s:_RTPWD_:${RTPWD}:g                     \
                     | sed s:_SCHEDULER_:${SCHEDULER}:g            \
                     | sed s:_NDAYS_:${NDAYS}:g   >  gfs_fcst_run

 chmod 755 gfs_fcst_run

 cp gfs_fcst_run ${RUNDIR}

################################################################################
# Copy init files
################################################################################

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
    export IC_DIR=${IC_DIR:-${RTPWD}/GFS_DFI_REDUCEDGRID_HYB}
    cp $IC_DIR/gfsanl.$CDATE ${RUNDIR}/.
    cp $IC_DIR/sfcanl.$CDATE ${RUNDIR}/.

#   cp ${RTPWD}/GFS_DFI_REDUCEDGRID_HYB/gfsanl.$CDATE ${RUNDIR}/.
#   cp ${RTPWD}/GFS_DFI_REDUCEDGRID_HYB/sfcanl.$CDATE ${RUNDIR}/.
#to run gfs test
#   if [ "$rungfstest" = ".true." ]; then
#     cp /climate/save/wx20wa/esmf/nems/20120913/data/nemsio/gfsanl.$CDATE ${RUNDIR}/.
#     cp /climate/save/wx20wa/esmf/nems/20120913/data/nemsio/sfcanl.$CDATE ${RUNDIR}/.
#   fi

#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62hyb/gfsanl.$CDATE ${RUNDIR}/.
#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62hyb/sfnanl.$CDATE ${RUNDIR}/sfcanl.$CDATE


  elif [ $IDVC = 3 ] ; then
    export IC_DIR=${IC_DIR:-${RTPWD}/GFS_NODFI}
    cp $IC_DIR/gfsanl.$CDATE ${RUNDIR}/.
    cp $IC_DIR/sfnanl.$CDATE ${RUNDIR}/.

#   cp ${RTPWD}/GFS_NODFI/gfsanl.$CDATE ${RUNDIR}/.
#   cp ${RTPWD}/GFS_NODFI/sfcanl.$CDATE ${RUNDIR}/.

#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62/gfsanl.$CDATE ${RUNDIR}/.
#  cp /climate/noscrub/wx20wa/esmf/nems/IC/nemsio_new/t62/sfnanl.$CDATE ${RUNDIR}/sfcanl.$CDATE
  fi

#                     NO NEMSIO INPUT
#                     ---------------
 else 
   if [ "$IDEA" = ".true." ]; then
     cp ${RTPWD}/WAM_gh_l150/*anl*${CDATE} ${RUNDIR}/.
   else
     export dprefix2=${dprefix2:-""}
     if [ "$rungfstest" = ".true." ] ; then
       if [ $MACHINE_ID = wcoss ] ; then
         IC_DIR=${IC_DIR:-$dprefix2/global/noscrub/Shrinivas.Moorthi/data}
       elif [ $MACHINE_ID = ccs ] ; then
         IC_DIR=${IC_DIR:-$dprefix2/global/noscrub/wx23sm/data}
       elif [ $MACHINE_ID = zeus ] ; then
         IC_DIR=${IC_DIR:-$dprefix2/global/noscrub/Shrinivas.Moorthi/data}
       fi
       cp $IC_DIR/siganl.$CDATE ${RUNDIR}/.
       cp $IC_DIR/sfcanl.$CDATE ${RUNDIR}/.
     fi
   fi
 fi

else

################################################################################
# For the concurrency ensemble GEFS regression test.
################################################################################

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

################################################################################
# Submit test
################################################################################

JBNME=RT_${TEST_NR}_$$

if [ $SCHEDULER = 'loadleveler' ]; then

 export TPN=$((32/THRD))
 cat gfs_ll.IN       | sed s:_JBNME_:${JBNME}:g   \
                     | sed s:_CLASS_:${CLASS}:g   \
                     | sed s:_GROUP_:${GROUP}:g   \
                     | sed s:_ACCNR_:${ACCNR}:g   \
                     | sed s:_WLCLK_:${WLCLK}:g   \
                     | sed s:_TASKS_:${TASKS}:g   \
                     | sed s:_RUND_:${RUNDIR}:g   \
                     | sed s:_THRDS_:${THRD}:g    >  gfs_ll
 
elif [ $SCHEDULER = 'moab' ]; then

 export TPN=$((32/THRD))
 cat gfs_msub.IN     | sed s:_JBNME_:${JBNME}:g   \
                     | sed s:_WLCLK_:${WLCLK}:g   \
                     | sed s:_TPN_:${TPN}:g       \
                     | sed s:_THRD_:${THRD}:g     >  gfs_msub

elif [ $SCHEDULER = 'pbs' ]; then

 export TPN=$((12/THRD))
 cat gfs_qsub.IN     | sed s:_JBNME_:${JBNME}:g   \
                     | sed s:_ACCNR_:${ACCNR}:g   \
                     | sed s:_WLCLK_:${WLCLK}:g   \
                     | sed s:_TASKS_:${TASKS}:g   \
                     | sed s:_THRD_:${THRD}:g     \
                     | sed s:_RUND_:${RUNDIR}:g   \
                     | sed s:_SCHED_:${SCHEDULER}:g   >  gfs_qsub

elif [ $SCHEDULER = 'lsf' ]; then

 export TPN=$((16/THRD))
 cat gfs_bsub.IN     | sed s:_JBNME_:${JBNME}:g   \
                     | sed s:_CLASS_:${CLASS}:g   \
                     | sed s:_WLCLK_:${WLCLK}:g   \
                     | sed s:_TPN_:${TPN}:g       \
                     | sed s:_TASKS_:${TASKS}:g   \
                     | sed s:_RUND_:${RUNDIR}:g   \
                     | sed s:_THRDS_:${THRD}:g    \
                     | sed s:_CDATE_:${CDATE}:g   \
                     | sed s:_SCHED_:${SCHEDULER}:g   >  gfs_bsub
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
elif [ $SCHEDULER = 'lsf' ]; then
  bsub <$PATHRT/gfs_bsub > /dev/null
fi

echo "Test ${TEST_NR}" >> ${REGRESSIONTEST_LOG}
echo "Test ${TEST_NR}"
echo ${TEST_DESCR} >> ${REGRESSIONTEST_LOG}
echo ${TEST_DESCR}
(echo "GFS, ${TASKS} proc, ${THRD} thread";echo;echo)>> ${REGRESSIONTEST_LOG}
 echo "GFS, ${TASKS} proc, ${THRD} thread";echo;echo

# wait for the job to enter the queue
job_running=0
until [ $job_running -eq 1 ] ; do
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
until [ $job_running -eq 0 ] ; do

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
 elif [ $SCHEDULER = 'lsf' ] ; then

  status=`bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep " dev " | awk '{print $3}'` ; status=${status:--}
  if [ $status != '-' ] ; then FnshHrs=`bpeek -J ${JBNME} | grep Finished | tail -1 | awk '{ print $9 }'` ; fi
  if [ -f ${RUNDIR}/err ] ; then FnshHrs=`grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $9 }'` ; fi
  FnshHrs=${FnshHrs:-0}
  if   [ $status = 'PEND' ];  then echo $n "min. TEST ${TEST_NR} is waiting in a queue, Status: " $status
  elif [ $status = 'RUN'  ];  then echo $n "min. TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
  else                             echo $n "min. TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
  fi


 fi

 sleep 60
 if [ $SCHEDULER = 'loadleveler' ]; then
  job_running=`llq -u ${USER} -f %st %jn | grep ${JBNME} | wc -l`
 elif [ $SCHEDULER = 'moab' -o $SCHEDULER = 'pbs' ]; then
  job_running=`showq -u ${USER} -n | grep ${JBNME} | wc -l`
 elif [ $SCHEDULER = 'lsf' ] ; then
  job_running=`bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep " dev " | wc -l`
 fi
  (( n=n+1 ))
done

################################################################################
# Check results
################################################################################

(echo;echo;echo "Checking test ${TEST_NR} results ....")>> ${REGRESSIONTEST_LOG}
 echo;echo;echo "Checking test ${TEST_NR} results ...."

#
if [ ${CREATE_BASELINE} = false ]; then
#
# --- regression test comparison ----
#

  for i in ${LIST_FILES} ; do
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

 for i in ${LIST_FILES} ; do
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
