#!/bin/ksh
#set -eu
#set -x

export GEFS_ENSEMBLE=${GEFS_ENSEMBLE:-0}
echo "GEFS_ENSEMBLE=" $GEFS_ENSEMBLE

source ./atparse.auto

mkdir -p ${RUNDIR}

export CDATE=${CDATE:-2012010100}
export NEMSIOIN=${NEMSIOIN:-.false.}
export SIGIOIN=${SIGIOIN:-.true.}
export SFCIOOUT=${SFCIOOUT:-.true.}
export NEMSIOOUT=${NEMSIOOUT:-.false.}
export SIGIOOUT=${SIGIOOUT:-.true.}
export MACHINE_ID=${MACHINE_ID:-wcoss}
if [ "$NEMSIOIN" = ".true." ]; then
  if [ $MACHINE_ID = wcoss ]; then
    export SIGHDR=/nwprod/exec/global_sighdr
  elif [ $MACHINE_ID = theia ]; then
    export SIGHDR=/scratch4/NCEPDEV/nems/save/Jun.Wang/nems/util/nemsio_get
  fi
else
  export SIGHDR=${SIGHDR:-/nwprod/exec/global_sighdr}
fi
export fcst_begin=${fcst_begin:-YES}

export SCHEDULER=${SCHEDULER:-lsf}
export SHOWQ=${SHOWQ:-/opt/moab/default/bin/showq}
export MSUB=${MSUB:-/opt/moab/default/bin/msub}
export pex=${pex:-1}

export IEMS=0
export ISOL=1
export ICO2=2
export IAER=111
export fcyc=24
export IALB=${IALB:-0}
export isot=${isot:-0}
export ivegsrc=${ivegsrc:-2}
export FHOUT_HF=${FHOUT_HF:-1}
export FHMAX_HF=${FHMAX_HF:-0}

if [ $GEFS_ENSEMBLE = 0 ] ; then

################################################################################
# For the stand alone GFS regression tests.
################################################################################

################################################################################
# Make configure and run files
################################################################################

## determine GOCART and TRACER from gocart_aerosol and passive_tracer 
 export gocart_aerosol=${gocart_aerosol:-NO}
 export passive_tracer=${passive_tracer:-YES}
 if [ $gocart_aerosol = 'YES' ] ; then
  export GOCART=1 
 else
  export GOCART=0 
 fi
 if  [ $passive_tracer = 'YES' ] ; then
  export TRACER=.true.
 else
  export TRACER=.false.
 fi
##

 export WAM_IPE_COUPLING=${WAM_IPE_COUPLING:-.false.}
 export HEIGHT_DEPENDENT_G=${HEIGHT_DEPENDENT_G:-.false.}
 export F107_KP_SKIP_SIZE=${F107_KP_SKIP_SIZE:-0}

 cd $PATHRT

 (   # Rename misnamed variables in gfs_fcst_run.IN's @[...] blocks:
     wrtdopost="$WRITE_DOPOST" ;
     postgrbvs="$POST_GRIBVERSION" ;
     aer2post="$GOCART_AER2POST" ;
     THRDS="$THRD" ;
     ISOT="$isot" ;
     IVEGSRC="$ivegsrc" ;
     FCST_BEGIN="$fcst_begin" ;
     NSTFCST="$NST_FCST" ;
     NSTSPINUP="$NST_SPINUP" ;
     NSTREV="$NST_RESERVED" ;
     SFCPRESSID="$SFCPRESS_ID" ;
     THERMODYNID="$THERMODYN_ID" ;
     SKEBVFILT="${SKEB_VFILT}" ;
     SPPTLOGIT="${SPPT_LOGIT}" ;
     SPPTSFCLIMIT="${SPPT_SFCLIMIT}" ;
     SPPTTAU="${SPPT_TAU}" ;
     SPPTLSCALE="${SPPT_LSCALE}" ;
     SHUMTAU="${SHUM_TAU}" ;
     SHUMLSCALE="${SHUM_LSCALE}" ;
     SKEBTAU="${SKEB_TAU}" ;
     SKEBLSCALE="${SKEB_LSCALE}" ;
     ISEEDSPPT="${ISEED_SPPT}" ;
     ISEEDSHUM="${ISEED_SHUM}" ;
     ISEEDSKEB="${ISEED_SKEB}" ;
     FHOUTHF="$FHOUT_HF" ;
     FHMAXHF="$FHMAX_HF" ;
     # Generate the script:
     cat gfs_fcst_run.IN | atparse > gfs_fcst_run )
 chmod 755 gfs_fcst_run

 # | sed s:_LDFI_SPECT_:${LDFI_SPECT}:g               \

 cp gfs_fcst_run ${RUNDIR}

if [ ${nems_configure}"x" == "x" ]; then
  nems_configure=atm_nostep
  atm_model=gsm
fi
 cat nems.configure.${nems_configure}.IN | atparse >  nems.configure
                         
 cp nems.configure ${RUNDIR}

################################################################################
# Copy init files
################################################################################

 cat atmos.configure_gfs | atparse >  atmos.configure
 cp atmos.configure ${RUNDIR}/atmos.configure
 cp MAPL.rc ${RUNDIR}/MAPL.rc
 cp Chem_Registry.rc ${RUNDIR}/Chem_Registry.rc

 if [ $GOCART = 1 ] ; then
  if [ $SCHEDULER = 'loadleveler' ]; then
    export EXTDIR=/global/save/wx23lu/NEMS/fix
    export RCSDIR=/global/save/wx23lu/NEMS/Chem_Registry
    cp -r ${EXTDIR}/ExtData ${RUNDIR}/.
  elif [ $SCHEDULER = 'pbs' ]; then
    export EXTDIR=@[RTPWD]/data_GOCART
    export RCSDIR=@[RTPWD]/data_GOCART
    cp -r ${EXTDIR}/ngac_fix ${RUNDIR}/.
  elif [ $SCHEDULER = 'lsf' ]; then
    export EXTDIR=/nwprod/ngac.v1.0.0/fix
    export RCSDIR=@[RTPWD]/data_GOCART
    cp -r ${EXTDIR}/ngac_fix ${RUNDIR}/.
  fi
 fi

 if [ "$NEMSIOIN" = ".true." ]; then
  export IC_DIR=${IC_DIR:-${RTPWD}/$DATAICDIR}
  if [ $MACHINE_ID = wcoss ] ; then
     export nemsioget=${nemsioget:-/nwprod/ngac.v1.0.0/exec/nemsio_get}
  elif [ $MACHINE_ID = theia ] ; then
     export nemsioget=${nemsioget:-/scratch3/NCEPDEV/nems/save/Jun.Wang/nems/util/nemsio_get}
  fi
  export SIGHDR=$nemsioget
  if [ "$IDEA" = ".true." ]; then
                                     
    cp $IC_DIR/*anl*${CDATE} ${RUNDIR}/.                                                      
                                       
                                                                                                               
      
  else
    if [ $fcst_begin = YES ]; then
      cp $IC_DIR/gfsanl.$CDATE $RUNDIR
      cp $IC_DIR/sfnanl.$CDATE $RUNDIR
      if [ $NST_FCST -gt 0 ] ; then
        cp $IC_DIR/nsnanl.$CDATE $RUNDIR
      fi
    fi
  fi

# These gfsanl and sfnanl data were copy from Moorthi's directory at
# /global/noscrub/Shrinivas.Moorthi/data on Surge machine. Weiyu.
#-------------------------------------------------------------------
# cp ${RTPWD}/GFS_SLG_NEMSIO_READ/gfsanl.$CDATE $RUNDIR
# cp ${RTPWD}/GFS_SLG_NEMSIO_READ/sfnanl.$CDATE $RUNDIR
#  cp $IC_DIR/gfsanl.$CDATE $RUNDIR
#  cp $IC_DIR/sfnanl.$CDATE $RUNDIR

#                     NO NEMSIO INPUT
#                     ---------------
 else 
   if [ "$IDEA" = ".true." ]; then
     cp ${RTPWD}/WAM_gh_l150/*anl*${CDATE} ${RUNDIR}/.
   else
     export dprefix=${dprefix:-""}
     if [ "$rungfstest" = ".true." ] ; then
       if [ $MACHINE_ID = wcoss ] ; then
         IC_DIR=${IC_DIR:-$dprefix/global/noscrub/Shrinivas.Moorthi/data}
       elif [ $MACHINE_ID = theia ] ; then
         IC_DIR=${IC_DIR:-$dprefix/global/noscrub/Shrinivas.Moorthi/data}
       fi
          cp $IC_DIR/siganl.$CDATE ${RUNDIR}/.
          cp $IC_DIR/sfcanl.$CDATE ${RUNDIR}/.
          if [ $NST_FCST -gt 0 ] ; then
            cp $IC_DIR/nstanl.$CDATE ${RUNDIR}/.
          fi
     fi
   fi
 fi

else

################################################################################
# For the concurrency ensemble GEFS regression test.
################################################################################

 cd $PATHRT

 cp ${RTPWD}/GEFS_data_2008082500/* $RUNDIR

 atparse SRCDIR="$PATHTR" < gfs_fcst_run_GEFS.IN > gfs_fcst_run
 
 cp gfs_fcst_run ${RUNDIR}
 chmod +x ${RUNDIR}/gfs_fcst_run
 cp Chem_Registry.rc ${RUNDIR}/Chem_Registry.rc
 cp atmos.configure_gfs ${RUNDIR}/atmos.configure

fi

################################################################################
# Submit test
################################################################################

if [ $SCHEDULER = 'moab' ]; then

 export TPN=$((32/THRD))
 cat gfs_msub.IN | atparse RUND=$RUNDIR >  gfs_msub


elif [ $SCHEDULER = 'pbs' ]; then

 export TPN=$((24/THRD))
 export QUEUE=${QUEUE:-batch}
 cat gfs_qsub.IN | atparse RUND="$RUNDIR" SCHED="$SCHEDULER" > gfs_qsub

elif [ $SCHEDULER = 'lsf' ]; then

 if [ $MACHINE_ID = wcoss ] ; then
   export pex=${pex:-1}
   export QUEUE=${QUEUE:-dev}
   if [ $pex = 2 ] ; then
     export TPN=${TPN:-$((24/THRD))}
     export QUEUE=${QUEUE:-dev$pex}
   else
     export TPN=${TPN:-$((16/THRD))}
   fi
 fi
 cat gfs_bsub.IN | atparse RUND=$RUNDIR THRDS=$THRD SCHED=$SCHEDULER \
     > gfs_bsub
fi

cp exglobal_fcst_nems.sh $RUNDIR

export RUNDIR=$RUNDIR

cd $PATHRT

if [ $SCHEDULER = 'moab' ]; then
  $MSUB gfs_msub
elif [ $SCHEDULER = 'pbs' ]; then
  rm -f $PATHRT/err $PATHRT/out
  qsub $PATHRT/gfs_qsub
elif [ $SCHEDULER = 'lsf' ]; then
  bsub <$PATHRT/gfs_bsub
fi

echo "Test ${TEST_NR}" >> ${REGRESSIONTEST_LOG}
echo "Test ${TEST_NR}"
echo ${TEST_DESCR} >> ${REGRESSIONTEST_LOG}
echo ${TEST_DESCR}
(echo "GFS, ${TASKS} proc, ${THRD} thread")>> ${REGRESSIONTEST_LOG}
 echo "GFS, ${TASKS} proc, ${THRD} thread"

# wait for the job to enter the queue
job_running=0

until [ $job_running -eq 1 ] ; do
 echo "TEST is waiting to enter the queue"
 if [ $SCHEDULER = 'moab' ]; then
  job_running=`$SHOWQ -u ${USER} -n | grep ${JBNME} | wc -l`;sleep 5
 elif [ $SCHEDULER = 'pbs' ]; then
  job_running=`qstat -u ${USER} -n | grep ${JBNME} | wc -l`;sleep 5
 elif [ $SCHEDULER = 'lsf' ]; then
  job_running=`bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep ${QUEUE} | wc -l`;sleep 5
 fi
done

job_running=1

# wait for the job to finish and compare results
n=1
until [ $job_running -eq 0 ] ; do

 sleep 30

 if [ $SCHEDULER = 'moab' ]; then

  status=`$SHOWQ -u ${USER} -n | grep ${JBNME} | awk '{print $3}'` ; status=${status:--}
  if [ -f ${RUNDIR}/err ] ; then FnshHrs=`grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $6 }'` ; fi
  FnshHrs=${FnshHrs:-0}
  if   [ $status = 'Idle' ];       then echo "$n/2min TEST ${TEST_NR} is waiting in a queue, Status: " $status
  elif [ $status = 'Running' ];    then echo "$n/2min TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
  elif [ $status = 'Starting' ];   then echo "$n/2min TEST ${TEST_NR} is ready to run,       Status: " $status  ", Finished " $FnshHrs "hours"
  elif [ $status = 'Completed' ];  then echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status
  else                                  echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
  fi

 elif [ $SCHEDULER = 'pbs' ]; then

   status=`qstat -u ${USER} -n | grep ${JBNME} | awk '{print $"10"}'` ; status=${status:--}
   if [ -f ${RUNDIR}/err ] ; then FnshHrs=`grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $10 }'` ; fi
   FnshHrs=${FnshHrs:-0}
   if   [ $status = 'Q' ];  then echo "$n/2min TEST ${TEST_NR} is waiting in a queue, Status: " $status
   elif [ $status = 'H' ];  then echo "$n/2min TEST ${TEST_NR} is held in a queue,    Status: " $status
   elif [ $status = 'R' ];  then echo "$n/2min TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
   elif [ $status = 'E' ];  then echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status ; job_running=0
   elif [ $status = 'C' ];  then echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status ; job_running=0
   else                          echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
   fi

 elif [ $SCHEDULER = 'lsf' ] ; then

  status=`bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep ${QUEUE} | awk '{print $3}'` ; status=${status:--}
#  if [ $status != '-' ] ; then FnshHrs=`bpeek -J ${JBNME} | grep Finished | tail -1 | awk '{ print $9 }'` ; fi
  if [ -f ${RUNDIR}/err ] ; then FnshHrs=`grep Finished ${RUNDIR}/err | tail -1 | awk '{ print $9 }'` ; fi
  FnshHrs=${FnshHrs:-0}
  if   [ $status = 'PEND' ];  then echo "$n/2min TEST ${TEST_NR} is waiting in a queue, Status: " $status
  elif [ $status = 'RUN'  ];  then echo "$n/2min TEST ${TEST_NR} is running,            Status: " $status  ", Finished " $FnshHrs "hours"
  else                             echo "$n/2min TEST ${TEST_NR} is finished,           Status: " $status  ", Finished " $FnshHrs "hours"
  fi

 fi

 if [ $SCHEDULER = 'moab' ]; then
  job_running=`$SHOWQ -u ${USER} -n | grep ${JBNME} | wc -l`
 elif [ $SCHEDULER = 'lsf' ] ; then
  job_running=`bjobs -u ${USER} -J ${JBNME} 2>/dev/null | grep ${QUEUE} | wc -l`
 fi
  (( n=n+1 ))
done

################################################################################
# Check results
################################################################################

test_status='PASS'

# Give 10 seconds for data to show up on file system
sleep 10

(echo;echo;echo "baseline dir = ${RTPWD}/${CNTL_DIR}")  >> ${REGRESSIONTEST_LOG}
           echo "working dir  = ${RUNDIR}"              >> ${REGRESSIONTEST_LOG}
           echo "Checking test ${TEST_NR} results ...." >> ${REGRESSIONTEST_LOG}
(echo;echo;echo "baseline dir = ${RTPWD}/${CNTL_DIR}")
           echo "working dir  = ${RUNDIR}"
           echo "Checking test ${TEST_NR} results ...."

#
if [ ${CREATE_BASELINE} = false ]; then
#
# --- regression test comparison ----
#

  for i in ${LIST_FILES} ; do
    printf %s " Comparing " $i "....." >> ${REGRESSIONTEST_LOG}
    printf %s " Comparing " $i "....."

    if [ ! -f ${RUNDIR}/$i ] ; then

     echo ".......MISSING file" >> ${REGRESSIONTEST_LOG}
     echo ".......MISSING file"
     test_status='FAIL'

    elif [ ! -f ${RTPWD}/${CNTL_DIR}/$i ] ; then

     echo ".......MISSING baseline" >> ${REGRESSIONTEST_LOG}
     echo ".......MISSING baseline"
     test_status='FAIL'

    else

     d=`cmp ${RTPWD}/${CNTL_DIR}/$i ${RUNDIR}/$i | wc -l`

     if [[ $d -ne 0 ]] ; then
       echo ".......NOT OK" >> ${REGRESSIONTEST_LOG}
       echo ".......NOT OK"
       test_status='FAIL'

     else

       echo "....OK" >> ${REGRESSIONTEST_LOG}
       echo "....OK"
     fi

    fi

  done

if [ $test_status = 'FAIL' ]; then echo $TEST_NAME >> fail_test ; fi

#
else
#
# --- create baselines
#

 sleep 60
 echo;echo;echo "Moving set ${TEST_NR} files ...."

 for i in ${LIST_FILES} ; do
  printf %s " Moving " $i "....."
  ls -ltr ${RUNDIR}/${i}
  if [ -f ${RUNDIR}/$i ] ; then
    cp ${RUNDIR}/${i} /${STMP}/${USER}/REGRESSION_TEST/${CNTL_DIR}/${i}
  else
    echo "Missing " ${RUNDIR}/$i " output file"
    echo;echo " Set ${TEST_NR} failed "
    exit 2
  fi
 done

# ---
fi
# ---

echo "Test ${TEST_NR} ${test_status} " >> ${REGRESSIONTEST_LOG}
(echo;echo;echo)                       >> ${REGRESSIONTEST_LOG}
echo "Test ${TEST_NR} ${test_status} "
(echo;echo;echo)

sleep 4
echo;echo

####################################################################################################
# End test
####################################################################################################

exit 0
