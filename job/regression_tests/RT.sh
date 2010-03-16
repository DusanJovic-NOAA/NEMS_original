#!/bin/ksh
set -ua

#########################################################################
#         USER DEFINED PART!!!!!!
# Setup machine related variables (CLASS, ACCNR, & DISKNM)
#########################################################################

export MACHINE_ID=`hostname | cut -c1`

if [ ${MACHINE_ID} = c -o ${MACHINE_ID} = s ]; then
  export CLASS=dev
  export ACCNR=NAM-T2O
  export DISKNM=meso
elif [ ${MACHINE_ID} = v ]; then 
  export CLASS=mtb
  export ACCNR=MTB003-RES
  export DISKNM=mtb
fi

############################################################
# RTPWD - Path to previously stored regression test answers
############################################################

export RTPWD=/${DISKNM}/noscrub/wx20rv/REGRESSION_TEST

#########################################################################
# Check if running regression test or creating baselines.
# If one argument is provided and argument is = create_baselines
# then this script will prepare new answers for regression tests
#########################################################################

argn=$#
  export CREATE_BASELINE=false
  CB_arg=all
if [ $argn -eq 1 ]; then
  export CREATE_BASELINE=true
  CB_arg=$1
  if [ ${CB_arg} != nmm -a ${CB_arg} != gfs -a ${CB_arg} != all ]; then
    echo "Wrong CB_arg choice: " $CB_arg
    echo "  Options are: "
    echo "          RT.sh  nmm  (create baselines for NMM)"
    echo "          RT.sh  gfs  (create baselines for GFS)"
    echo "          RT.sh  all  (create baselines for both NMM & GFS)"
    exit
  fi
  #
  # prepare new regression test directory
  #
  rm -rf /stmp/${LOGIN}/REGRESSION_TEST
  cp -r /${DISKNM}/noscrub/wx20rv/REGRESSION_TEST_baselines \
	/stmp/${LOGIN}/REGRESSION_TEST
  if [ ${CB_arg} = nmm ]; then
    cp ${RTPWD}/GFS_DFI_REDUCEDGRID/*  /stmp/${LOGIN}/REGRESSION_TEST/GFS_DFI_REDUCEDGRID/.
    cp ${RTPWD}/GFS_NODFI/*            /stmp/${LOGIN}/REGRESSION_TEST/GFS_NODFI/.
    cp ${RTPWD}/GFS_OPAC/*             /stmp/${LOGIN}/REGRESSION_TEST/GFS_OPAC/.
    cp ${RTPWD}/GEFS_data_2008082500/* /stmp/${LOGIN}/REGRESSION_TEST/GEFS_data_2008082500/.
    cp ${RTPWD}/GEFS_m4/*              /stmp/${LOGIN}/REGRESSION_TEST/GEFS_m4/.
  elif [ ${CB_arg} = gfs ]; then
    cp ${RTPWD}/NMMB_gfsP_glob/*      /stmp/${LOGIN}/REGRESSION_TEST/NMMB_gfsP_glob/.
    cp ${RTPWD}/NMMB_gfsP_reg/*       /stmp/${LOGIN}/REGRESSION_TEST/NMMB_gfsP_reg/.
    cp ${RTPWD}/NMMB_glob/*           /stmp/${LOGIN}/REGRESSION_TEST/NMMB_glob/.
    cp ${RTPWD}/NMMB_nests/*          /stmp/${LOGIN}/REGRESSION_TEST/NMMB_nests/.
    cp ${RTPWD}/NMMB_reg/*            /stmp/${LOGIN}/REGRESSION_TEST/NMMB_reg/.
    cp ${RTPWD}/NMMB_reg_pcpadj/*     /stmp/${LOGIN}/REGRESSION_TEST/NMMB_reg_pcpadj/.
  fi
fi

################################################
# List of variables in use:
################################################
#
# TEST_NR     - test number
# TEST_DESCR  - test description
# RTPWD       - path with previous stored data
# PATHRT      - regression test script path
# PATHTR      - NEMS path
# RUNDIR_ROOT - temporary run directory
# RUNDIR      - current test run directory
# CNTL_DIR    - control directory name for current test
# LIST_FILES  - list of files for comparison
# TPN         - number of tasks per node
# TASKS       - total number of tasks
# PE1         - number of computing tasks
# WTPG        - number of write tasks
# THRD        - number of threads
# GS          - read global sums (' '), or not (#)
# INPES       - number od PE's on x direction
# FCSTL       - forecast length in hours
# NDAYS       - forecast length in days
# NEMSI       - NEMSIO as input file
# RSTRT       - restarted run
# gfsP        - GFS physics suite
# RGS         - read global sums
# WGS         - write global sums
# GBRG        - NMMB global/regional option
# QUILT       - quilting ON/OFF
# NSOUT       - number of timesteps for output
# CLASS       - job class/group (LoadLeveler)
# ACCNR       - account number (LoadLeveler)
# DISKNM      - disk name ( /meso or /mtb)
# CP2         - 2-copy option
#
################################################

date > Compile.log
date > RegressionTests.log
echo "Start Regression test" >> RegressionTests.log
(echo;echo;echo)             >> RegressionTests.log

###################################
# PATHRT - Path to regression test
###################################

export PATHRT=`pwd`
cd ../../
export PATHTR=`pwd`

export RUNDIR_ROOT=/ptmp/${LOGIN}/RT_$$
typeset -Z3 TEST_NR
export TEST_NR=0

clear;echo;echo

#########################################################################
# Clean and compile NMMB core
#########################################################################
if [ ${CB_arg} != gfs ]; then
  echo "Preparing NMMB core for regression tests"
  printf %s "Compiling NMMB core (this will take some time)......."
  cd ${PATHTR}/ush

  ./clean_stub.sh                        >> ${PATHRT}/Compile.log 2>&1
  ./clean.sh                             >> ${PATHRT}/Compile.log 2>&1
  ./compile_configure.sh nmm gfs_physics >> ${PATHRT}/Compile.log 2>&1
  ./compile.sh                           >> ${PATHRT}/Compile.log 2>&1

  if [ -f ../exe/NMM_NEMS.x ] ; then
    echo "   NMMB core Compiled";echo;echo
  else
    echo "   NMMB core is NOT compiled" >> ${PATHRT}/RegressionTests.log
    echo "   NMMB core is NOT compiled"
    exit
  fi
fi

cd $PATHRT

####################################################################################################
#
# TEST   - Global NMM-B with pure binary input
#        - 6x5 compute  tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs ]; then

export TEST_DESCR="Compare NMMB-global results with previous trunk version"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_CNTRL
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmm_b_history.000         nmm_b_history.003         nmm_b_history.006         nmm_b_history.012        \
nmm_b_history.024         nmm_b_history.027         nmm_b_history.030         nmm_b_history.036        \
nmm_b_history.048         nmm_b_history_nemsio.000  nmm_b_history_nemsio.003  nmm_b_history_nemsio.006 \
nmm_b_history_nemsio.012  nmm_b_history_nemsio.024  nmm_b_history_nemsio.048  nmm_b_restart.024        \
nmm_b_restart_nemsio.024"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=48
export NEMSI=false  ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=true
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

export timing1=`grep total_integration_tim $PATHRT/err | tail -1 | awk '{ print $5 }'`
export timingc=`cat ${RTPWD}/NMMB_glob/timing.txt`
(echo " Original timing: " $timingc " , test_glob timing: " $timing1;echo;echo)>> RegressionTests.log
 echo " Original timing: " $timingc " , test_glob timing: " $timing1;echo;echo

fi

####################################################################################################
#
# TEST   - Global NMM-B with NEMSIO input
#        - 6x5 compute tasks / 1 thread / opnl physics / free fcst / nemsio input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-global NEMSIO as input file"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_NEMSIO
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=48
export NEMSI=true   ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Global NMM-B restart from pure binary input
#        - 6x5 compute tasks / 1 thread / opnl physics / restart / pure binary input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-global restart run"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REST
export CNTL_DIR=NMMB_glob
export LIST_FILES="nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=48
export NEMSI=false  ; export RSTRT=true  ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Global NMM-B restart from NEMSIO input
#        - 6x5 compute tasks / 1 thread / opnl physics / restart / nemsio input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-global restart run from NEMSIO file"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REST_NIO
export CNTL_DIR=NMMB_glob
export LIST_FILES="nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=48
export NEMSI=true   ; export RSTRT=true  ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Global NMM-B with different domain decomposition
#        - 3x5 compute tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-global different decomposition"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_DECOMP
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024   \
nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                  \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024"
#---------------------
export TPN=16       ; export THRD=1      ; export GS=''      ; export GBRG=glob
export INPES=03     ; export JNPES=05    ; export WTPG=1     ; export FCSTL=24
export NEMSI=false  ; export RSTRT=false ; export gfsP=false ; export RGS=true  ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Global NMM-B with multiple threads
#        - 6x5 compute tasks / 2 threads / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-global threading "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_THREAD
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024"
#---------------------
export TPN=32       ; export THRD=2      ; export GS=''      ; export GBRG=glob
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=48
export NEMSI=false  ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Global NMM-B with GFS physics
#        - 6x5 compute tasks / 1 thread / GFS physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs ]; then
 
export TEST_DESCR="Test NMMB-global with GFS physics package "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_gfsP
export CNTL_DIR=NMMB_gfsP_glob
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024         \
nmm_b_restart.012 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006      \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_restart_nemsio.012"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=24
export NEMSI=false  ; export RSTRT=false ; export gfsP=true  ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Regional NMM-B with pure binary input
#        - 6x5 compute tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs ]; then

export TEST_DESCR="Compare NMMB-regional results with previous trunk version"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_CTL
export CNTL_DIR=NMMB_reg
export LIST_FILES=" \
nmm_b_history.000         nmm_b_history.003         nmm_b_history.006         nmm_b_history.012        \
nmm_b_history.024         nmm_b_history.027         nmm_b_history.030         nmm_b_history.036        \
nmm_b_history.048         nmm_b_history_nemsio.000  nmm_b_history_nemsio.003  nmm_b_history_nemsio.006 \
nmm_b_history_nemsio.012  nmm_b_history_nemsio.024  nmm_b_history_nemsio.048  nmm_b_restart.024        \
nmm_b_restart_nemsio.024"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=48
export NEMSI=false  ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=true
export PCPFLG=false ; export WPREC=true  ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

export timing1=`grep total_integration_tim $PATHRT/err | tail -1 | awk '{ print $5 }'`
export timingc=`cat ${RTPWD}/NMMB_reg/timing.txt`
(echo " Original timing: " $timingc " , test_reg timing: " $timing1;echo;echo)>> RegressionTests.log
 echo " Original timing: " $timingc " , test_reg timing: " $timing1;echo;echo

fi

####################################################################################################
# 
# TEST   - Regional NMM-B with NEMSIO input
#        - 6x5 compute tasks / 1 thread / opnl physics / free fcst / nemsio input
# 
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-regional NEMSIO as input file"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_NEMSIO
export CNTL_DIR=NMMB_reg
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012      \
nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006   \
nmm_b_history_nemsio.012"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=12
export NEMSI=true   ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Regional NMM-B restart with pure binary input
#        - 6x5 compute tasks / 1 thread / opnl physics / restart / pure binary input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-regional restart run"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_RST
export CNTL_DIR=NMMB_reg
export LIST_FILES="nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=48
export NEMSI=false  ; export RSTRT=true  ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Regional NMM-B restart with NEMSIO input
#        - 6x5 compute tasks / 1 thread / opnl physics / restart / nemsio input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-regional restart run with NEMSIO file "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_RST_NIO
export CNTL_DIR=NMMB_reg
export LIST_FILES="nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=48
export NEMSI=true   ; export RSTRT=true  ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Regional NMM-B with different domain decomposition
#        - 3x5 compute tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-regional different decomposition"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_DECOMP
export CNTL_DIR=NMMB_reg
export LIST_FILES="nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012"
#---------------------
export TPN=16       ; export THRD=1      ; export GS=''      ; export GBRG=reg
export INPES=03     ; export JNPES=05    ; export WTPG=1     ; export FCSTL=12
export NEMSI=false  ; export RSTRT=false ; export gfsP=false ; export RGS=true  ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Regional NMM-B with multiple threads
#        - 6x5 compute tasks / 2 threads / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-regional threading "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_THREAD
export CNTL_DIR=NMMB_reg
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024"
#---------------------
export TPN=32       ; export THRD=2      ; export GS=''      ; export GBRG=reg
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=48
export NEMSI=false  ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Regional NMM-B with GFS physics
#        - 6x5 compute tasks / 1 thread / GFS physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs ]; then

export TEST_DESCR="Test NMMB-regional with GFS physics package "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_gfsP
export CNTL_DIR=NMMB_gfsP_reg
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024     \
nmm_b_restart.012 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006  \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_restart_nemsio.012"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=24
export NEMSI=false  ; export RSTRT=false ; export gfsP=true  ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - NMM-B nesting: Regional parent with two children and one grandchild
#        - Compute tasks - Upper parent 2x3 | Child #1 4x8 | Child #2 2x4 | Grandchild 7x10
#        - 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs ]; then

export TEST_DESCR="Test NMMB-regional with nesting"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_nests
export CNTL_DIR=NMMB_nests
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024                \
nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                               \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024                                                        \
nmm_b_restart.012 nmm_b_restart.024 nmm_b_restart_nemsio.012 nmm_b_restart_nemsio.024                    \
nmm_b_history.02.000 nmm_b_history.02.003 nmm_b_history.02.006 nmm_b_history.02.012 nmm_b_history.02.024 \
nmm_b_history.02_nemsio.000 nmm_b_history.02_nemsio.003 nmm_b_history.02_nemsio.006                      \
nmm_b_history.02_nemsio.012 nmm_b_history.02_nemsio.024                                                  \
nmm_b_restart.02.012 nmm_b_restart.02.024 nmm_b_restart.02_nemsio.012 nmm_b_restart.02_nemsio.024        \
nmm_b_history.03.000 nmm_b_history.03.003 nmm_b_history.03.006 nmm_b_history.03.012 nmm_b_history.03.024 \
nmm_b_history.03_nemsio.000 nmm_b_history.03_nemsio.003 nmm_b_history.03_nemsio.006                      \
nmm_b_history.03_nemsio.012 nmm_b_history.03_nemsio.024                                                  \
nmm_b_restart.03.012 nmm_b_restart.03.024 nmm_b_restart.03_nemsio.012 nmm_b_restart.03_nemsio.024        \
nmm_b_history.04.000 nmm_b_history.04.003 nmm_b_history.04.006 nmm_b_history.04.012 nmm_b_history.04.024 \
nmm_b_history.04_nemsio.000 nmm_b_history.04_nemsio.003 nmm_b_history.04_nemsio.006                      \
nmm_b_history.04_nemsio.012 nmm_b_history.04_nemsio.024                                                  \
nmm_b_restart.04.012 nmm_b_restart.04.024 nmm_b_restart.04_nemsio.012 nmm_b_restart.04_nemsio.024"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=nests
export INPES=02     ; export JNPES=03    ; export WTPG=1     ; export FCSTL=24
export NEMSI=false  ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=false ; export WPREC=false ; export CPPCP=#    ; export NCHILD=02
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Regional NMM-B with precipitation adjustment on
#        - 6x5 compute tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs ]; then

export TEST_DESCR="Test NMMB-regional with precipitation adjustment on"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_PCPADJ
export CNTL_DIR=NMMB_reg_pcpadj
export LIST_FILES=" \
nmm_b_history.000  nmm_b_history.012         nmm_b_history_nemsio.006  nmm_b_restart_nemsio.012 \
nmm_b_history.003  nmm_b_history_nemsio.000  nmm_b_history_nemsio.012                           \
nmm_b_history.006  nmm_b_history_nemsio.003  nmm_b_restart.012"
#---------------------
export TPN=32       ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06     ; export JNPES=05    ; export WTPG=2     ; export FCSTL=12
export NEMSI=false  ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
export PCPFLG=true  ; export WPREC=false ; export CPPCP=''   ; export NCHILD=0
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
# Clean and compile GFS core
####################################################################################################

if [ ${CB_arg} != nmm ]; then

  echo "Preparing GFS core for regression tests" >> RegressionTests.log
  echo "Preparing GFS core for regression tests"
  printf %s "Compiling GFS core (this will take some time)......."
  cd ${PATHTR}/ush

  ./clean_stub.sh            >> ${PATHRT}/Compile.log 2>&1
  ./clean.sh                 >> ${PATHRT}/Compile.log 2>&1
  ./compile_configure.sh gfs >> ${PATHRT}/Compile.log 2>&1
  ./compile.sh               >> ${PATHRT}/Compile.log 2>&1

  if [ -f ../exe/GFS_NEMS.x ] ; then
    echo "   GFS core Compiled";echo;echo
  else
    echo "   GFS core is NOT compiled" >> ${PATHRT}/RegressionTests.log
    echo "   GFS core is NOT compiled"
    exit
  fi
fi

cd $PATHRT

####################################################################################################
# 
# TEST   - GFS 
#        - 30 compute tasks / 1 thread 
#
####################################################################################################

if [ ${CB_arg} != nmm ]; then

export TEST_DESCR="Compare GFS results with previous trunk version"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24 "
#---------------------
export TASKS=32    ; export THRD=1       ; export NSOUT=0     ; export QUILT=.true.
export PE1=30      ; export WTPG=2       ; export NDAYS=1     ; export CP2=#
export WRTGP=1     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS as two copies
#        - 30 compute tasks / 1 thread  2 copy restart
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test GFS with 2-copy option"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=32    ; export THRD=1       ; export NSOUT=0     ; export QUILT=.true.
export PE1=30      ; export WTPG=2       ; export NDAYS=2     ; export CP2=''
export WRTGP=1     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS with different decomposition
#        - 58 compute tasks / 1 thread 
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test GFS different decomposition"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_60_16
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export TASKS=60    ; export THRD=1       ; export NSOUT=0     ; export QUILT=.true.
export PE1=58      ; export WTPG=2       ; export NDAYS=1     ; export CP2=#
export WRTGP=1     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS with multiple threads
#        - 12 compute tasks / 2 threads / 2WrtGrp & 2WrtPePerGrp
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test GFS threads"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_60_16
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
       sigf03 sigf06 sigf12 sigf24 sigf48 \
       sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
       flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=16    ; export THRD=2       ; export NSOUT=0     ; export QUILT=.true.
export PE1=12      ; export WTPG=2       ; export NDAYS=2     ; export CP2=#
export WRTGP=2     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS on a single processor
#        - 1 task / 1 thread
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test GFS single processor"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_01_NSOUT
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export TASKS=1     ; export THRD=1       ; export NSOUT=0     ; export QUILT=.false.
export PE1=1       ; export WTPG=1       ; export NDAYS=1     ; export CP2=#
export WRTGP=1     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS on a single processor with no quilting
#        - 1 task / 1 thread
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test GFS, 1 proc, 1 thread, no quilting,nsout=1"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_01_NSOUT
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=1     ; export THRD=1       ; export NSOUT=1     ; export QUILT=.false.
export PE1=1       ; export WTPG=1       ; export NDAYS=2     ; export CP2=#
export WRTGP=1     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS with multiple threads, no quilting, and frequent output
#        - 16 tasks / 2 threads
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test GFS, 16 proc, 2 threads,no quilt, output every 2 time steps"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16_60_NOQUILT_NSOUT
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export TASKS=16    ; export THRD=2       ; export NSOUT=2     ; export QUILT=.false.
export PE1=16      ; export WTPG=1       ; export NDAYS=1     ; export CP2=#
export WRTGP=1     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS with multiple tasks and no quilting
#        - 60 tasks / 1 thread
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test GFS, 60 proc, 1 thread, no quilt"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16_60_NOQUILT_NSOUT
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=60    ; export THRD=1       ; export NSOUT=1     ; export QUILT=.false.
export PE1=60      ; export WTPG=1       ; export NDAYS=2     ; export CP2=#
export WRTGP=1     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS with multiple tasks, no quilting, and frequent output
#        - 32 tasks / 1 thread
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="GFS, 32 proc, 1 thread, no quilt, output every 4 timestep"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_NOQUILT
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=32    ; export THRD=1       ; export NSOUT=4     ; export QUILT=.false.
export PE1=32      ; export WTPG=1       ; export NDAYS=2     ; export CP2=#
export WRTGP=1     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS digital filter
#        - 30 compute tasks / 1 thread
#
####################################################################################################

if [ ${CB_arg} != nmm ]; then

export TEST_DESCR="GFS,32 proc, 1 thread, quilt, digital filter on reduced grid"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_16_dfi
export CNTL_DIR=GFS_DFI_REDUCEDGRID
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export TASKS=32    ; export THRD=1       ; export NSOUT=0     ; export QUILT=.true.
export PE1=30      ; export WTPG=2       ; export NDAYS=1     ; export CP2=#
export WRTGP=1     ; export FDFI=3      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS digital filter
#        - 12 compute tasks / 2 thread ,2 WrtGrp x 2 WrtPePerGrp
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="GFS,16 proc, 2 thread, quilt,2x2 wrt pe, digital filter on reduced grid"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_16_dfi
export CNTL_DIR=GFS_DFI_REDUCEDGRID
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
        sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
        flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=16    ; export THRD=2       ; export NSOUT=0     ; export QUILT=.true.
export PE1=12      ; export WTPG=2       ; export NDAYS=2     ; export CP2=''
export WRTGP=2     ; export FDFI=3      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS digital filter
#        - 1pe nsout=1
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="GFS,1 proc, no quilt, digital filter on reduced grid"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_1_dfi
export CNTL_DIR=GFS_DFI_REDUCEDGRID
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export TASKS=1     ; export THRD=1       ; export NSOUT=0     ; export QUILT=.false.
export PE1=1       ; export WTPG=1       ; export NDAYS=1     ; export CP2=#
export WRTGP=1     ; export FDFI=3      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=0      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
# 
# TEST   - GFS 
#        - OPAC aerosols
#
####################################################################################################

if [ ${CB_arg} != nmm ]; then

export TEST_DESCR="GFS, use the OPAC climo scheme for SW and LW"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_OPAC
export CNTL_DIR=GFS_OPAC
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=32    ; export THRD=1       ; export NSOUT=0     ; export QUILT=.true.
export PE1=30      ; export WTPG=2       ; export NDAYS=2     ; export CP2=#
export WRTGP=1     ; export FDFI=0      ; export ADIAB=.false.; export REDUCEDGRID=.true.
export NUMFILE=3   ; export IAER=11      ; export FHRES=24
export wave=62     ; export lm=64       ; export lsoil=4      ; export MEMBER_NAMES=c00
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Concurrency GEFS
#        - 4 members, every 6 hours, couple and add stochastic perturbations, T190L28.
#
####################################################################################################

if [ ${CB_arg} != nmm ]; then

export TEST_DESCR="Concurrency GEFS, stochastic perturbations, 4 members, T190L28."

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GEFS_Concurrency_Run
export CNTL_DIR=GEFS_m4
export LIST_FILES=" \
        SIG.F06_01 SIG.F06_02 SIG.F06_03 SIG.F06_04 \
        SIG.F12_01 SIG.F12_02 SIG.F12_03 SIG.F12_04 \
        SIG.F18_01 SIG.F18_02 SIG.F18_03 SIG.F18_04 \
        SIG.F24_01 SIG.F24_02 SIG.F24_03 SIG.F24_04 \
        SFC.F06_01 SFC.F06_02 SFC.F06_03 SFC.F06_04 \
        SFC.F12_01 SFC.F12_02 SFC.F12_03 SFC.F12_04 \
        SFC.F18_01 SFC.F18_02 SFC.F18_03 SFC.F18_04 \
        SFC.F24_01 SFC.F24_02 SFC.F24_03 SFC.F24_04 \
        FLX.F06_01 FLX.F06_02 FLX.F06_03 FLX.F06_04 \
        FLX.F12_01 FLX.F12_02 FLX.F12_03 FLX.F12_04 \
        FLX.F18_01 FLX.F18_02 FLX.F18_03 FLX.F18_04 \
        FLX.F24_01 FLX.F24_02 FLX.F24_03 FLX.F24_04"
#---------------------
export GEFS_ENSEMBLE=1
export TASKS=64    ; export THRD=1
echo 'PATHTR=' $PATHTR
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
# Finalize
####################################################################################################

rm -f err out configure_file nmm_ll gfs_fcst_run  gfs_ll
cd ${PATHTR}/ush
./clean_stub.sh      > /dev/null 2>&1
./clean.sh           > /dev/null 2>&1
rm -rf ${RUNDIR_ROOT}

date >> ${PATHRT}/RegressionTests.log

exit
