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
# CNTL_PTH    - control path for current test
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

#########################################################################
# Setup directory names
#########################################################################
export RTPWD=/${DISKNM}/noscrub/wx20rv/REGRESSION_TEST
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

cd $PATHRT

####################################################################################################
export TEST_DESCR="Compare NMMB-global results with previous trunk version"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_CNTRL
export CNTL_PTH=${RTPWD}/NMMB_glob
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06    ; export WTPG=2      ; export FCSTL=48
export NEMSI=false ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=true
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

export timing1=`grep total_tim $PATHRT/err | tail -1 | awk '{ print $5 }'`
export timingc=`cat ${RTPWD}/NMMB_glob/timing.txt`
(echo " Original timing: " $timingc " , test_glob timing: " $timing1;echo;echo)>> RegressionTests.log
 echo " Original timing: " $timingc " , test_glob timing: " $timing1;echo;echo

####################################################################################################
export TEST_DESCR="Test NMMB-global NEMSIO as input file"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_NEMSIO
export CNTL_PTH=${RTPWD}/NMMB_glob
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06    ; export WTPG=2      ; export FCSTL=48
export NEMSI=true  ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-global restart run"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REST
export CNTL_PTH=${RTPWD}/NMMB_glob
export LIST_FILES="nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06    ; export WTPG=2      ; export FCSTL=48
export NEMSI=false ; export RSTRT=true  ; export gfsP=false ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-global restart run from NEMSIO file"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REST_NIO
export CNTL_PTH=${RTPWD}/NMMB_glob
export LIST_FILES="nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06    ; export WTPG=2      ; export FCSTL=48
export NEMSI=true  ; export RSTRT=true  ; export gfsP=false ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-global different decomposition"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_DECOMP
export CNTL_PTH=${RTPWD}/NMMB_glob
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024   \
nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                  \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024"
#---------------------
export TPN=16      ; export THRD=1      ; export GS=''      ; export GBRG=glob
export INPES=03    ; export WTPG=1      ; export FCSTL=24
export NEMSI=false ; export RSTRT=false ; export gfsP=false ; export RGS=true  ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-global threading "
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_THREAD
export CNTL_PTH=${RTPWD}/NMMB_glob
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=''      ; export GBRG=glob
export INPES=06    ; export WTPG=2      ; export FCSTL=48
export NEMSI=false ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-global with GFS physics package "
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_gfsP
export CNTL_PTH=${RTPWD}/NMMB_gfsP_glob
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024         \
nmm_b_restart.012 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006      \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_restart_nemsio.012"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=glob
export INPES=06    ; export WTPG=2      ; export FCSTL=24
export NEMSI=false ; export RSTRT=false ; export gfsP=true  ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Compare NMMB-regional results with previous trunk version"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_CTL
export CNTL_PTH=${RTPWD}/NMMB_reg
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06    ; export WTPG=2      ; export FCSTL=48
export NEMSI=false ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=true
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

export timing1=`grep total_tim $PATHRT/err | tail -1 | awk '{ print $5 }'`
export timingc=`cat ${RTPWD}/NMMB_reg/timing.txt`
(echo " Original timing: " $timingc " , test_reg timing: " $timing1;echo;echo)>> RegressionTests.log
 echo " Original timing: " $timingc " , test_reg timing: " $timing1;echo;echo

####################################################################################################
export TEST_DESCR="Test NMMB-regional NEMSIO as input file"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NEM_REG_NEMSIO
export CNTL_PTH=${RTPWD}/NMMB_reg
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.009 nmm_b_history.012    \
nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                   \
nmm_b_history_nemsio.009 nmm_b_history_nemsio.012"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06    ; export WTPG=2      ; export FCSTL=12
export NEMSI=true  ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-regional restart run"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_RST
export CNTL_PTH=${RTPWD}/NMMB_reg
export LIST_FILES="nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06    ; export WTPG=2      ; export FCSTL=48
export NEMSI=false ; export RSTRT=true  ; export gfsP=false ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-regional restart run with NEMSIO file "
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_RST_NIO
export CNTL_PTH=${RTPWD}/NMMB_reg
export LIST_FILES="nmm_b_history.027 nmm_b_history.030 nmm_b_history.036 nmm_b_history.048"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06    ; export WTPG=2      ; export FCSTL=48
export NEMSI=true  ; export RSTRT=true  ; export gfsP=false ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-regional different decomposition"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_DECOMP
export CNTL_PTH=${RTPWD}/NMMB_reg
export LIST_FILES="nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.009 nmm_b_history.012"
#---------------------
export TPN=16      ; export THRD=1      ; export GS=''      ; export GBRG=reg
export INPES=03    ; export WTPG=1      ; export FCSTL=12
export NEMSI=false ; export RSTRT=false ; export gfsP=false ; export RGS=true  ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-regional threading "
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_THREAD
export CNTL_PTH=${RTPWD}/NMMB_reg
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024 nmm_b_history.048 \
nmm_b_restart.024 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006                \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_history_nemsio.048 nmm_b_restart_nemsio.024"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=''      ; export GBRG=reg
export INPES=06    ; export WTPG=2      ; export FCSTL=48
export NEMSI=false ; export RSTRT=false ; export gfsP=false ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test NMMB-regional with GFS physics package "
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_gfsP
export CNTL_PTH=${RTPWD}/NMMB_gfsP_reg
export LIST_FILES=" \
nmm_b_history.000 nmm_b_history.003 nmm_b_history.006 nmm_b_history.012 nmm_b_history.024     \
nmm_b_restart.012 nmm_b_history_nemsio.000 nmm_b_history_nemsio.003 nmm_b_history_nemsio.006  \
nmm_b_history_nemsio.012 nmm_b_history_nemsio.024 nmm_b_restart_nemsio.012"
#---------------------
export TPN=32      ; export THRD=1      ; export GS=#       ; export GBRG=reg
export INPES=06    ; export WTPG=2      ; export FCSTL=24
export NEMSI=false ; export RSTRT=false ; export gfsP=true  ; export RGS=false ; export WGS=false
#---------------------
./rt_nmm.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

#########################################################################
# Clean and compile GFS core
#########################################################################
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

cd $PATHRT

####################################################################################################
export TEST_DESCR="Compare GFS results with previous trunk version"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32
export CNTL_PTH=${RTPWD}/GFS
export LIST_FILES=" \
       sigf03 sigf06 sigf12 sigf24 sigf48 \
       sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
       flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=32    ; export THRD=1       ; export NSOUT=0     ; export QUILT=.true.
export PE1=30      ; export WTPG=2       ; export NDAYS=2     ; export CP2=#
#---------------------
./rt_gfs.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test GFS with 2-copy option"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_2copy
export CNTL_PTH=${RTPWD}/GFS
export LIST_FILES=" \
       sigf03 sigf06 sigf12 sigf24 sigf48 \
       sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
       flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=32    ; export THRD=1       ; export NSOUT=0     ; export QUILT=.true.
export PE1=30      ; export WTPG=2       ; export NDAYS=2     ; export CP2=''
#---------------------
./rt_gfs.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test GFS different decomposition"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_60
export CNTL_PTH=${RTPWD}/GFS
export LIST_FILES=" \
       sigf03 sigf06 sigf12 sigf24 sigf48 \
       sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
       flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=60    ; export THRD=1       ; export NSOUT=0     ; export QUILT=.true.
export PE1=58      ; export WTPG=2       ; export NDAYS=2     ; export CP2=#
#---------------------
./rt_gfs.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test GFS threads"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16
export CNTL_PTH=${RTPWD}/GFS
export LIST_FILES=" \
       sigf03 sigf06 sigf12 sigf24 sigf48 \
       sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
       flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=16    ; export THRD=2       ; export NSOUT=0     ; export QUILT=.true.
export PE1=14      ; export WTPG=2       ; export NDAYS=2     ; export CP2=#
#---------------------
./rt_gfs.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test GFS single processor"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_01
export CNTL_PTH=${RTPWD}/GFS
export LIST_FILES=" \
         sigf03 sigf06 sigf12 sigf24 \
         sfcf03 sfcf06 sfcf12 sfcf24 \
         flxf03 flxf06 flxf12 flxf24"
#---------------------
export TASKS=1     ; export THRD=1       ; export NSOUT=0     ; export QUILT=.false.
export PE1=1       ; export WTPG=1       ; export NDAYS=1     ; export CP2=#
#---------------------
./rt_gfs.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test GFS, 1 proc, 1 threads,no quilting,nsout=1"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_01_NSOUT
export CNTL_PTH=${RTPWD}/GFS
export LIST_FILES=" \
         sigf03 sigf06 sigf09 sigf12 sigf24 \
         sfcf03 sfcf06 sigf09 sfcf12 sfcf24 \
         flxf03 flxf06 sigf09 flxf12 flxf24"
#---------------------
export TASKS=1     ; export THRD=1       ; export NSOUT=1     ; export QUILT=.false.
export PE1=1       ; export WTPG=1       ; export NDAYS=1     ; export CP2=#
#---------------------
./rt_gfs.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test GFS, 16 proc, 2 threads,no quilt, output every 2 time steps"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16_NOQUILT_NSOUT
export CNTL_PTH=${RTPWD}/GFS
export LIST_FILES=" \
       sigf03 sigf06 sigf12 sigf24 sigf48 \
       sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
       flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=16    ; export THRD=2       ; export NSOUT=2     ; export QUILT=.false.
export PE1=16      ; export WTPG=1       ; export NDAYS=2     ; export CP2=#
#---------------------
./rt_gfs.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="Test GFS, 60 proc, 1 thread, no quilt"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_60_NOQUILT
export CNTL_PTH=${RTPWD}/GFS
export LIST_FILES=" \
       sigf03 sigf06 sigf12 sigf24 sigf48 \
       sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
       flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=60    ; export THRD=1       ; export NSOUT=1     ; export QUILT=.false.
export PE1=60      ; export WTPG=1       ; export NDAYS=2     ; export CP2=#
#---------------------
./rt_gfs.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
export TEST_DESCR="GFS, 32 proc, 1 thread, no quilt, output every 4 timestep"
####################################################################################################
#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_NOQUILT
export CNTL_PTH=${RTPWD}/GFS
export LIST_FILES=" \
       sigf03 sigf06 sigf12 sigf24 sigf48 \
       sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
       flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export TASKS=32    ; export THRD=1       ; export NSOUT=4     ; export QUILT=.false.
export PE1=32      ; export WTPG=1       ; export NDAYS=2     ; export CP2=#
#---------------------
./rt_gfs.sh
 if [ $? = 2 ]; then exit ; fi
#---------------------

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
