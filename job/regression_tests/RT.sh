#!/bin/ksh


#########################################################################
#         USER DEFINED PART!!!!!!
# Setup machine related variables (CLASS, ACCNR, & DISKNM)
#########################################################################

. ./detect_machine.sh

if [ ${MACHINE_ID} = c -o ${MACHINE_ID} = s ]; then
  export CLASS=dev
  export GROUP=dev
  export ACCNR=GFS-T2O
  export DISKNM=/meso
  export STMP=/stmp
  export PTMP=/ptmp
  export SCHEDULER=loadleveler
  STMP=/stmp
elif [ ${MACHINE_ID} = v ]; then 
  export CLASS=mtb
  export ACCNR=MTB003-RES
  export DISKNM=/mtb
  export STMP=/stmp
  export PTMP=/ptmp
  export SCHEDULER=loadleveler
  STMP=/stmp
elif [ ${MACHINE_ID} = g ]; then 
  export DISKNM=/lustre/ltfs/scratch/Ratko.Vasic
  export STMP=/lustre/fs/scratch
  export PTMP=/lustre/fs/scratch
  export SCHEDULER=moab
elif [ ${MACHINE_ID} = h ]; then 
  export DISKNM=/tds_scratch2/users/NCEPDEV/meso
  export STMP=/tds_scratch2/users/NCEPDEV/stmp
  export PTMP=/tds_scratch2/users/NCEPDEV/ptmp
  export SCHEDULER=pbs
elif [ ${MACHINE_ID} = z ]; then 
  export ACCNR=rm
  export DISKNM=/scratch2/portfolios/NCEPDEV/meso
  export STMP=/scratch2/portfolios/NCEPDEV/stmp
  export PTMP=/scratch2/portfolios/NCEPDEV/ptmp
  export SCHEDULER=pbs
else
  echo "Unknown machine ID, please edit detect_machine.sh file"
  exit
fi

############################################################
# RTPWD - Path to previously stored regression test answers
############################################################

   export RTPWD=${DISKNM}/noscrub/wx20rv/REGRESSION_TEST
#  export RTPWD=${STMP}/${USER}/REGRESSION_TEST

#########################################################################
# Check if running regression test or creating baselines.
# If one argument is provided and argument is = nmm, gfs, gen, or all,
# then this script will prepare new answers for regression tests
#########################################################################

argn=$#
  export CREATE_BASELINE=false
  CB_arg=all
if [ $argn -eq 1 ]; then
  export CREATE_BASELINE=true
  CB_arg=$1
  if [ ${CB_arg} != nmm -a ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a  ${CB_arg} != post -a ${CB_arg} != all ]; then
    echo "Wrong CB_arg choice: " $CB_arg
    echo "  Options are: "
    echo "          RT.sh  nmm  (create baselines for NMM)"
    echo "          RT.sh  gfs  (create baselines for GFS)"
    echo "          RT.sh  gen  (create baselines for GEN)"
    echo "          RT.sh  fim  (create baselines for FIM)"
    echo "          RT.sh  all  (create baselines for NMM & GFS & FIM)"
    echo "          RT.sh  post  (create baselines for post)"
    exit
  fi
  #
  # prepare new regression test directory
  #
  rm -rf ${STMP}/${USER}/REGRESSION_TEST
  cp -r ${DISKNM}/noscrub/wx20rv/REGRESSION_TEST_baselines \
	${STMP}/${USER}/REGRESSION_TEST
  CP_nmm=false
  CP_gfs=false
  CP_fim=false
  CP_post=false
  if [ ${CB_arg} = nmm ]; then
    CP_gfs=true
    CP_fim=true
    CP_post=true
  elif [ ${CB_arg} = gfs ]; then
    CP_nmm=true
    CP_fim=true
    CP_post=true
  elif [ ${CB_arg} = fim ]; then
    CP_nmm=true
    CP_gfs=true
    CP_post=true
  elif [ ${CB_arg} = post ]; then
    CP_nmm=true
    CP_gfs=true
    CP_fim=true
  fi
  if [ ${CP_gfs} = true ]; then
    cp ${RTPWD}/GEFS_data_2008082500/*      ${STMP}/${USER}/REGRESSION_TEST/GEFS_data_2008082500/.
    cp ${RTPWD}/GEFS_m4/*                   ${STMP}/${USER}/REGRESSION_TEST/GEFS_m4/.
    cp ${RTPWD}/GFS_DFI_REDUCEDGRID/*       ${STMP}/${USER}/REGRESSION_TEST/GFS_DFI_REDUCEDGRID/.
    cp ${RTPWD}/GFS_DFI_REDUCEDGRID_HYB/*   ${STMP}/${USER}/REGRESSION_TEST/GFS_DFI_REDUCEDGRID_HYB/.
    cp ${RTPWD}/GFS_DFI_REDUCEDGRID_NDSL/*  ${STMP}/${USER}/REGRESSION_TEST/GFS_DFI_REDUCEDGRID_NDSL/.
    cp ${RTPWD}/GFS_DFI_hyb_2loop/*         ${STMP}/${USER}/REGRESSION_TEST/GFS_DFI_hyb_2loop/.
    cp ${RTPWD}/GFS_DFI_hyb_2loop_nst/*     ${STMP}/${USER}/REGRESSION_TEST/GFS_DFI_hyb_2loop_nst/.
    cp ${RTPWD}/GFS_NODFI/*                 ${STMP}/${USER}/REGRESSION_TEST/GFS_NODFI/.
    cp ${RTPWD}/GFS_OPAC/*                  ${STMP}/${USER}/REGRESSION_TEST/GFS_OPAC/.
    cp ${RTPWD}/GFS_ADIAB/*                  ${STMP}/${USER}/REGRESSION_TEST/GFS_ADIAB/.
  elif [ ${CP_nmm} = true ]; then
    cp ${RTPWD}/NMMB_gfsP_glob/*            ${STMP}/${USER}/REGRESSION_TEST/NMMB_gfsP_glob/.
    cp ${RTPWD}/NMMB_gfsP_reg/*             ${STMP}/${USER}/REGRESSION_TEST/NMMB_gfsP_reg/.
    cp ${RTPWD}/NMMB_glob/*                 ${STMP}/${USER}/REGRESSION_TEST/NMMB_glob/.
    cp ${RTPWD}/NMMB_mvg_nests/*            ${STMP}/${USER}/REGRESSION_TEST/NMMB_mvg_nests/.
    cp ${RTPWD}/NMMB_nests/*                ${STMP}/${USER}/REGRESSION_TEST/NMMB_nests/.
    cp ${RTPWD}/NMMB_reg/*                  ${STMP}/${USER}/REGRESSION_TEST/NMMB_reg/.
    cp ${RTPWD}/NMMB_reg_filt/*             ${STMP}/${USER}/REGRESSION_TEST/NMMB_reg_filt/.
    cp ${RTPWD}/NMMB_reg_pcpadj/*           ${STMP}/${USER}/REGRESSION_TEST/NMMB_reg_pcpadj/.
    cp ${RTPWD}/NMMB_reg_sel_phy/*          ${STMP}/${USER}/REGRESSION_TEST/NMMB_reg_sel_phy/.
    cp ${RTPWD}/NMMB_reg_timesr/*           ${STMP}/${USER}/REGRESSION_TEST/NMMB_reg_timesr/.
  elif [ ${CP_fim} = true ]; then
    mkdir -p ${STMP}/${USER}/REGRESSION_TEST/FIMdata
    cp ${RTPWD}/FIMdata/*                   ${STMP}/${USER}/REGRESSION_TEST/FIMdata/.
    #TODO:  generalize with $GLVL (see below)
    mkdir -p ${STMP}/${USER}/REGRESSION_TEST/FIM_G4L38_24hr
    cp ${RTPWD}/FIM_G4L38_24hr/*        ${STMP}/${USER}/REGRESSION_TEST/FIM_G4L38_24hr/.
  elif [ ${CP_post} = true ]; then
    cp ${RTPWD}/GFS_DFI_POST/*              ${STMP}/${USER}/REGRESSION_TEST/GFS_DFI_POST/.
    cp ${RTPWD}/NMMB_reg_post/*             ${STMP}/${USER}/REGRESSION_TEST/NMMB_reg_post/.
  fi
fi

#########################################################################
# If two arguments are provided ("full test") this will run complete set
# of regression tests, otherwise it will run only mandatory (standard) test.
#########################################################################

  RT_FULL=false
  if [ ${CREATE_BASELINE} = true ] ; then RT_FULL=true ; fi
if [ $argn -eq 2 ]; then
  RT_arg1=$1
  RT_arg2=$2
  if [ ${RT_arg1} != full -o ${RT_arg2} != test ]; then
    echo "Wrong RT_arg choice: " $RT_arg1 $RT_arg2
    echo "  Option is: "
    echo "          RT.sh  full test  (run complete regression test)"
    exit
  fi
  RT_FULL=true
fi

################################################
# List of variables in use:
################################################
#
# RT_FULL     - true: full test; false: standard test
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
# NODE        - number of nodes
# TASKS       - total number of tasks
# PE1         - number of computing tasks
# WRTGP       - number of write groups
# WTPG        - number of write tasks per group
# THRD        - number of threads
# WLCLK       - wall clock limit
# TS          - write time series
# INPES       - number od PE's on x direction
# JNPES       - number od PE's on y direction
# FCSTL       - forecast length in hours
# NDAYS       - forecast length in days
# NEMSI       - NEMSIO as input file
# RSTRT       - restarted run
# gfsP        - GFS physics suite
# PCPFLG      - precipitation adjustment (true/false)
# WPREC       - write output for precipitation adjustment (true/false)
# CPPCP       - read precipitation adjustment files
# NCHILD      - number of direct nested domains
# CONVC       - convective scheme (bmj, sas)
# MICRO       - microphysics scheme (fer, gfs)
# TURBL       - PBL scheme (myj,gfs)
# GBRG        - NMMB global/regional/nest/filter option
# QUILT       - quilting ON/OFF
# NSOUT       - number of timesteps for output
# CLASS       - job class (LoadLeveler)
# GROUP       - job group (LoadLeveler)
# ACCNR       - account number (LoadLeveler)
# MACHINE_ID  - =c (cirrus), =s (stratus), =v (vapor)
# DISKNM      - disk name ( /meso or /mtb)
# CREATE_BASELINE - true/false
# CB_arg      - baseline arguments:
#                 =nmm - create nmm baselines only
#                 =gfs - create gfs baselines only
#                 =fim - create fim baselines only
#                 =all - create all baselines
# CP2         - 2-copy option
# FDFI	      - half of digital filter window in hours (e.g. 3)
# ADIAB       - logical, .true. if the run is adiabatic (i.e. no physics)
# REDUCEDGRID - turn off/on to run on reduced grid
# IAER        - 3-digit aerosol flag (for volc, lw, sw)
#		  =0: turn all aeros effects off (sw, lw, volc)
#		  =1: use clim tropospheric aerosol for sw only
#		  =10: use clim tropospheric aerosol for lw only
#		  =11: use clim tropospheric aerosol for both sw/lw
#		  =100: volc aerosol only for both sw/lw
#		  =101: volc and clim trops aerosol for sw only
#		  =110: volc and clim trops aerosol for lw only
#		  =111: volc and clim trops aerosol for both sw/lw
# FHRES       - restart interval in hours
# wavecoef    - spectral truncation of the model wave (i.e. JCAP  like T878)
# wavegrid    - spectral truncation of the model grid (i.e. JCAPG like T574)
# lm          - number of model levels
# lsoil       - number of soil levels in the land model
# IDVC        - integer describing the vertical coord - 2 means hybrid, 3 means
#		generalized hybrid
# THERMODYN_ID- integer identifying the thermodynamic variable - virtual T
#		or enthalpy
# SFCPRESS_ID - integer identifying whether lnps or ps is the variable used
# SPECTRALLOOP- option for one loop or two loop
# NST_FCST    - integer determining whether NST model is run in forecast mode 
#		or not (values 0,1 or 2)
# NDSLFV      - option be true or false for non-iteration dimensional-split
#               semi-Lagrangian advection
# MEMBER_NAMES  - ensemble member names (c00 if only one member)
# GEFS_ENSEMBLE - GEFS ensemble (=1 true, =0 false)
# GEN_ENSEMBLE  - GEN ensemble (=1 true, =0 false)
# STMP/PTMP     - temporary directory names 
# SCHEDULER     - machine scheduler name
# AFFN          - IBM specific - task affinity (core or cpu)
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

export RUNDIR_ROOT=${PTMP}/${USER}/RT_$$
mkdir -p ${RUNDIR_ROOT}
typeset -Z3 TEST_NR
export TEST_NR=0

clear;echo;echo

####################################################################################################
#
# Export variables to the default values
#  - first common variables, then model specific ones
#  - different machines, different defaults:
#
####################################################################################################

if [ ${MACHINE_ID} = c -o ${MACHINE_ID} = s -o ${MACHINE_ID} = v ]; then
  WTPGm=2
  TPNm=32
  TPNn=64
  INPm=06
  JNPm=05
elif [ ${MACHINE_ID} = g -o ${MACHINE_ID} = h -o ${MACHINE_ID} = z ]; then 
  WTPGm=3
  TPNm=48
  TPNn=128
  INPm=05
  JNPm=09
fi

export_common ()
{
export THRD=1
export WTPG=$WTPGm
export WLCLK=05
export GEFS_ENSEMBLE=0
export GEN_ENSEMBLE=0
export WRITE_DOPOST=.false.
export POST_GRIBVERSION='"grib1"'
}

export_nmm ()
{
export_common
export GBRG=reg     ; export TPN=$TPNm   ; export INPES=$INPm ; export JNPES=$JNPm
export AFFN=core    ; export NODE=1
export NEMSI=false  ; export RSTRT=false ; export gfsP=false  ; export FCSTL=48
export PCPFLG=false ; export WPREC=false ; export CPPCP=#     ; export TS=#
export NCHILD=0     ; export CONVC=bmj   ; export MICRO=fer   ; export TURBL=myj
}

export_gfs ()
{
export_common
export TASKS=32    ; export PE1=30          ; export NSOUT=0       ; export QUILT=.true.
export NDAYS=2     ; export CP2=.false.     ; export IAER=0        ; export FHRES=24
export WRTGP=1     ; export FDFI=0          ; export ADIAB=.false. ; export REDUCEDGRID=.true.
export wavecoef=62 ; export wavegrid=62
export lm=64       ; export lsoil=4         ; export MEMBER_NAMES=c00
export IDVC=3      ; export THERMODYN_ID=3  ; export SFCPRESS_ID=2 ; export SPECTRALLOOP=1
export NST_FCST=0  ; export NDSLFV=.false.
export GOCART_AER2POST=.false.
}

export_fim ()
{
export_common
export FIM_USE_NEMS=true
}

export_nmm ; export_gfs

############################################################################
# not run post
############################################################################

 if [ ${CB_arg} != post ]; then

############################################################################

############################################################################
# Clean and compile both NMMB & GFS cores, using ESMF 3.1.0rp2 library.
############################################################################

echo "Preparing model code for regression tests"
echo "Using the ESMF 3.1.0rp2 library"
printf %s "Using the ESMF 3.1.0rp2 library.   "
printf %s "Compiling model code (this will take some time)......."
cd ${PATHTR}/src

date                                     >> ${PATHRT}/RegressionTests.log
echo "Compilation ALL"                   >> ${PATHRT}/RegressionTests.log
#rm -f ../exe/NEMS.x
#gmake clean                              >  ${PATHRT}/Compile.log 2>&1
#esmf_version 3                           >> ${PATHRT}/Compile.log 2>&1
if [ ${MACHINE_ID} = c -o ${MACHINE_ID} = s -o ${MACHINE_ID} = v ]; then
  gmake nmm_gfs_gen GOCART_MODE=full     >> ${PATHRT}/Compile.log 2>&1
elif [ ${MACHINE_ID} = g -o ${MACHINE_ID} = h -o ${MACHINE_ID} = z ]; then 
  gmake nmm                              >> ${PATHRT}/Compile.log 2>&1
fi
date                                     >> ${PATHRT}/RegressionTests.log

if [ -f ../exe/NEMS.x ] ; then
  echo "   Model code Compiled";echo;echo
else
  echo "   Model code is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   Model code is NOT compiled"
  exit
fi

cd $PATHRT

####################################################################################################
#
# TEST   - Global NMM-B with pure binary input
#        - 6x5 compute  tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post ]; then

export TEST_DESCR="Compare NMMB-global results with previous trunk version"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_CNTRL
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s nmmb_hst_01_bin_0048h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s nmmb_hst_01_nio_0048h_00m_00.00s \
nmmb_rst_01_bin_0024h_00m_00.00s nmmb_rst_01_nio_0024h_00m_00.00s"
#---------------------
export_nmm
export GBRG=glob ; export WLCLK=04
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

if [ ${MACHINE_ID} = c -o ${MACHINE_ID} = s -o ${MACHINE_ID} = v ]; then
  export timing1=`grep total_integration_tim $RUNDIR/err | tail -1 | awk '{ print $5 }'`
elif [ ${MACHINE_ID} = g -o ${MACHINE_ID} = h -o ${MACHINE_ID} = z ]; then 
  export timing1=`grep total_integration_tim $RUNDIR/err | tail -1 | awk '{ print $4 }'`
fi
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

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-global NEMSIO as input file"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_NEMSIO
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0048h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0048h_00m_00.00s \
nmmb_rst_01_bin_0024h_00m_00.00s nmmb_rst_01_nio_0024h_00m_00.00s"
#---------------------
export_nmm
export GBRG=glob ; export NEMSI=true
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

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-global restart run"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REST
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmmb_hst_01_bin_0048h_00m_00.00s"
#---------------------
export_nmm
export GBRG=glob ; export RSTRT=true
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
export LIST_FILES=" \
nmmb_hst_01_bin_0048h_00m_00.00s"
#---------------------
export_nmm
export GBRG=glob ; export NEMSI=true ; export RSTRT=true
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
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s"
#---------------------
export_nmm
export GBRG=glob   ; export FCSTL=24
export INPES=$JNPm ; export JNPES=$INPm
export WLCLK=04
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
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0048h_00m_00.00s        \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0048h_00m_00.00s \
nmmb_rst_01_bin_0024h_00m_00.00s nmmb_rst_01_nio_0024h_00m_00.00s"
#---------------------
export_nmm
export GBRG=glob ; export THRD=2
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

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then
 
export TEST_DESCR="Test NMMB-global with GFS physics package "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_gfsP
export CNTL_DIR=NMMB_gfsP_glob
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s \
nmmb_rst_01_bin_0012h_00m_00.00s nmmb_rst_01_nio_0012h_00m_00.00s"
#---------------------
export_nmm
export GBRG=glob ; export gfsP=true ; export FCSTL=24
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

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="Compare NMMB-regional results with previous trunk version"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_CTL
export CNTL_DIR=NMMB_reg
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_nio_0000h_00m_00.00s \
nmmb_hst_01_bin_0012h_00m_00.00s nmmb_hst_01_nio_0012h_00m_00.00s \
nmmb_hst_01_bin_0048h_00m_00.00s nmmb_hst_01_nio_0048h_00m_00.00s \
nmmb_rst_01_bin_0024h_00m_00.00s nmmb_rst_01_nio_0024h_00m_00.00s \
fort.41  fort.42  fort.43  fort.44  fort.45  fort.46  fort.47"
#---------------------
export_nmm
export GBRG=reg ; export WLCLK=06 ; export WPREC=true
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

if [ ${MACHINE_ID} = c -o ${MACHINE_ID} = s -o ${MACHINE_ID} = v ]; then
  export timing1=`grep total_integration_tim $RUNDIR/err | tail -1 | awk '{ print $5 }'`
elif [ ${MACHINE_ID} = g -o ${MACHINE_ID} = h -o ${MACHINE_ID} = z ]; then 
  export timing1=`grep total_integration_tim $RUNDIR/err | tail -1 | awk '{ print $4 }'`
fi
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
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0012h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0012h_00m_00.00s"
#---------------------
export_nmm
export GBRG=reg ; export NEMSI=true ; export FCSTL=12
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
export LIST_FILES=" \
nmmb_hst_01_bin_0048h_00m_00.00s"
#---------------------
export_nmm
export GBRG=reg ; export RSTRT=true
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

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-regional restart run with NEMSIO file "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_RST_NIO
export CNTL_DIR=NMMB_reg
export LIST_FILES=" \
nmmb_hst_01_bin_0048h_00m_00.00s"
#---------------------
export_nmm
export GBRG=reg ; export NEMSI=true ; export RSTRT=true
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

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-regional different decomposition"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_DECOMP
export CNTL_DIR=NMMB_reg
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0012h_00m_00.00s"
#---------------------
export_nmm
export GBRG=reg    ; export FCSTL=12
export INPES=$JNPm ; export JNPES=$INPm
export WLCLK=04
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

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-regional threading "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_THREAD
export CNTL_DIR=NMMB_reg
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0048h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0048h_00m_00.00s \
nmmb_rst_01_bin_0024h_00m_00.00s nmmb_rst_01_nio_0024h_00m_00.00s"
#---------------------
export_nmm
export GBRG=reg ; export THRD=2
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

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post ]; then

export TEST_DESCR="Test NMMB-regional with GFS physics package "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_gfsP
export CNTL_DIR=NMMB_gfsP_reg
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s \
nmmb_rst_01_bin_0012h_00m_00.00s nmmb_rst_01_nio_0012h_00m_00.00s"
#---------------------
export_nmm
export GBRG=reg ; export gfsP=true ; export FCSTL=24
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------
 
fi

####################################################################################################
#
# TEST   - Regional NMM-B with selected GFS physics schemes
#        - 6x5 compute tasks / 1 thread / selected GFS physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-regional with selected GFS physics schemes "

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_SEL_PHY
export CNTL_DIR=NMMB_reg_sel_phy
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s \
nmmb_rst_01_bin_0012h_00m_00.00s nmmb_rst_01_nio_0012h_00m_00.00s"
#---------------------
export_nmm
export GBRG=reg  ; export FCSTL=24
export CONVC=sas ; export MICRO=gfs ; export TURBL=gfs
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

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-regional with precipitation adjustment on"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_PCPADJ
export CNTL_DIR=NMMB_reg_pcpadj
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0012h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0012h_00m_00.00s \
nmmb_rst_01_bin_0012h_00m_00.00s nmmb_rst_01_nio_0012h_00m_00.00s"
#---------------------
export_nmm
export GBRG=reg    ; export FCSTL=12
export PCPFLG=true ; export CPPCP=''
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Regional NMM-B writing time series
#        - 6x5 compute tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-regional writing time series"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_TIMESR
export CNTL_DIR=NMMB_reg_timesr
export LIST_FILES=" \
nmmb_hst_01_bin_0006h_00m_00.00s ts_p01_d01.bin ts_p02_d01.bin"
#---------------------
export_nmm
export GBRG=reg ; export TS='' ; export FCSTL=06
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - NMM-B static nests: Regional parent with two children and one grandchild
#        - Compute tasks - Upper parent 2x3 | Child #1 4x8 | Child #2 2x4 | Grandchild 7x10
#        - 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-regional with static nests"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_nests
export CNTL_DIR=NMMB_nests
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s \
nmmb_rst_01_bin_0012h_00m_00.00s nmmb_rst_01_nio_0012h_00m_00.00s \
nmmb_hst_02_bin_0000h_00m_00.00s nmmb_hst_02_bin_0024h_00m_00.00s \
nmmb_hst_02_nio_0000h_00m_00.00s nmmb_hst_02_nio_0024h_00m_00.00s \
nmmb_rst_02_bin_0012h_00m_00.00s nmmb_rst_02_nio_0012h_00m_00.00s \
nmmb_hst_03_bin_0000h_00m_00.00s nmmb_hst_03_bin_0024h_00m_00.00s \
nmmb_hst_03_nio_0000h_00m_00.00s nmmb_hst_03_nio_0024h_00m_00.00s \
nmmb_rst_03_bin_0012h_00m_00.00s nmmb_rst_03_nio_0012h_00m_00.00s \
nmmb_hst_04_bin_0000h_00m_00.00s nmmb_hst_04_bin_0024h_00m_00.00s \
nmmb_hst_04_nio_0000h_00m_00.00s nmmb_hst_04_nio_0024h_00m_00.00s \
nmmb_rst_04_bin_0012h_00m_00.00s nmmb_rst_04_nio_0012h_00m_00.00s"
#---------------------
export_nmm
export GBRG=nests ; export TPN=$TPNn ; export FCSTL=24
export AFFN=cpu   ; export NODE=2
export INPES=02   ; export JNPES=03  ; export WTPG=1
export WLCLK=20   ; export NCHILD=02
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - NMM-B restart static nests: Regional parent with two children and one grandchild
#        - Compute tasks - Upper parent 2x3 | Child #1 4x8 | Child #2 2x4 | Grandchild 7x10
#        - 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

if [ ${CREATE_BASELINE} = false ]; then

export TEST_DESCR="Test NMMB-regional static nests with restart"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_nest_rest
export CNTL_DIR=NMMB_nests
export LIST_FILES=" \
nmmb_hst_01_bin_0024h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s \
nmmb_hst_02_bin_0024h_00m_00.00s nmmb_hst_02_nio_0024h_00m_00.00s \
nmmb_hst_03_bin_0024h_00m_00.00s nmmb_hst_03_nio_0024h_00m_00.00s \
nmmb_hst_04_bin_0024h_00m_00.00s nmmb_hst_04_nio_0024h_00m_00.00s"
#---------------------
export_nmm
export GBRG=nests ; export TPN=$TPNn ; export FCSTL=24
export AFFN=cpu   ; export NODE=2
export INPES=02   ; export JNPES=03  ; export WTPG=1
export RSTRT=true ; export WLCLK=12  ; export NCHILD=02
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - Regional NMM-B static nests with filter
#        - Compute tasks - Upper parent 2x2 | Child #1 3x5 | Grandchild 6x7
#        - 1 thread / opnl physics / free fcst / nemsio binary input
#
#
####################################################################################################

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-regional digital filter with static nests"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_REG_FILT
export CNTL_DIR=NMMB_reg_filt
export LIST_FILES=" \
nmmb_hst_01_bin_0003h_00m_00.00s nmmb_hst_01_nio_0003h_00m_00.00s \
nmmb_hst_02_bin_0003h_00m_00.00s nmmb_hst_02_nio_0003h_00m_00.00s \
nmmb_hst_03_bin_0003h_00m_00.00s nmmb_hst_03_nio_0003h_00m_00.00s"
#---------------------
export_nmm
export GBRG=fltr  ; export TPN=64   ; export FCSTL=03
export INPES=02   ; export JNPES=02 ; export WTPG=1
export NEMSI=true ; export WLCLK=06 ; export NCHILD=01
export AFFN=cpu
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - NMM-B moving nests: Regional parent with two children and one grandchild
#        - Compute tasks - Upper parent 8x8 | Child #1 2x6 | Child #2 2x6 | Grandchild 5x6
#        - 1 thread / opnl physics / free fcst / nemsio and pure binary input
#
####################################################################################################

if [ ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test NMMB-regional with moving nests"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_mvg_nests
export CNTL_DIR=NMMB_mvg_nests
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s \
nmmb_rst_01_bin_0012h_00m_00.00s nmmb_rst_01_nio_0012h_00m_00.00s \
nmmb_hst_02_bin_0000h_00m_00.00s nmmb_hst_02_bin_0024h_00m_00.00s \
nmmb_hst_02_nio_0000h_00m_00.00s nmmb_hst_02_nio_0024h_00m_00.00s \
nmmb_rst_02_bin_0012h_00m_00.00s nmmb_rst_02_nio_0012h_00m_00.00s \
nmmb_hst_03_bin_0000h_00m_00.00s nmmb_hst_03_bin_0024h_00m_00.00s \
nmmb_hst_03_nio_0000h_00m_00.00s nmmb_hst_03_nio_0024h_00m_00.00s \
nmmb_rst_03_bin_0012h_00m_00.00s nmmb_rst_03_nio_0012h_00m_00.00s \
nmmb_hst_04_bin_0000h_00m_00.00s nmmb_hst_04_bin_0024h_00m_00.00s \
nmmb_hst_04_nio_0000h_00m_00.00s nmmb_hst_04_nio_0024h_00m_00.00s \
nmmb_rst_04_bin_0012h_00m_00.00s nmmb_rst_04_nio_0012h_00m_00.00s"
#---------------------
export_nmm
export GBRG=mnests ; export TPN=$TPNn ; export FCSTL=24
export AFFN=cpu    ; export NODE=2
export INPES=08    ; export JNPES=08  ; export WTPG=4
export NEMSI=true  ; export WLCLK=10  ; export NCHILD=02
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

cd $PATHRT

if [ ${MACHINE_ID} = z ]; then
  exit
fi

####################################################################################################
# 
# TEST   - GFS 
#        - 30 compute tasks / 1 thread 
#
####################################################################################################

if [ ${CB_arg} != nmm -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post ]; then

export TEST_DESCR="Compare GFS results with previous trunk version"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf00 sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf00 sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf00 flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
# 
# TEST   - GFS  adiabatic run
#        - 30 compute tasks / 1 thread 
#
####################################################################################################

if [ ${CB_arg} = gfs -o ${RT_FULL} = true ]; then

export TEST_DESCR="Compare GFS adiabatic results with previous trunk version"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_ADIAB
export CNTL_DIR=GFS_ADIAB
export LIST_FILES=" sigf03 sigf06 sigf12 sigf24 "
#---------------------
export_gfs
export NDAYS=1 ; export ADIAB=.true.
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS as two copies
#        - 30 compute tasks / 1 thread  2 copy
#
####################################################################################################

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test GFS with 2-copy option"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_48
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export_gfs
export NDAYS=1 ; export CP2=.true.
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS with different decomposition
#        - 58 compute tasks / 1 thread  restart
#
####################################################################################################

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test GFS different decomposition and restart"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_48
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
export TASKS=48 ; export PE1=46
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
export RUNDIR=${RUNDIR_ROOT}/GFS_16_32
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
       sigf03 sigf06 sigf12 sigf24 \
       sfcf03 sfcf06 sfcf12 sfcf24 \
       flxf03 flxf06 flxf12 flxf24"
#---------------------
export_gfs
export TASKS=16 ; export PE1=12 ; export THRD=2
export WRTGP=2  ; export NDAYS=1
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

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="GFS, 32 proc, 1 thread, no quilt, output every 4 timestep"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16_32
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
        sigf03 sigf06 sigf12 sigf24 sigf48 \
        sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
        flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
export NSOUT=4 ; export QUILT=.false.
export PE1=32  ; export WTPG=1
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

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

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
export_gfs
export TASKS=1 ; export PE1=1 ; export WTPG=1
export QUILT=.false. ; export NDAYS=1 ; export WLCLK=20
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

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

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
export_gfs
export TASKS=1 ; export PE1=1 ; export WTPG=1
export NSOUT=1 ; export WLCLK=20
export QUILT=.false.
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

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test GFS, 16 proc, 2 threads,no quilt, output every 2 time steps"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16_48_NOQUILT_NSOUT
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export_gfs
export TASKS=16 ; export PE1=16  ; export WTPG=1
export THRD=2   ; export NSOUT=2 ; export NDAYS=1
export QUILT=.false.
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS with multiple tasks and no quilting
#        - 48 tasks / 1 thread
#
####################################################################################################

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="Test GFS, 48 proc, 1 thread, no quilt"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16_48_NOQUILT_NSOUT
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
export TASKS=48 ; export PE1=46 ; export WTPG=1
export NSOUT=1  ; export QUILT=.false.
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

if [ ${CB_arg} != nmm -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="GFS,32 total proc (tasks), 1 thread, quilt, digital filter on reduced grid"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_16_dfi
export CNTL_DIR=GFS_DFI_REDUCEDGRID
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf03 flxf06 flxf12 flxf24 flxf48" 
#---------------------
export_gfs
export FDFI=3
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

export TEST_DESCR="GFS,16 total proc (tasks), 2 thread, quilt,2x2 wrt pe, digital filter on reduced grid"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16_1_dfi
export CNTL_DIR=GFS_DFI_REDUCEDGRID
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export_gfs
export TASKS=16 ; export PE1=12 ; export WRTGP=2 ; export THRD=2
export NDAYS=1  ; export FDFI=3 ; export CP2=.true.
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS digital filter
#        - 1pe nsout=0
#
####################################################################################################

if [ ${CREATE_BASELINE} = false -a ${RT_FULL} = true ]; then

export TEST_DESCR="GFS,1 proc, no quilt, digital filter on reduced grid"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16_1_dfi
export CNTL_DIR=GFS_DFI_REDUCEDGRID
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 sigf48 \
        sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
        flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
export TASKS=1 ; export PE1=1    ; export WTPG=1
export FDFI=3  ; export WLCLK=20 ; export QUILT=.false.
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

if [ ${CB_arg} != nmm -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="GFS, use the OPAC climo scheme for SW and LW"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_OPAC
export CNTL_DIR=GFS_OPAC
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export_gfs
export IAER=11 ; export NDAYS=1
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS 
#        - NDSL 12 compute tasks / 2 thread
#
####################################################################################################

if [ ${CB_arg} != nmm -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post ]; then

export TEST_DESCR="GFS, 16tasks, 2threads, quilt, dfi3hr, reduced grid, NDSL"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_dfi_ndsl
export CNTL_DIR=GFS_DFI_REDUCEDGRID_NDSL
export LIST_FILES=" \
        sigf03 sigf06 sigf12 sigf24 \
        sfcf03 sfcf06 sfcf12 sfcf24 \
        flxf03 flxf06 flxf12 flxf24"
#---------------------
export_gfs
export TASKS=16 ; export PE1=12 ; export THRD=2 ; export WRTGP=2
export FDFI=3   ; export NDSLFV=.true.
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

#
####################################################################################################
#
# TEST   - GFS hyb 2 loop with digital filter
#        - 12 compute tasks / 2 thread ,2 WrtGrp x 2 WrtPePerGrp
#
####################################################################################################

if [ ${CB_arg} != nmm -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="GFS,16 total proc (tasks), 2 thread, quilt,2x2 wrt pe, HYB 2loop digital filter on reduced grid"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_hyb_2loop
export CNTL_DIR=GFS_DFI_hyb_2loop
export LIST_FILES=" \
	sigf03 sigf06 sigf12 sigf24 \
	sfcf03 sfcf06 sfcf12 sfcf24 \
	flxf03 flxf06 flxf12 flxf24"
#---------------------
export_gfs
export TASKS=16 ; export PE1=12 ; export WRTGP=2 ; export THRD=2
export FDFI=3   ; export CP2=.true.
export IDVC=2   ; export THERMODYN_ID=0  ; export SFCPRESS_ID=0 ; export SPECTRALLOOP=2
export NDAYS=1
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi
#
###################################################################################################
#
#
# TEST   - GFS hyb 2 loop with digital filter restart
#        - 12 compute tasks / 2 thread ,2 WrtGrp x 2 WrtPePerGrp
#
###################################################################################################
#

if [ ${CB_arg} != nmm -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="GFS,16 total proc (tasks), 2 thread, quilt,2x2 wrt pe, HYB 2loop digital filter
on reduced grid"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_hyb_2loop
export CNTL_DIR=GFS_DFI_hyb_2loop
export LIST_FILES=" \
        sigf03 sigf06 sigf12 sigf24 sigf48 \
        sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
        flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
export TASKS=16 ; export PE1=12 ; export WRTGP=2 ; export THRD=2
export FDFI=3   ; export CP2=.true.
export IDVC=2   ; export THERMODYN_ID=0  ; export SFCPRESS_ID=0 ; export SPECTRALLOOP=2
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS hyb 2 loop with digital filter, with nst
#        - 12 compute tasks / 2 thread ,2 WrtGrp x 2 WrtPePerGrp
#
####################################################################################################

if [ ${CB_arg} != nmm -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post ]; then

export TEST_DESCR="GFS,16 total proc (tasks), 2 thread, quilt,2x2 wrt pe, HYB 2loop digital filter on reduced grid with nst"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_hyb_2loop_nst
export CNTL_DIR=GFS_DFI_hyb_2loop_nst
export LIST_FILES=" \
        sigf03 sigf06 sigf12 sigf24 sigf48 \
        sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
        flxf03 flxf06 flxf12 flxf24 flxf48 \
        nstf03 nstf06 nstf12 nstf24 nstf48"
#---------------------
export_gfs
export TASKS=16 ; export PE1=12 ; export THRD=2 ; export WRTGP=2
export FDFI=3   ; export NST_FCST=1 ; export CP2=.true.
export IDVC=2   ; export THERMODYN_ID=0  ; export SFCPRESS_ID=0 ; export SPECTRALLOOP=2
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

####################################################################################################
#
# TEST   - GFS digital filter  HYb 1loop
#        - 12 compute tasks / 2 thread ,2 WrtGrp x 2 WrtPePerGrp
#
####################################################################################################

if [ ${CB_arg} != nmm -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post ]; then

export TEST_DESCR="GFS,16 total proc (tasks), 2 thread, quilt,2x2 wrt pe, HYB 1loop digital filter on reduced grid"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_16_dfi_hyb_1loop
export CNTL_DIR=GFS_DFI_REDUCEDGRID_HYB
export LIST_FILES=" \
	sigf00 sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf00 sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf00 flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
export TASKS=16 ; export PE1=12 ; export THRD=2 ; export WRTGP=2
export FDFI=3 ; export CP2=.true.
export IDVC=2 ; export THERMODYN_ID=0  ; export SFCPRESS_ID=0 ; export SPECTRALLOOP=1
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

if [ ${CB_arg} != nmm -a ${CB_arg} != gen -a ${CB_arg} != fim -a ${CB_arg} != post ]; then

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
export_gfs
export GEFS_ENSEMBLE=1
export TASKS=64 ; export WLCLK=20
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

cd $PATHRT

####################################################################################################
#
# TEST   - GEN, 16 PEs, 1 node.
#        - 1 members.
#
####################################################################################################

if [ ${CB_arg} != nmm -a ${CB_arg} != gfs -a ${CB_arg} != fim -a ${CB_arg} != post ]; then

export TEST_DESCR="GEN, 1 members."

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GEN_Run_m1
#---------------------
export_common
export TASKS=16 ; export WLCLK=02
#---------------------
  ./rt_gen.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

cd $PATHRT

####################################################################################################
#
# TEST   - Concurrency GEN, 64PEs, 1 node.
#        - 4 members.
#
####################################################################################################

if [ ${CB_arg} != nmm -a ${CB_arg} != gfs -a ${CB_arg} != fim -a ${CB_arg} != post -a ${RT_FULL} = true ]; then

export TEST_DESCR="Concurrency GEN, 4 members."

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GEN_Concurrency_Run_m4
#---------------------
export_common
export GEN_ENSEMBLE=1
export TASKS=64 ; export WLCLK=02
#---------------------
  ./rt_gen.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

################################################################################################
# run nmm gfs gen without post 
################################################################################################
  fi


################################################################################################
  if [ ${RT_FULL} = true ]; then
#
################################################################################################

#########################################################################
#########################################################################
#
# FIM-only tests
#
#########################################################################
#########################################################################

################################################################################################
if [ ${CB_arg} != nmm -a ${CB_arg} != gfs -a ${CB_arg} != gen -a ${CB_arg} != post ]; then
#
################################################################################################

#########################################################################
#########################################################################
#
# Clean and compile only FIM core
#
#########################################################################
#########################################################################

echo "Preparing model code for regression tests"
printf %s "Compiling model code (this will take some time)......."
cd ${PATHTR}/src

date                                     >> ${PATHRT}/RegressionTests.log
echo "Compilation only FIM"              >> ${PATHRT}/RegressionTests.log
rm -f ../exe/NEMS.x
gmake clean                              >> ${PATHRT}/Compile.log 2>&1
esmf_version 3                           >> ${PATHRT}/Compile.log 2>&1
gmake fim                                >> ${PATHRT}/Compile.log 2>&1
date                                     >> ${PATHRT}/RegressionTests.log

if [ -f ../exe/NEMS.x ] ; then
  echo "   Model code Compiled";echo;echo
else
  echo "   Model code is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   Model code is NOT compiled"
  exit
fi

cd $PATHRT


####################################################################################################
#
# TEST   - FIM
#        - 40 compute tasks / 1 thread, 1 node.  
#
####################################################################################################

export TEST_DESCR="Compare FIM results with previous trunk version, only FIM compiled"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export TASKS=40
export GLVL=4
export RUNDIR=${RUNDIR_ROOT}/FIM_G${GLVL}L38P${TASKS}_24hr
export CNTL_DIR=FIM_G${GLVL}L38_24hr
export LIST_FILES=" \
        fim_out_2D__000000hr fim_out_2D__000003hr fim_out_2D__000006hr \
        fim_out_2D__000009hr fim_out_2D__000012hr fim_out_2D__000015hr \
        fim_out_2D__000018hr fim_out_2D__000021hr fim_out_2D__000024hr \
        fim_out_da3D000024hr fim_out_db3D000024hr fim_out_dp3D000024hr \
        fim_out_hgtP000024hr fim_out_mp3D000024hr fim_out_oz3D000024hr \
        fim_out_ph3D000024hr fim_out_pr3D000024hr fim_out_qv3D000024hr \
        fim_out_qw3D000024hr fim_out_rh3D000024hr fim_out_rp3P000024hr \
        fim_out_td3D000024hr fim_out_th3D000024hr fim_out_tk3D000024hr \
        fim_out_tmpP000024hr fim_out_up3P000024hr fim_out_us3D000024hr \
        fim_out_vo3D000024hr fim_out_vp3P000024hr fim_out_vs3D000024hr \
        fim_out_ws3D000024hr"
#---------------------
export_fim
export WLCLK=15
#---------------------
  ./rt_fim.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi  # ${CB_arg} "if" statement



#
#   Now test post options for nmm and gfs
#     - nmm_gfs_gen_post GOCART_MODE=full (ESMF4)
#
################################################################################################
# Clean and compile both NMMB & GFS cores, using ESMF 3.1.0rp2 and POST library.
################################################################################################

echo "Preparing model code for regression tests"
echo "Using the ESMF 3.1.0rp2 library"
printf %s "Using the ESMF 3.1.0rp2 and POST library.   "
printf %s "Compiling model code (this will take some time)......."
cd ${PATHTR}/src

date                                     >> ${PATHRT}/RegressionTests.log
echo "Compilation ALL"                   >> ${PATHRT}/RegressionTests.log
rm -f ../exe/NEMS.x
gmake clean                              >  ${PATHRT}/Compile.log 2>&1
esmf_version 3                           >> ${PATHRT}/Compile.log 2>&1
gmake nmm_gfs_gen_post GOCART_MODE=full  >> ${PATHRT}/Compile.log 2>&1
date                                     >> ${PATHRT}/RegressionTests.log

cd $PATHRT
#################################################################################################
#
# TEST   - Regional NMM-B with pure binary input and post
#        - 6x5 compute tasks / 1 thread / opnl physics / free fcst / pure binary input
#
#################################################################################################

if [ ${CREATE_BASELINE} = true -a ${CB_arg} != gfs -a ${CB_arg} != fim -o ${RT_FULL} = true -a $argn -eq 2 ]; then

export TEST_DESCR="NMMB-regional run with post on quilt"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_reg_post
export CNTL_DIR=NMMB_reg_post
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0048h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0048h_00m_00.00s \
nmmb_rst_01_bin_0024h_00m_00.00s nmmb_rst_01_nio_0024h_00m_00.00s \
BGDAWP.GrbF48 BGRD3D.GrbF48 BGRDSF.GrbF48"
#---------------------
export_nmm
export GBRG=reg ; export WLCLK=10 ; export WPREC=true
export WRITE_DOPOST=.true.
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi

cd $PATHRT
#################################################################################################
#
# TEST   - GFS_post
#        - 30 compute tasks / 1 thread
#
#################################################################################################

if [ ${CREATE_BASELINE} = true -a ${CB_arg} != nmm -a ${CB_arg} != fim -o ${RT_FULL} = true -a $argn -eq 2 ]; then

export TEST_DESCR="Compare GFS results with previous trunk version"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_DFI_POST
export CNTL_DIR=GFS_DFI_POST
export LIST_FILES=" \
        sigf00 sigf03 sigf06 sigf12 sigf24 sigf48 \
        sfcf00 sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
        flxf00 flxf03 flxf06 flxf12 flxf24 flxf48 \
        GFSPRS.GrbF03 GFSPRS.GrbF06 GFSPRS.GrbF12 \
        GFSPRS.GrbF24 GFSPRS.GrbF48"
#---------------------
export_gfs
#
export FDFI=3
export WRITE_DOPOST=.true.
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

fi


####################################################################
fi #end regression test for POST
#####################################################################




####################################################################################################
#
#   Now test different compilation options (only as part of FULL test)
#     - nmm_gfs_gen GOCART_MODE=full (ESMF4)
#     - nmm
#     - gfs
#     - nmm (traps ON)
#
####################################################################################################

if [ ${RT_FULL} = true -a $argn -eq 2 ]; then

#########################################################################
#########################################################################
#
# Clean and compile both NMMB & GFS cores, using ESMF 4.0.0rp2 library.
#
#########################################################################
#########################################################################

echo "Preparing model code for regression tests"
echo "Using the ESMF 4.0.0rp2 library"
printf %s "Using the ESMF 4.0.0rp2 library.   "
printf %s "Compiling model code (this will take some time)......."
cd ${PATHTR}/src

date                                     >> ${PATHRT}/RegressionTests.log
echo "Compilation ESMF 4"                >> ${PATHRT}/RegressionTests.log
rm -f ../exe/NEMS.x
gmake clean                              >> ${PATHRT}/Compile.log 2>&1
esmf_version 4                           >> ${PATHRT}/Compile.log 2>&1
gmake nmm_gfs_gen                        >> ${PATHRT}/Compile.log 2>&1
date                                     >> ${PATHRT}/RegressionTests.log

if [ -f ../exe/NEMS.x ] ; then
  echo "   Model code Compiled";echo;echo
else
  echo "   Model code is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   Model code is NOT compiled"
  exit
fi

cd $PATHRT

####################################################################################################
#
# TEST   - Global NMM-B with pure binary input
#        - 6x5 compute  tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

export TEST_DESCR="Compare NMMB-global results with previous trunk version ESMF4"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_CNTRL_ESMF4
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s \
nmmb_rst_01_bin_0024h_00m_00.00s nmmb_rst_01_nio_0024h_00m_00.00s "
#---------------------
export_nmm
export GBRG=glob ; export FCSTL=24 ; export WLCLK=02
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
# 
# TEST   - GFS 
#        - 30 compute tasks / 1 thread 
#
####################################################################################################

export TEST_DESCR="Compare GFS results with previous trunk version ESMF4"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_ESMF4
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf00 sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf00 sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf00 flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
#
# TEST   - Concurrency GEFS
#        - 4 members, every 6 hours, couple and add stochastic perturbations, T190L28.
#
####################################################################################################

export TEST_DESCR="Concurrency GEFS, stochastic perturbations, 4 members, T190L28. ESMF4"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GEFS_Concurrency_Run_ESMF4
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
export_gfs
export GEFS_ENSEMBLE=1
export TASKS=64 ; export WLCLK=20
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
#
# TEST   - Concurrency GEN, 64PEs, 1 node.
#        - 4 members.
#
####################################################################################################

export TEST_DESCR="Concurrency GEN, 4 members. ESMF4"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GEN_Concurrency_Run_m4_ESMF4
#---------------------
export_common
export GEN_ENSEMBLE=1
export TASKS=64 ; export WLCLK=02
#---------------------
  ./rt_gen.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------


#########################################################################
#########################################################################
#
# Clean and compile both NMMB & GFS cores, using ESMF 5.1.0 library.
#
#########################################################################
#########################################################################

echo "Preparing model code for regression tests"
echo "Using the ESMF 5.1.0 library"
printf %s "Using the ESMF 5.1.0 library.   "
printf %s "Compiling model code (this will take some time)......."
cd ${PATHTR}/src

date                                     >> ${PATHRT}/RegressionTests.log
echo "Compilation ESMF 5"                >> ${PATHRT}/RegressionTests.log
rm -f ../exe/NEMS.x
gmake clean                              >> ${PATHRT}/Compile.log 2>&1
esmf_version 5                           >> ${PATHRT}/Compile.log 2>&1
gmake nmm_gfs_gen                        >> ${PATHRT}/Compile.log 2>&1
date                                     >> ${PATHRT}/RegressionTests.log

if [ -f ../exe/NEMS.x ] ; then
  echo "   Model code Compiled";echo;echo
else
  echo "   Model code is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   Model code is NOT compiled"
  exit
fi

cd $PATHRT

####################################################################################################
#
# TEST   - Global NMM-B with pure binary input
#        - 6x5 compute  tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

export TEST_DESCR="Compare NMMB-global results with previous trunk version ESMF5"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_CNTRL_ESMF5
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s \
nmmb_rst_01_bin_0024h_00m_00.00s nmmb_rst_01_nio_0024h_00m_00.00s "
#---------------------
export_nmm
export GBRG=glob ; export FCSTL=24 ; export WLCLK=02
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
# 
# TEST   - GFS 
#        - 30 compute tasks / 1 thread 
#
####################################################################################################

export TEST_DESCR="Compare GFS results with previous trunk version ESMF5"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_ESMF5
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
	sigf00 sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf00 sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf00 flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
#
# TEST   - Concurrency GEFS
#        - 4 members, every 6 hours, couple and add stochastic perturbations, T190L28.
#
####################################################################################################

export TEST_DESCR="Concurrency GEFS, stochastic perturbations, 4 members, T190L28. ESMF5"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GEFS_Concurrency_Run_ESMF5
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
export_gfs
export GEFS_ENSEMBLE=1
export TASKS=64 ; export WLCLK=20
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
#
# TEST   - Concurrency GEN, 64PEs, 1 node.
#        - 4 members.
#
####################################################################################################

export TEST_DESCR="Concurrency GEN, 4 members. ESMF5"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GEN_Concurrency_Run_m4_ESMF5
#---------------------
export_common
export GEN_ENSEMBLE=1
export TASKS=64 ; export WLCLK=02
#---------------------
  ./rt_gen.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

#########################################################################
#########################################################################
#
# Clean and compile both NMMB & GFS cores, 
# using ESMF 5.2.0rp1 library.
#
#########################################################################
#########################################################################

echo "Preparing model code for regression tests"
echo "Using the ESMF 5.2.0rp1 library"
printf %s "Using the ESMF 5.2.0rp1 library.   "
printf %s "Compiling model code (this will take some time)......."
cd ${PATHTR}/src

date                                     >> ${PATHRT}/RegressionTests.log
echo "Compilation ESMF 5.2.0rp1"         >> ${PATHRT}/RegressionTests.log
rm -f ../exe/NEMS.x
gmake clean                              >> ${PATHRT}/Compile.log 2>&1
esmf_version 5.2                         >> ${PATHRT}/Compile.log 2>&1
gmake gfs                                >> ${PATHRT}/Compile.log 2>&1
date                                     >> ${PATHRT}/RegressionTests.log

if [ -f ../exe/NEMS.x ] ; then
  echo "   Model code Compiled";echo;echo
else
  echo "   Model code is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   Model code is NOT compiled"
  exit
fi

cd $PATHRT

####################################################################################################
# 
# TEST   - GFS 
#        - 30 compute tasks / 1 thread 
#
####################################################################################################

export TEST_DESCR="Compare GFS results with previous trunk version ESMF5.2.0rp1"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_ESMF5.2.0rp1
export CNTL_DIR=GFS_NODFI_5.2.0rp1
export LIST_FILES=" \
	sigf00 sigf03 sigf06 sigf12 sigf24 sigf48 \
	sfcf00 sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
	flxf00 flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
export WLCLK=120
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

####################################################################################################
#
# TEST   - Concurrency GEFS
#        - 4 members, every 6 hours, couple and add stochastic perturbations, T190L28.
#
####################################################################################################

export TEST_DESCR="Concurrency GEFS, stochastic perturbations, 4 members, T190L28. ESMF5.2.0rp1"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GEFS_Concurrency_Run_ESMF5.2.0rp1
export CNTL_DIR=GEFS_m4_5.2.0rp1
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
export_gfs
export GEFS_ENSEMBLE=1
export TASKS=64 ; export WLCLK=158
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------

#########################################################################
#########################################################################
#
# Clean and compile only NMMB core
#
#########################################################################
#########################################################################

echo "Preparing model code for regression tests"
printf %s "Compiling model code (this will take some time)......."
cd ${PATHTR}/src

date                                     >> ${PATHRT}/RegressionTests.log
echo "Compilation only NMM"              >> ${PATHRT}/RegressionTests.log
rm -f ../exe/NEMS.x
gmake clean                              >> ${PATHRT}/Compile.log 2>&1
esmf_version 3                           >> ${PATHRT}/Compile.log 2>&1
gmake nmm                                >> ${PATHRT}/Compile.log 2>&1
date                                     >> ${PATHRT}/RegressionTests.log

if [ -f ../exe/NEMS.x ] ; then
  echo "   Model code Compiled";echo;echo
else
  echo "   Model code is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   Model code is NOT compiled"
  exit
fi

cd $PATHRT

####################################################################################################
#
# TEST   - Global NMM-B with pure binary input
#        - 6x5 compute  tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

export TEST_DESCR="Compare NMMB-global results with previous trunk version, only NMM compiled"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_CNTRL_NMM_only
export CNTL_DIR=NMMB_glob
export LIST_FILES=" \
nmmb_hst_01_bin_0000h_00m_00.00s nmmb_hst_01_bin_0024h_00m_00.00s \
nmmb_hst_01_nio_0000h_00m_00.00s nmmb_hst_01_nio_0024h_00m_00.00s \
nmmb_rst_01_bin_0024h_00m_00.00s nmmb_rst_01_nio_0024h_00m_00.00s "
#---------------------
export_nmm
export GBRG=glob ; export FCSTL=24 ; export WLCLK=02
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------


#########################################################################
#########################################################################
#
# Clean and compile only NMMB core with TRAPS turned on
#
#########################################################################
#########################################################################

echo "Preparing model code for regression tests"
printf %s "Compiling model code (this will take some time)......."
cd ${PATHTR}/src

date                                     >> ${PATHRT}/RegressionTests.log
echo "Compilation NMM with TRAPS on"     >> ${PATHRT}/RegressionTests.log
rm -f ../exe/NEMS.x
gmake clean                              >> ${PATHRT}/Compile.log 2>&1
esmf_version traps_on                    >> ${PATHRT}/Compile.log 2>&1
gmake nmm                                >> ${PATHRT}/Compile.log 2>&1
date                                     >> ${PATHRT}/RegressionTests.log

if [ -f ../exe/NEMS.x ] ; then
  echo "   Model code Compiled";echo;echo
else
  echo "   Model code is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   Model code is NOT compiled"
  exit
fi

cd $PATHRT

####################################################################################################
#
# TEST   - Global NMM-B with pure binary input
#        - 6x5 compute  tasks / 1 thread / opnl physics / free fcst / pure binary input
#
####################################################################################################

export TEST_DESCR="Compare NMMB-global results with previous trunk version, TRAPS on"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/NMM_CNTRL_TRAPS_on
export CNTL_DIR=NMMB_glob
export LIST_FILES=" nmmb_hst_01_bin_0000h_00m_00.00s "
#---------------------
export_nmm
export GBRG=glob ; export FCSTL=12 ; export WLCLK=02
#---------------------
  ./rt_nmm.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------


#########################################################################
#########################################################################
#
# Clean and compile only GFS core
#
#########################################################################
#########################################################################

echo "Preparing model code for regression tests"
printf %s "Compiling model code (this will take some time)......."
cd ${PATHTR}/src

date                                     >> ${PATHRT}/RegressionTests.log
echo "Compilation only GFS"              >> ${PATHRT}/RegressionTests.log
rm -f ../exe/NEMS.x
gmake clean                              >> ${PATHRT}/Compile.log 2>&1
esmf_version 3                           >> ${PATHRT}/Compile.log 2>&1
gmake gfs                                >> ${PATHRT}/Compile.log 2>&1
date                                     >> ${PATHRT}/RegressionTests.log

if [ -f ../exe/NEMS.x ] ; then
  echo "   Model code Compiled";echo;echo
else
  echo "   Model code is NOT compiled" >> ${PATHRT}/RegressionTests.log
  echo "   Model code is NOT compiled"
  exit
fi

cd $PATHRT


####################################################################################################
#
# TEST   - GFS
#        - 30 compute tasks / 1 thread
#
####################################################################################################

export TEST_DESCR="Compare GFS results with previous trunk version, only GFS compiled"

#---------------------
(( TEST_NR=TEST_NR+1 ))
export RUNDIR=${RUNDIR_ROOT}/GFS_32_GFS_only
export CNTL_DIR=GFS_NODFI
export LIST_FILES=" \
        sigf00 sigf03 sigf06 sigf12 sigf24 sigf48 \
        sfcf00 sfcf03 sfcf06 sfcf12 sfcf24 sfcf48 \
        flxf00 flxf03 flxf06 flxf12 flxf24 flxf48"
#---------------------
export_gfs
#---------------------
  ./rt_gfs.sh
  if [ $? = 2 ]; then exit ; fi
#---------------------




####################################################################################################
fi # end FULL test
####################################################################################################




####################################################################################################
# Finalize
####################################################################################################

rm -f err out configure_file* nmm_ll nmm_msub nmm_run gfs_fcst_run gfs_ll gen_fcst_run gen_ll fim_fcst_run fim_ll

cd ${PATHTR}/src
gmake clean          > /dev/null 2>&1
esmf_version 3       > /dev/null 2>&1

rm -rf ${RUNDIR_ROOT}

date >> ${PATHRT}/RegressionTests.log

banner REGRESSION TEST WAS SUCCESSFUL


exit
