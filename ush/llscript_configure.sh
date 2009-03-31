#! /bin/ksh

################################################################################
# Script name:         prepare
# Script description:  Sets up Load Leveler commands for serial batch jobs
#                      in /scripts needed to run the WRF-NMM launcher
#
# Authors:    Ed Colon
################################################################################
#
#  Used to set the following Load Leveler commands in serial batch jobs in the
#  /scripts subdirectory
#
#
#-- *** IMPORTANT *** IMPORTANT *** IMPORTANT *** IMPORTANT ***
#
#   Make sure that the settings used here are consistent with those set in prelaunch!!!
#

###############################################################################################################

# ------ NCEP IBM SP ----------------
cd ..
SRCDIR="$PWD"
cd ush
CORE1='core: nmm'
CORE2='core: gfs'
GLOBAL='global:    true'
RESTART='restart:        true'
for FILE     #loop over files pointed by the for variable FILE
do
    FOUND=`grep "$CORE1" $FILE`
    if [[ $? = 0 ]]
    then
    NMM_CORE=1
    GFS_CORE=0
    fi
    FOUND=`grep "$CORE2" $FILE`
    if [[ $? = 0 ]]
    then
    NMM_CORE=0
    GFS_CORE=1
    fi
    FOUND=`grep "$GLOBAL" $FILE`
    if [[ $? = 0 ]]
    then
     BOCO=" "
    else
     BOCO="cp \$DATADIR/boco* ."
    fi
    FOUND=`grep "$RESTART" $FILE`
    if [[ $? = 0 ]]
    then
     MAIN_INPUT=${RESTART}
     IC="restart_file"
    else
     MAIN_INPUT=${MAIN_INPUT_FILENAME}
     IC="main_input_filename"
    fi
done

#=======================================================================

if [ ${NMM_CORE} -eq 1 ]; then

 #####################################NMM SUBMISSION SCRIPT OPTIONS######################################

  JOB_TYPE=parallel
  TASK_PER_NODE=32
  NODES=1
  NODE_USAGE=not_shared
  CONSUMABLE_MEMORY=500MB
  CLASS=debug
  OUTPUT=/meso/save/wx22ec/NEMS/exp/nmm
  ERROR=/meso/save/wx22ec/NEMS/exp/gfs
  WALL_CLOCK_TIME="00\:10\:00"
  ACCT_NO=NAM-T2O
  MP_SINGLE_THREAD=yes
  THREADS=1
  RUNDIR=/ptmp/\${USER}/test
  DATADIR=/stmp/wx20py/nmmb_init
  MAIN_INPUT_FILENAME=test_input_umo_regional.d01
  main_input_file=main_input_filename
  GWD=GWD.bin
  RESTART_FILE=/ptmp/wx22ec/nmm.1/nmm_b_restart.024
  ETA=/nwprod/fix/nam_micro_lookup.dat
  RRT=/meso/save/wx22ec/test
  TBL=/meso/save/wx22ec/test
  TIMING_FILE=nmm.0.test
  HYBRID_LAUNCH=/u/wx15ja/bin/hybrid_launch
  EXEDIR=$SRCDIR/exe

#####################################NMM SUBMISSION SCRIPT OPTIONS######################################

 echo "Setup for NMM core"

cat ll_templates/nmm_submission_script_tmp | sed s:_OUTPUT_:${OUTPUT}:g \
	     | sed s:_ERROR_:${ERROR}:g \
             | sed s:_JOB_TYPE_:${JOB_TYPE}:g \
             | sed s:_TASK_PER_NODE_:${TASK_PER_NODE}:g \
             | sed s:_NODES_:${NODES}:g \
	     | sed s:_NODE_USAGE_:${NODE_USAGE}:g \
	     | sed s:_CONSUMABLE_MEMORY_:${CONSUMABLE_MEMORY}:g \
	     | sed s:_CLASS_:${CLASS}:g \
	     | sed s:_WALL_CLOCK_TIME_:${WALL_CLOCK_TIME}:g \
	     | sed s:_ACCT_NO_:${ACCT_NO}:g \
	     | sed s:_MP_SINGLE_THREAD_:${MP_SINGLE_THREAD}:g \
	     | sed s:_THREADS_:${THREADS}:g \
	     | sed s:_RUNDIR_:${RUNDIR}:g \
             | sed s:_DATADIR_:${DATADIR}:g \
	     | sed s:_SRCDIR_:${SRCDIR}:g \
	     | sed s:_MAIN_INPUT_FILENAME_:${MAIN_INPUT_FILENAME}:g \
             | sed s:_main_input_file_:${IC}:g \
	     | sed s:_CONFIG_FILE_:${FILE}:g \
             | sed s:_BOCO_:"${BOCO}":g \
             | sed s:_ETA_:${ETA}:g \
	     | sed s:_RRT_:${RRT}:g \
	     | sed s:_TBL_:${TBL}:g \
             | sed s:_TIMING_FILE_:${TIMING_FILE}:g \
             | sed s:_HYBRID_LAUNCH_:${HYBRID_LAUNCH}:g \
	     | sed s:_EXEDIR_:${EXEDIR}:g \
             > ../job/run_nmm_script

elif [ ${GFS_CORE} -eq 1 ]; then

#####################################GFS PREP SUBMISSION SCRIPT OPTIONS##################################

  JOB_NAME1=prep_atmos
  JOB_TYPE=parallel
  OUTPUT1=../exp/gfs/out_prep
  ERROR1=../exp/gfs/err_prep
  ACCT_NO=GFS-T2O
  GROUP=dev
  TASK_PER_NODE1=1
  EXPERIMENT_NAME=gfs.1
  NODE1=1
  NODE_USAGE=not_shared
  CONSUMABLE_MEMORY=4000MB
  CLASS1=debug
  GROUP=dev
  WALL_CLOCK_TIME1="00\:10\:00"
  CDATE=2007070200
  SFCPRESS_ID=2
  THERMODYN_ID=3
  HYBRID=NO
  GEN_COORD_HYBRID=YES
  WAVE=62
  LM=64
  LSOIL=4
  ENS_NUM=1
  TASKS1=1
  PE1=1
  MP_SHARED_MEMORY=yes
  MP_STDOUTMODE=ordered
  MP_SINGLE_THREAD=yes
  MEMORY_AFFINITY=MCM
  THREADS1=1
  SPECTRAL_LOOP=1
  XLSMPOPTS=parthds=\$THREADS1:spins=0:yields=0:stack=512000000
  SPINLOOPTIME=500
  YIELDLOOPTIME=500
  AIXTHREAD_SCOPE=S
  MALLOCMULTIHEAP=true
  NSCDIR=/global/noscrub
  TOPDIR=/global/save
  DUMPDIR=/meso/noscrub/wx22ec
  MP_COREFILE_FORMAT=lite
  RUNDIR=/ptmp/\${USER}/test
  SCRIPTS=${SCRIPTS:-$TOPDIR/wx23sm/gsm/scripts}
  COLD_SFC=${cold_sfc:-NO}
  NTRC=3  ;  VARID=21  ;  NUMCLD=1
  CHANGE_RES=YES
  IC_TYPE=prx
  NTRAC=${NTRAC:-3}
  NTOZ=${NTOZ:-2}
  NTCW=${NTCW:-3}
  NCLD=${NCLD:-1}
  NMTVR=${NMTVR:-14}
  LSM=${lsm:-1}
  SIGLEVEL1=\$TOPDIR/wx23sm/2jif/f2006/fix/global_hyblev.l\$lm.txt
  SIGLEVEL2=\$TOPDIR/wx23hh/00wkgfs/fix/global_hyblev3.ipa\$Apercent.txt
  Apercent=\${Apercent:-100}
  IDVM1=1
  IDSL1=1
  IDVC1=2
  NVCOORD1=2
  IVSSFC1=200509
  IDVC2=3
  IDSL2=2
  IVSSFC2=200509
  IDVM2=\$THERMODYN_ID\$SFCPRESS_ID
  IVSSIG=200509
  NVCOORD2=3
  LATCH=48
  IDVC3=1
  IVSSFC3=200509
  FIXGLOBAL=/nwprod/fix
  CHGRESDIR=\$TOPDIR/wx23hh/00wkgfs/src/global_chgres-moorthi.fd
  CHGRESEXEC=\$CHGRESDIR/global_chgres
  CHGRESSH1=\$TOPDIR/wx23hh/00wkgfs/ush/global_chgres-moorthi.sh
  LANDICE_OPT=2
  CLIMO_FIELDS_OPT=2
  OROGRAPHY=\$TOPDIR/wx23ys/fix.all/global_orography.t\${wave}.grb
  MTNVAR=\$TOPDIR/wx23my/fix.all/global_mtnvar.t\${wave}.f77
  CDUMP=gfs
  CHGRESVARS1=\"ntrac=\$ntrc,idvt=\$varid,ncldt=\$numcld,idvc=\$IDVC,IVSSIG=\$ivssig,NVCOORD=\$nvcoord,IVSSFC=\$ivssfc,idvm=\$IDVM,idsl=\$IDSL,OUTTYP=1,\"
  SIGLEVEL3=${SIGLEVEL:-/nwprod/fix/global_siglevel.l\$lm.txt}
  LONSPERLAT=\$TOPDIR/wx20mi/newres/global_lonsperlat.t\$JCAP.txt
  CHGRESSH2=${CHGRESSH:-/nwprod/ush/global_chgres.sh}
  CHGRESVARS2=\$CHGRESVARS"\"RI=\$RIlist,CPI=\$CPIlist\""

 #####################################GFS RUN SUBMISSION SCRIPT OPTIONS##################################

  JOB_NAME2=run_atmos
  OUTPUT2=../exp/gfs/out_run
  ERROR2=../exp/gfs/err_run
  TASK_PER_NODE2=32
  NODE2=1
  TASKS2=$(($TASK_PER_NODE2*$NODE2))
  PE2=${TASKS2}
  ADIAB=.false.
  REDUCED_GRID=.false.
  CLASS2=dev
  WALL_CLOCK_TIME2="02\:30\:00"
  THREADS2=2
  FHROT=0
  BIND_TASKS=no
  NSCDIR=/global/noscrub
  TOPDIR=/global/save
  DUMPDIR=/meso/noscrub/wx22ec
  MP_COREFILE_FORMAT=lite
  MP_LABELIO=yes
  FCSTBEGIN=${FCSTBEGIN:-YES}
  NDAYS=2
  IDVC=3
  NHOURSB1=0
  NHOURSB2="/nwprod/exec/global_sighdr \$RUNDIR/sigr1i\$FM ifhr"
  MP_LABELIO=yes
  XLSMPOPTS2=parthds=\$THREADS2:spins=0:yields=0:stack=512000000
  fmax=$nhours
  FOUT=3
  FZER=3
  FCYC=0
  FDFI=0
  FHRES=$nhours
  NSOUT=${NSOUT:-0}
  GFSIO_IN=".true."
  GFSIO_OUT=".true."
  DYNVARS=liope=.F.,gfsio_in=\$gfsio_in,gfsio_out=\$gfsio_out
  PHYVARS=liope=.F.,gfsio_in=\$gfsio_in,gfsio_out=\$gfsio_out
  TRACERVARS=RI=\$RIlist,CPI=\$CPIlist,
  NGPTC=12
  NGPTC=${NGPTC:-$((wave/10))}
  LEVR=${levr:-0}
  COMOUT=$RUNDIR
  POSTGPDIR=$TOPDIR/wx23hh/00wkgfs/src/global_postgp.fd
  POSTGPEXEC=$POSTGPDIR/global_postgp
  POSTGPSH=/nwprod/ush/global_postgp.sh

 #####################################GFS RUN SUBMISSION SCRIPT OPTIONS##################################
cat ll_templates/gfs_prep_script_tmp | sed s:_JOB_NAME1_:${JOB_NAME1}:g \
             | sed s:_ERROR_:${ERROR1}:g \
	     | sed s:_OUTPUT_:${OUTPUT1}:g \
             | sed s:_JOB_TYPE_:${JOB_TYPE}:g \
             | sed s:_TASK_PER_NODE_:${TASK_PER_NODE1}:g \
             | sed s:_NODE_:${NODE1}:g \
	     | sed s:_GROUP_:${GROUP}:g \
             | sed s:_NODEUSAGE_:${NODE_USAGE}:g \
             | sed s:_CONSUMABLE_MEMORY_:${CONSUMABLE_MEMORY}:g \
             | sed s:_CLASS_:${CLASS1}:g \
             | sed s:_WALL_CLOCK_TIME_:${WALL_CLOCK_TIME1}:g \
             | sed s:_ACCT_NO_:${ACCT_NO}:g \
	     | sed s:_CDATE_:${CDATE}:g \
	     | sed s:_SFCPRESS_ID_:${SFCPRESS_ID}:g \
	     | sed s:_THERMODYN_ID_:${THERMODYN_ID}:g \
	     | sed s:_HYBRID_:${HYBRID}:g \
             | sed s:_SPECTRAL_LOOP_:${SPECTRAL_LOOP}:g \
             | sed s:_GENCOORDHYBRID_:${GEN_COORD_HYBRID}:g \
             | sed s:_ENS_NUM_:${ENS_NUM}:g \
             | sed s:_TASKS_:${TASKS}:g \
	     | sed s:_PE1_:${PE2}:g \
             | sed s:_HYBRID_:${HYBRID}:g \
	     | sed s:_WAVE_:${WAVE}:g \
	     | sed s:_LM_:${LM}:g \
	     | sed s:_LSOIL_:${LSOIL}:g \
             | sed s:_MP_SHARED_MEMORY_:${MP_SHARED_MEMORY}:g \
	     | sed s:_MP_STDOUTMODE_:${MP_STDOUTMODE}:g \
	     | sed s:_MEMORY_AFFINITY_:${MEMORY_AFFINITY}:g \
             | sed s:_THREADS_:${THREADS1}:g \
	     | sed s:_SPINLOOPTIME_:${SPINLOOPTIME}:g \
             | sed s:_YIELDLOOPTIME_:${YIELDLOOPTIME}:g \
	     | sed s:_AIXTHREAD_SCOPE_:${AIXTHREAD_SCOPE}:g \
	     | sed s:_MALLOCMULTIHEAP_:${MALLOCMULTIHEAP}:g \
	     | sed s:_EXPERIMENT_NAME_:${EXPERIMENT_NAME}:g \
             | sed s:_NSCDIR_:${NSCDIR}:g \
             | sed s:_TOPDIR_:${TOPDIR}:g \
             | sed s:_NDAYS_:${NDAYS}:g \
	     | sed s:_DUMPDIR_:${DUMPDIR}:g \
	     | sed s:_MP_COREFILE_FORMAT_:${MP_COREFILE_FORMAT}:g \
	     | sed s:_RUNDIR_:${RUNDIR}:g \
	     | sed s:_SCRIPTS_:${SCRIPTS}:g \
	     | sed s:_COLD_SFC_:${COLD_SFC}:g \
	     | sed s:_NTRC_:${NTRC}:g \
	     | sed s:_VARID_:${VARID}:g \
	     | sed s:_NUMCLD_:${NUMCLD}:g \
	     | sed s:_CHANGE_RES_:${CHANGE_RES}:g \
	     | sed s:_IC_TYPE_:${IC_TYPE}:g \
	     | sed s:_NTRAC_:${NTRAC}:g \
	     | sed s:_NTOZ_:${NTOZ}:g \
             | sed s:_NTCW_:${NTCW}:g \
	     | sed s:_NCLD_:${NCLD}:g \
	     | sed s:_NMTVR_:${NMTVR}:g \
             | sed s:_LSM_:${LSM}:g \
	     | sed s:_SIGLEVEL1_:${SIGLEVEL1}:g \
	     | sed s:_APERCENT_:${APERCENT}:g \
             | sed s:_SIGLEVEL2_:${SIGLEVEL2}:g \
             | sed s:_IDVM1_:${IDVM1}:g \
	     | sed s:_IDSL1_:${IDSL1}:g \
	     | sed s:_IDVC1_:${IDVC1}:g \
	     | sed s:_NVCOORD1_:${NVCOORD1}:g \
	     | sed s:_IVSSFC1_:${IVSSFC1}:g \
	     | sed s:_IDVC2_:${IDVC2}:g \
             | sed s:_IDSL2_:${IDSL2}:g \
             | sed s:_IVSSFC2_:${IVSSFC2}:g \
	     | sed s:_IDVM2_:${IDVM2}:g \
	     | sed s:_IVSSIG_:${IVSSIG}:g \
             | sed s:_NVCOORD2_:${NVCOORD2}:g \
	     | sed s:_LATCH_:${LATCH}:g \
	     | sed s:_IDVC3_:${IDVC3}:g \
	     | sed s:_IVSSFC3_:${IVSSFC3}:g \
             | sed s:_FIXGLOBAL_:${FIXGLOBAL}:g \
	     | sed s:_CHGRESDIR_:${CHGRESDIR}:g \
	     | sed s:_CHGRESEXEC_:${CHGRESEXEC}:g \
	     | sed s:_CHGRESSH1_:${CHGRESSH1}:g \
	     | sed s:_LANDICE_OPT_:${LANDICE_OPT}:g \
             | sed s:_CLIMO_FIELDS_OPT_:${CLIMO_FIELDS_OPT}:g \
             | sed s:_OROGRAPHY_:${OROGRAPHY}:g \
	     | sed s:_MTNVAR_:${MTNVAR}:g \
	     | sed s:_CDUMP_:${CDUMP}:g \
             | sed s:_CHGRESVARS1_:${CHGRESVARS1}:g \
             | sed s:_SIGLEVEL3_:${SIGLEVEL3}:g \
	     | sed s:_LONSPERLAT_:${LONSPERLAT}:g \
	     | sed s:_CHGRESSH2_:${CHGRESSH2}:g \
	     | sed s:_CHGRESVARS2_:${CHGRESVARS2}:g \
             > ../job/prep_gfs_script
cat ll_templates/gfs_run_script_tmp | sed s:_JOB_NAME_:${JOB_NAME1}:g \
             | sed s:_ERROR_:${ERROR2}:g \
             | sed s:_OUTPUT_:${OUTPUT2}:g \
             | sed s:_JOB_TYPE_:${JOB_TYPE}:g \
             | sed s:_TASK_PER_NODE_:${TASK_PER_NODE2}:g \
             | sed s:_IDVC_:${IDVC}:g \
             | sed s:_NODE_:${NODE2}:g \
             | sed s:_NODEUSAGE_:${NODE_USAGE}:g \
             | sed s:_CONSUMABLE_MEMORY_:${CONSUMABLE_MEMORY}:g \
             | sed s:_CLASS_:${CLASS2}:g \
	     | sed s:_SPECTRAL_LOOP_:${SPECTRAL_LOOP}:g \
             | sed s:_GROUP_:${GROUP}:g \
             | sed s:_WALL_CLOCK_TIME_:${WALL_CLOCK_TIME2}:g \
             | sed s:_ACCT_NO_:${ACCT_NO}:g \
             | sed s:_CDATE_:${CDATE}:g \
             | sed s:_ADIAB_:${ADIAB}:g \
             | sed s:_NDAYS_:${NDAYS}:g \
             | sed s:__FCSTBEGIN_:${FCSTBEGIN}:g \
             | sed s:_NHOURSB1_:${NHOURSB1}:g \
             | sed s:_NHOURSB2_:"${NHOURSB2}":g \
             | sed s:_GFSIO_IN_:${GFSIO_IN}:g \
             | sed s:_GFSIO_OUT_:${GFSIO_OUT}:g \
	     | sed s:_REDUCED_GRID_:${REDUCED_GRID}:g \
             | sed s:_SFCPRESS_ID_:${SFCPRESS_ID}:g \
             | sed s:_THERMODYN_ID_:${THERMODYN_ID}:g \
             | sed s:_HYBRID_:${HYBRID}:g \
             | sed s:_GENCOORDHYBRID_:${GEN_COORD_HYBRID}:g \
             | sed s:_ENS_NUM_:${ENS_NUM}:g \
             | sed s:_TASKS_:${TASKS2}:g \
             | sed s:_PE_:${PE2}:g \
             | sed s:_HYBRID_:${HYBRID}:g \
             | sed s:_WAVE_:${WAVE}:g \
             | sed s:_LM_:${LM}:g \
             | sed s:_FOUT_:${FOUT}:g \
             | sed s:_FZER_:${FZER}:g \
             | sed s:_FCYC_:${FCYC}:g \
             | sed s:_FDFI_:${FDFI}:g \
             | sed s:_NSOUT_:${NSOUT}:g \
             | sed s:_LSOIL_:${LSOIL}:g \
             | sed s:_MP_SHARED_MEMORY_:${MP_SHARED_MEMORY}:g \
             | sed s:_MP_STDOUTMODE_:${MP_STDOUTMODE}:g \
             | sed s:_MEMORY_AFFINITY_:${MEMORY_AFFINITY}:g \
             | sed s:_THREADS_:${THREADS2}:g \
             | sed s:_SPINLOOPTIME_:${SPINLOOPTIME}:g \
             | sed s:_YIELDLOOPTIME_:${YIELDLOOPTIME}:g \
             | sed s:_AIXTHREAD_SCOPE_:${AIXTHREAD_SCOPE}:g \
             | sed s:_MALLOCMULTIHEAP_:${MALLOCMULTIHEAP}:g \
             | sed s:_EXPERIMENT_NAME_:${EXPERIMENT_NAME}:g \
             | sed s:_NSCDIR_:${NSCDIR}:g \
             | sed s:_TOPDIR_:${TOPDIR}:g \
             | sed s:_DUMPDIR_:${DUMPDIR}:g \
             | sed s:_MP_COREFILE_FORMAT_:${MP_COREFILE_FORMAT}:g \
             | sed s:_RUNDIR_:${RUNDIR}:g \
             | sed s:_SCRIPTS_:${SCRIPTS}:g \
             | sed s:_COLD_SFC_:${COLD_SFC}:g \
             | sed s:_NTRC_:${NTRC}:g \
             | sed s:_VARID_:${VARID}:g \
             | sed s:_NUMCLD_:${NUMCLD}:g \
             | sed s:_CHANGE_RES_:${CHANGE_RES}:g \
             | sed s:_IC_TYPE_:${IC_TYPE}:g \
             | sed s:_NTRAC_:${NTRAC}:g \
             | sed s:_NTOZ_:${NTOZ}:g \
             | sed s:_NTCW_:${NTCW}:g \
             | sed s:_NCLD_:${NCLD}:g \
             | sed s:_NMTVR_:${NMTVR}:g \
             | sed s:_LSM_:${LSM}:g \
             | sed s:_SIGLEVEL1_:${SIGLEVEL1}:g \
             | sed s:_APERCENT_:${APERCENT}:g \
             | sed s:_SIGLEVEL2_:${SIGLEVEL2}:g \
             | sed s:_IDVM1_:${IDVM1}:g \
             | sed s:_IDSL1_:${IDSL1}:g \
             | sed s:_IDVC1_:${IDVC1}:g \
             | sed s:_NVCOORD1_:${NVCOORD1}:g \
             | sed s:_IVSSFC1_:${IVSSFC1}:g \
             | sed s:_IDVC2_:${IDVC2}:g \
             | sed s:_IDSL2_:${IDSL2}:g \
             | sed s:_IVSSFC2_:${IVSSFC2}:g \
             | sed s:_IDVM2_:${IDVM2}:g \
             | sed s:_IVSSIG_:${IVSSIG}:g \
             | sed s:_NVCOORD2_:${NVCOORD2}:g \
             | sed s:_LATCH_:${LATCH}:g \
             | sed s:_IDVC3_:${IDVC3}:g \
             | sed s:_IVSSFC3_:${IVSSFC3}:g \
             | sed s:_FIXGLOBAL_:${FIXGLOBAL}:g \
             | sed s:_CHGRESDIR_:${CHGRESDIR}:g \
             | sed s:_CHGRESEXEC_:${CHGRESEXEC}:g \
             | sed s:_CHGRESSH1_:${CHGRESSH1}:g \
             | sed s:_LANDICE_OPT_:${LANDICE_OPT}:g \
             | sed s:_CLIMO_FIELDS_OPT_:${CLIMO_FIELDS_OPT}:g \
             | sed s:_OROGRAPHY_:${OROGRAPHY}:g \
             | sed s:_MTNVAR_:${MTNVAR}:g \
             | sed s:_CDUMP_:${CDUMP}:g \
             | sed s:_GFSIO_IN_:${GFSIO_IN}:g \
             | sed s:_DYNVARS__:${DYNVARS}:g \
             | sed s:_PHYVARS_:${PHYVARS}:g \
             | sed s:_TRACEVARS_:${TRACEVARS}:g \
             | sed s:_GFSIO_OUT_:${GFSIO_OUT}:g \
             | sed s:_CHGRESVARS1_:${CHGRESVARS1}:g \
             | sed s:_SIGLEVEL3_:${SIGLEVEL3}:g \
             | sed s:_LONSPERLAT_:${LONSPERLAT}:g \
             | sed s:_CHGRESSH2_:${CHGRESSH2}:g \
             | sed s:_CHGRESVARS2_:${CHGRESVARS2}:g \
             | sed s:_FCSTSCRIPT_:${FCSTSCRIPT}:g \
             | sed s:_FCSTEXEC_:${FCSTEXEC}:g \
             > ../job/run_gfs_script
fi
exit
