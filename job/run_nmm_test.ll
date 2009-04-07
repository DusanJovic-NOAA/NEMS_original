#!/bin/ksh
#
#@ output =out
#@ error = err
#@ job_type = parallel
##@ node = 2
#@ total_tasks = 64
#@ blocking=unlimited
#@ task_affinity = core(1)
#@ node_usage=shared
#@ resources=ConsumableCPUs(1)ConsumableMemory(3000MB)
#@ class = dev
#@ wall_clock_limit = 00:30:00
#@ preferences = Feature == "dev"
#@ network.MPI = csss,shared,us
#@ account_no=NAM-T2O
#@ queue
#
#
set -x
#
### For MPI
#
export BIND_TASKS=yes
export MP_EAGER_LIMIT=65536
export MP_SHARED_MEMORY=yes
export MP_SINGLE_THREAD=yes
export MP_LABELIO=yes
export MP_STDOUTMODE=ordered
export MEMORY_AFFINITY=mcm
export MP_COREFILE_FORMAT=core.txt
export TARGET_CPU_RANGE="-1"
export XLSMPOPTS="parthds=1:stack=128000000"
export AIXTHREAD_SCOPE=S


RUNDIR=/ptmp/${USER}/trunk ### change here
DATADIR=/meso/save/${USER}/tomscase
SRCDIR=/meso/save/${USER}/trunk   ##### change here

mkdir -p $RUNDIR
cd $RUNDIR
ls

echo $DATADIR
ls $DATADIR

cp $DATADIR/input_nmmb_regional.d01 main_input_filename
cp $DATADIR/boco* . 
cp $SRCDIR/exp/configure_files/configfile_regional_new $RUNDIR/configure_file

# cp $DATADIR/tr* .

cp /nwprod/fix/nam_micro_lookup.dat ETAMPNEW_DATA
# cp /meso/save/wx22tb/WRF_repository/run/RRT* .
# cp /meso/noscrub/wx20rv/NOAH_wrf_data/*.TBL .

cp /meso/save/wx20py/WRFV2_ijk/run/RRT* .
cp /meso/save/wx20py/WRFV2_ijk/run/*.TBL .
cp /meso/save/wx20py/WRFV2_ijk/run/tr* .


echo "Model started:  " `date`

poe  $SRCDIR/exe/NMM_NEMS.x

echo "Model ended:    " `date`

exit
