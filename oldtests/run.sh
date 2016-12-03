#!/bin/ksh

mkdir -p ${RUNDIR}

source ./atparse.ksh

if [ ${nems_configure}"x" != "x" ]; then
  cat nems.configure.${nems_configure}.IN | atparse > nems.configure
  cp nems.configure ${RUNDIR}
fi

cp ${PATHTR}/exe/NEMS.x ${RUNDIR}

cd ${RUNDIR}

INI_YEAR=$(echo $CDATE | cut -c1-4)
INI_MONTH=$(echo $CDATE | cut -c5-6)
INI_DAY=$(echo $CDATE | cut -c7-8)
INI_HOUR=$(echo $CDATE | cut -c9-10)

NDAYS=${NDAYS:-0}
nhours=`expr $NDAYS \* 24`
FHMAX=${NHRS:-$nhours}

cat << EOF > model_configure
print_esmf:             .true.
total_member:           1
PE_MEMBER01:            $TASKS
ENS_SPS:                .false.
RUN_CONTINUE:           .false.
start_year:             $INI_YEAR
start_month:            $INI_MONTH
start_day:              $INI_DAY
start_hour:             $INI_HOUR
start_minute:           0
start_second:           0
nhours_fcst:            $FHMAX
EOF

echo "Execute NEMS.x ..."

$MPIEXEC -np $TASKS ./NEMS.x > out 2>&1

echo "...DONE!"
