#!/bin/ksh

mkdir -p ${RUNDIR}

if [ ${nems_configure}"x" != "x" ]; then
  cat nems.configure.${nems_configure}.IN                                   \
    | sed s:_atm_model_:${atm_model}:g                                      \
    | sed s:_atm_petlist_bounds_:"${atm_petlist_bounds}":g                  \
    | sed s:_lnd_model_:${lnd_model}:g                                      \
    | sed s:_lnd_petlist_bounds_:"${lnd_petlist_bounds}":g                  \
    | sed s:_ice_model_:${ice_model}:g                                      \
    | sed s:_ice_petlist_bounds_:"${ice_petlist_bounds}":g                  \
    | sed s:_ocn_model_:${ocn_model}:g                                      \
    | sed s:_ocn_petlist_bounds_:"${ocn_petlist_bounds}":g                  \
    | sed s:_wav_model_:${wav_model}:g                                      \
    | sed s:_wav_petlist_bounds_:"${wav_petlist_bounds}":g                  \
    | sed s:_ipm_model_:${ipm_model}:g                                      \
    | sed s:_ipm_petlist_bounds_:"${ipm_petlist_bounds}":g                  \
    | sed s:_hyd_model_:${hyd_model}:g                                      \
    | sed s:_hyd_petlist_bounds_:"${hyd_petlist_bounds}":g                  \
    | sed s:_med_model_:${med_model}:g                                      \
    | sed s:_med_petlist_bounds_:"${med_petlist_bounds}":g                  \
    | sed s:_atm_coupling_interval_sec_:"${atm_coupling_interval_sec}":g    \
    | sed s:_ocn_coupling_interval_sec_:"${ocn_coupling_interval_sec}":g    \
    | sed s:_coupling_interval_sec_:"${coupling_interval_sec}":g            \
    | sed s:_coupling_interval_slow_sec_:"${coupling_interval_slow_sec}":g  \
    | sed s:_coupling_interval_fast_sec_:"${coupling_interval_fast_sec}":g  \
    >  nems.configure
                         
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
