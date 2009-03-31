#!/bin/ksh
set -x 
#
# this script makes the right files before compilation
#
#  ---  select lw radiation package
#QLW=rrtm0
 QLW=rrtm1
#QLW=rrtm2
#QLW=gfdl0
#  ---  select sw radiation package
#QSW=ncep0
#QSW=ncep1
#QSW=rrtm1
 QSW=rrtm2
#
 export COMPILE_RAD=YES
#
#CP1=/u/wx20mi/bin/ncp
CP1=cp
#
export COMPILE_RAD=${COMPILE_RAD:-NO}
if [ $COMPILE_RAD = YES ] ; then
 $CP1 radlw_${QLW}_datatb.f  radlw_datatb.f
 $CP1 radlw_${QLW}_param.f   radlw_param.f
 $CP1 radlw_${QLW}_main.f    radlw_main.f
 $CP1 radsw_${QSW}_datatb.f  radsw_datatb.f
 $CP1 radsw_${QSW}_param.f   radsw_param.f
 $CP1 radsw_${QSW}_main.f    radsw_main.f
fi
#
