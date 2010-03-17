#!/usr/bin/ksh

set -e

OS=`uname`
echo "platform is "${OS}
settings_file=../compile_settings
if [[ -a "$settings_file" ]]
then
echo "modifying compile settings file"
rm $settings_file
fi
# reading content of file into arrays
DYN=0
PHY=0
CHEM=0
read
read
INDEX=0
while [[ $INDEX -le 2 ]];do
   read NAMES NUMBERS
   if [[ $NUMBERS -eq 1 ]] then
   NEMS_DYN[$DYN]=$NAMES
   ((DYN=$DYN+1))
   fi  
  ((INDEX=$INDEX+1))
done 
((DYN=$DYN-1))
read
read
read
INDEX=0
while [[ $INDEX -le 1 ]];do
   read NAMES NUMBERS
   if [[ $NUMBERS -eq 1 ]] then
   NEMS_PHY[$PHY]=$NAMES
   ((PHY=$PHY+1))
   fi
  ((INDEX=$INDEX+1))
done
((PHY=$PHY-1))
read
read
read
INDEX=0
while [[ $INDEX -le 0 ]];do
   read NAMES NUMBERS
   if [[ $NUMBERS -eq 1 ]] then
   NEMS_CHEM[$CHEM]=$NAMES
   ((CHEM=$CHEM+1))
   fi
  ((INDEX=$INDEX+1))
done
((CHEM=$CHEM-1))

COMPILE_TEMPLATE_DYN_NMM=0
COMPILE_TEMPLATE_DYN_GFS=0
COMPILE_TEMPLATE_DYN_FIM=0
COMPILE_TEMPLATE_PHY_NMM=0
COMPILE_TEMPLATE_PHY_GFS=0


COMPILE_TEMPLATE_DYN=0
COMPILE_TEMPLATE_PHY=0
COMPILE_TEMPLATE_CHEM=0

if [[ ${#NEMS_DYN[*]} -eq 0 ]] then
   echo "At least one dynamical and one physics package must be selected."
   echo "Exiting."
   exit
fi
if [[ ${#NEMS_PHY[*]} -eq 0 ]] then
   echo "At least one dynamical and one physics package must be selected."
   echo "Exiting."
   exit
fi
INDEX=0
while [[ $INDEX -le $DYN ]];do
if [[ ${NEMS_DYN[$INDEX]} == 'fim_dyn' ]] then
   echo "fim dynamical core not yet integrated into NEMS."
   echo "Exiting."
fi
((INDEX=INDEX+1))
done
INDEX=0
while [[ $INDEX -le $DYN ]];do
if [[ ${NEMS_DYN[$INDEX]} == 'nmm_dyn' ]] then
   echo "nmm dynamics package to be compiled"
   COMPILE_TEMPLATE_DYN_NMM=1
fi
((INDEX=INDEX+1))
done
INDEX=0
while [[ $INDEX -le $DYN ]];do
if [[ ${NEMS_DYN[$INDEX]} == 'gfs_dyn' ]] then
   echo "gfs dynamics package to be compiled"
   COMPILE_TEMPLATE_DYN_GFS=2
fi
((INDEX=INDEX+1))
done
((COMPILE_TEMPLATE_DYN=$COMPILE_TEMPLATE_DYN_NMM+$COMPILE_TEMPLATE_DYN_GFS))
INDEX=0
while [[ $INDEX -le $PHY ]];do
if [[ ${NEMS_PHY[$INDEX]} == 'nmm_phys' ]] then
   echo "nmm physics package to be compiled"
   COMPILE_TEMPLATE_PHY_NMM=1
fi
((INDEX=INDEX+1))
done
INDEX=0
while [[ $INDEX -le $PHY ]];do
if [[ ${NEMS_PHY[$INDEX]} == 'gfs_phys' ]] then
   echo "gfs physics package to be compiled"
   COMPILE_TEMPLATE_PHY_GFS=2
fi
((INDEX=INDEX+1))
done
((COMPILE_TEMPLATE_PHY=$COMPILE_TEMPLATE_PHY_NMM+$COMPILE_TEMPLATE_PHY_GFS))
if [[ $COMPILE_TEMPLATE_DYN -eq 2 && $COMPILE_TEMPLATE_PHY -eq 3 ]] ; then
((COMPILE_TEMPLATE_DYN=3))
fi
if [[ ${NEMS_CHEM} == 'gocart' ]] then
   echo "gocart chemistry package to be compiled"
   COMPILE_TEMPLATE_CHEM=1
fi

NEMSTOPDIR=`cd ..; pwd`

echo $COMPILE_TEMPLATE_DYN $COMPILE_TEMPLATE_PHY > build_setup/configure_flags

. build_setup/build_settings.sh $COMPILE_TEMPLATE_DYN $COMPILE_TEMPLATE_PHY $COMPILE_TEMPLATE_CHEM $OS

. build_setup/build_makefiles.sh $COMPILE_TEMPLATE_DYN $COMPILE_TEMPLATE_PHY 

