#!/bin/ksh

if [ ! -a configure_settings ] ; then
   echo Cannot find file configure_settings needed to build NEMS ... exiting
   exit
fi

if [ ! -a ../exe/makedepf90 ] ; then
cd ../util/makedepf90-2.8.8
./build_makedep
mv makedepf90 ../../exe
cd ../../build
fi


. build_setup/build_driver.sh < configure_settings
cd build_setup
var1=`awk -F" " 'NR==1{print $1}' configure_flags`
var2=`awk -F" " 'NR==1{print $2}' configure_flags`
rm configure_flags
cd ..

if [[ $var1 -eq 3 && $var2 -eq 3 ]] ; then
./compile.sh
else
./compile_stub.sh
./compile.sh
fi

echo "done building NEMS."
exit
