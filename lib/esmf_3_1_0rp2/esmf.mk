# ESMF application makefile fragment
#
# Use the following ESMF_ variables to compile and link
# your ESMF application against this ESMF build.
#
# !!! VERY IMPORTANT: If the location of this ESMF build is   !!!
# !!! changed, e.g. libesmf.a is copied to another directory, !!!
# !!! this file - esmf.mk - must be edited to adjust to the   !!!
# !!! correct new path                                        !!!
#
# Please see end of file for options used on this ESMF build
#

ESMF_F90COMPILER=mpxlf90_r
ESMF_F90LINKER=mpxlf90_r

ESMF_F90COMPILEOPTS=-O -q64
ESMF_F90COMPILEPATHS=-I/meso/save/wx20rv/ESMF_libs/esmf_3_1_0rp2/mod/modO/AIX.default.64.mpi.default -I/meso/save/wx20rv/ESMF_libs/esmf_3_1_0rp2/src/include
ESMF_F90COMPILECPPFLAGS= -WF,-DS64=1 -WF,-DESMF_OS_AIX=1

ESMF_F90LINKOPTS= -q64
ESMF_F90LINKPATHS=-L/meso/save/wx20rv/ESMF_libs/esmf_3_1_0rp2/lib/libO/AIX.default.64.mpi.default
ESMF_F90LINKRPATHS=
ESMF_F90LINKLIBS= -lmpi_r -lxlf90_r -lC_r
ESMF_F90ESMFLINKLIBS=-lesmf  -lmpi_r -lxlf90_r -lC_r

ESMF_CXXCOMPILER=mpCC_r
ESMF_CXXLINKER=mpCC_r

ESMF_CXXCOMPILEOPTS=-O -DNDEBUG -q64
ESMF_CXXCOMPILEPATHS=-I/meso/save/wx20rv/ESMF_libs/esmf_3_1_0rp2/src/include
ESMF_CXXCOMPILECPPFLAGS=-DS64=1 -DESMF_OS_AIX=1 -D__SDIR__=''

ESMF_CXXLINKOPTS= -q64
ESMF_CXXLINKPATHS=-L/meso/save/wx20rv/ESMF_libs/esmf_3_1_0rp2/lib/libO/AIX.default.64.mpi.default
ESMF_CXXLINKRPATHS=
ESMF_CXXLINKLIBS= -lmpi_r -lm_r -lxlf90_r -lC_r
ESMF_CXXESMFLINKLIBS=-lesmf  -lmpi_r -lm_r -lxlf90_r -lC_r

#
# !!! The following options were used on this ESMF build !!!
#
# ESMF_DIR: /meso/save/wx20rv/ESMF_libs/esmf_3_1_0rp2
# ESMF_OS: AIX
# ESMF_MACHINE: default
# ESMF_ABI: 64
# ESMF_COMPILER: default
# ESMF_BOPT: O
# ESMF_COMM: mpi
# ESMF_SITE: default
# ESMF_PTHREADS: ON
# ESMF_ARRAY_LITE: FALSE
# ESMF_NO_INTEGER_1_BYTE: FALSE
# ESMF_NO_INTEGER_2_BYTE: FALSE
# ESMF_FORTRANSYMBOLS: default
# ESMF_TESTEXHAUSTIVE: OFF
# ESMF_TESTWITHTHREADS: OFF
# ESMF_TESTMPMD: OFF
# 
# ESMF environment variables pointing to 3rd party software:
