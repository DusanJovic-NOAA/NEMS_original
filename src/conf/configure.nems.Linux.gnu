## NEMS configuration file
##
## Platform: generic Linux
## Compiler: GNU compilers through wrappers determined by ESMF

SHELL           = /bin/ksh

################################################################################
## Include the common configuration parts

include         $(TOP)/conf/configure.nems.NUOPC

################################################################################
## Settings to be fed into the NEMS makefile. Only the following variables are
## used by NEMS makefile:
##  * EXTLIBS       ( for non-post targets)
##  * EXTLIBS_POST  ( for post targets)
##  * FFLAGS
##  * FFLAGS_NMM
##  * FFLAGS_GFS
##  * FFLAGS_GFSF
##  * FFLAGS_GEN
##  * FFLAGS_FIM
##  * FC
##  * CC
##  * CPP
##  * CPPFLAGS
##  * RM
################################################################################

NETCDF_LIB  = 

LIBDIR      = /twixhome/gerhard/WorkESMF/NUOPC-NEMS/TWIX/nemslibs

NEMSIO_INC  = -I${LIBDIR}/incmod/nemsio
NEMSIO_LIB  = -L${LIBDIR} -lnemsio

W3_LIB      = -L${LIBDIR} -lw3nco_d -lw3emc_d
BACIO_LIB   = -L${LIBDIR} -lbacio_4
SP_LIB      = -L${LIBDIR} -lsp_d
SYS_LIB     =

EXTLIBS     = $(NEMSIO_LIB) \
              $(W3_LIB) \
              $(BACIO_LIB) \
              $(SP_LIB) \
              $(ESMF_LIB) \
              $(NETCDF_LIB) \
              $(SYS_LIB)

## for the post quilting option
POSTDIR     = ${LIBDIR}
POST_LIB    = -L${POSTDIR} -lnceppost
W3_POST_LIB = -L${LIBDIR} -lw3nco_4 -lw3emc_4
CRTM_LIB    = -L${LIBDIR} -lcrtm_v2.0.6
G2_LIB      = -L${LIBDIR} -lg2tmpl -lg2_4 -ljasper -lpng -lz
XML_LIB     = -L${LIBDIR} -lxmlparse
SIGIO_LIB   = -L${LIBDIR} -lsigio_4
SFCIO_LIB   = -L${LIBDIR} -lsfcio_4

EXTLIBS_POST = $(POST_LIB) \
               $(NEMSIO_LIB) \
               $(W3_POST_LIB) \
               $(XML_LIB) \
               $(G2_LIB) \
               $(BACIO_LIB) \
               $(SIGIO_LIB) \
               $(SFCIO_LIB) \
               $(SP_LIB) \
               $(CRTM_LIB) \
               $(ESMF_LIB) \
               $(NETCDF_LIB) \
               $(SYS_LIB)
###
FC          = $(ESMF_F90COMPILER) -fopenmp
FPP         = -cpp
CC          = $(ESMF_CXXCOMPILER) -fopenmp
FREE        = -ffree-form
FIXED       = -ffixed-form
R8          = -fdefault-real-8 -fdefault-double-8

FINCS       = $(ESMF_INC) $(NEMSIO_INC)
TRAPS       =
#TRAPS       = -ftrapuv -check all -fp-stack-check
#TRAPS       = -ftrapuv -check bounds -check format -check output_conversion -check pointers -check uninit -fp-stack-check

FFLAGS      = $(TRAPS) $(FINCS) -fno-range-check  -ffree-line-length-none

OPTS_NMM    = -O3
OPTS_GFS    = -O3
OPTS_GEN    = -O3
OPTS_FIM    = -O3

FFLAGS_NMM  = $(OPTS_NMM) $(FFLAGS)
FFLAGS_GFS  = $(OPTS_GFS) $(FFLAGS) $(FREE)
FFLAGS_GFSF = $(OPTS_GFS) $(FFLAGS) $(FIXED)
FFLAGS_GEN  = $(OPTS_GEN) $(FFLAGS)
FFLAGS_FIM  = $(OPTS_FIM) $(FFLAGS)

CPP         = /lib/cpp -P -traditional
CPPFLAGS    = -DCHNK_RRTM=8 -DENABLE_SMP -DTHREAD_2D

AR          = ar
ARFLAGS     = -r

RM          = rm
