################################################################################
################################################################################
# 
#     Makefile for Global NMM Model 
#
#     Use:
#     make         -  build the executable
#     make clean   -  start with a clean slate
#
#     The following macros will be of interest:
#
#         TARGET   - name of the executable
#         FC       - name of Fortran compiler
#         CPP      - name of CPP
#         ARCH     - architecture
#         CPPFLAGS - CPP flags
#         OPTS     - compiler code optimizations
#         LIST     - source listing
#         SMP      - threading
#         TRAPS    - runtime traps for floating point exceptions
#         PROFILE  - source code profiling ( -pg )
#         DEBUG    - -g
#         MEM      - user data area and stack size
#         MAP      - load map
#         W3LIB    - w3lib
#         ESSL     - ESSL library
#         MASS     - MASS library
#         TRACE    - MPI tracing
#
#################################################################################

# Define the name of the executable
#
#
# CPP, Compiler, and Linker Options
#
EXEC     = atmos.a
NMM_UTIL = ../nmm_util
NMM_PHY = ../nmm_phys
NMM_DYN = ../nmm_dyn
GFS_DYN = ../gfs_dyn
GFS_PHY = ../gfs_phys
GFS_UTIL = ../gfs_util
NMM_UTIL_LIB = $(NMM_UTIL)/nmm_util.a
NMM_PHY_LIB = $(NMM_PHY)/nmm_phys.a
NMM_DYN_LIB = $(NMM_DYN)/nmm_dyn.a
GFS_DYN_LIB = $(GFS_DYN)/gfs_dyn.a
GFS_PHY_LIB = $(GFS_PHY)/gfs_phys.a
GFS_UTIL_LIB = $(GFS_UTIL)/gfs_util.a
GCINCS  = -I$(NMM_UTIL) -I$(NMM_PHY) -I$(NMM_DYN) -I$(GFS_PHY) -I$(GFS_DYN) -I$(GFS_UTIL)
GCLIBS = $(NMM_UTIL_LIB) $(NMM_PHY_LIB) $(NMM_DYN_LIB) $(GFS_PHY_LIB) $(GFS_DYN_LIB) $(GFS_UTIL_LIB)
FC       = mpxlf95_r
PROF     = -g -pg
PROF     =
ARCH     = AIX
CPP      = /lib/cpp -P
CPPFLAGS = -DENABLE_SMP -DIBM
#CPPFLAGS = -DENABLE_SMP
SMP      = -qsmp=omp
MAF      = -qfloat=nomaf
OPTS     = -O0 -qmaxmem=-1 -qarch=auto -qtune=auto $(PROF) -NS2000
#OPTS     = -O3 -qnostrict -qmaxmem=-1 -qarch=auto -qtune=auto $(PROF)
########     O3 failed in producing bit-identical results with diff. decomposition
TRAPS    = -NS2000 -qinitauto=FF911299 -qflttrap=ov:zero:inv:en -qsigtrap -C 
TRACE    = -L /usr/lib -lmpitrace
PROFILE  = 
DEBUG    =
MEM      =
MAP      = -bmap:map -bloadmap:lm
MASS     = -L/usrx/local/mass -lmass -lmassvp5
ESSL     = -L/usr/lpp/essl.rte.rs1/lib -lessl_r
#
# Assemble Options
#
FINCO    = -I/meso/save/wx20rv/ESMF_libs/esmf_test_310r/mod/modO/AIX.default.64.mpi.default
FINCS	 = $(FINCO) $(GCINCS)
FFLAGS   = $(OPTS) $(FINCS) $(LIST) $(PROFILE) $(DEBUG)
FFLAGSDB   = $(OPTS) $(FINCS) $(LIST) $(PROFILE) $(DEBUG) $(TRAPS)
FFLAGSNT = $(OPTS) $(FINCS) $(LIST) $(PROFILE) $(DEBUG)
LDFLAGS  = $(MEM) $(MAP) $(SMP) $(PROF) $(TRACE)
ESMFLIB  = -L/meso/save/wx20rv/ESMF_libs/esmf_test_310r/lib/libO/AIX.default.64.mpi.default                  \
             /meso/save/wx20rv/ESMF_libs/esmf_test_310r/lib/libO/AIX.default.64.mpi.default/libesmf.a
LIBS     = -lC ${ESSL} ${MASS} ${ESMFLIB} $(GCLIBS)
#

MODULE_AGC=     module_ATM_GRID_COMP.o

INCLUDES=	kind.inc


DEPS= $(INCLUDES)

AR = ar
ARFLAGS = -r


.SUFFIXES:	.F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

$(EXEC): module_ATM_GRID_COMP.o
	$(AR) $(ARFLAGS) $(EXEC) module_ATM_GRID_COMP.o

module_ATM_GRID_COMP.o: $(DEPS)
	$(FC) $(FFLAGSDB) -c module_ATM_GRID_COMP.F90

clean:	
	/bin/rm -f  $(TARGET) *.f90 *.lst *.o *.mod *.a lm map

