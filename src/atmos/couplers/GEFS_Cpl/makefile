##  This is an elementary grid component or say independent grid component ##
include ../../../../build/compile_settings
MAKEFILE=makefile
TARGET = GEFS_Cpl.a

#
# *****************************************************************
#
OBJS =  Lib_ESMFStateAddGetMod.o            \
	machine.o                           \
	GEFS_Cpl_InternalState_ESMFMod.o    \
	GEFS_bcst_global.o                  \
	GEFS_Sto_Per_Scheme_Step2.o         \
	DistributeForStep1.o                \
	GEFS_CplState_ESMFMod.o             \
	GEFS_Cpl_Run_ESMFMod.o              \
	GEFS_Cpl_ESMFMod.o                  \
	GEFS_GetParameterFromStateMod.o     \
	GEFS_CplComp_ESMFMod.o              \
	GEFS_Sto_Per_Scheme_Step1.o         

OBJS_90=kinds.o                             \
	peuc.o                              \
	pran.o                              \
	prana.o                             \
	Cal_Sto_Coef.o

#
#
# *****************************************************************
#

OBJS_f = ${OBJS:.o=.f}
OBJS_90_f90 = ${OBJS_90:.o=.f90}

.SUFFIXES:    .o .f .f90

default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)

$(LIBDIR)/lib/$(TARGET):  $(OBJS) $(OBJS_90)
	$(AR) $(ARFLAGS) $@ $?

.f.o: %.o: %.f
	$(F77) $(FFLAGS) -c $*.f

.f90.o: %.o: %.f90
	$(F90) $(FFLAGS) -c $*.f90

depend: $(OBJS_f) $(OBJS_90_f90)
	$(MAKEDEPEXEC)

clean:
	rm -f  $(LIBDIR)/lib/$(TARGET)
	rm -f *.o 
	rm -f *.mod 
	rm -f depend


Lib_ESMFStateAddGetMod.o: Lib_ESMFStateAddGetMod.f
	$(F77) $(FFLAG90) -c Lib_ESMFStateAddGetMod.f

machine.o: machine.f
	$(F77) $(FFLAGB) -c machine.f

GEFS_Cpl_InternalState_ESMFMod.o: GEFS_Cpl_InternalState_ESMFMod.f
	$(F77) $(FFLAG90) -c GEFS_Cpl_InternalState_ESMFMod.f

DistributeForStep1.o: DistributeForStep1.f
	$(F77) $(FFLAG90) -c DistributeForStep1.f

GEFS_CplState_ESMFMod.o: GEFS_CplState_ESMFMod.f
	$(F77) $(FFLAG90) -c GEFS_CplState_ESMFMod.f

GEFS_Sto_Per_Scheme_Step1.o: GEFS_Sto_Per_Scheme_Step1.f
	$(F77) $(FFLAG90) -c GEFS_Sto_Per_Scheme_Step1.f

GEFS_Sto_Per_Scheme_Step2.o: GEFS_Sto_Per_Scheme_Step2.f
	$(F77) $(FFLAG90) -c GEFS_Sto_Per_Scheme_Step2.f

GEFS_bcst_global.o: GEFS_bcst_global.f
	$(F77) $(FFLAG90) -c GEFS_bcst_global.f

GEFS_Cpl_Run_ESMFMod.o:	GEFS_Cpl_Run_ESMFMod.f
	$(F77) $(FFLAG90) -c GEFS_Cpl_Run_ESMFMod.f

GEFS_Cpl_ESMFMod.o: GEFS_Cpl_ESMFMod.f
	$(F77) $(FFLAG90) -c GEFS_Cpl_ESMFMod.f

GEFS_GetParameterFromStateMod.o: GEFS_GetParameterFromStateMod.f
	$(F77) $(FFLAG90) -c GEFS_GetParameterFromStateMod.f

GEFS_CplComp_ESMFMod.o:	GEFS_CplComp_ESMFMod.f
	$(F77) $(FFLAG90) -c GEFS_CplComp_ESMFMod.f

kinds.o: kinds.f90
	$(F90) $(FFLAG90) -c kinds.f90

peuc.o: peuc.f90
	$(F90) $(FFLAG90) -c peuc.f90

pran.o: pran.f90
	$(F90) $(FFLAG90) -c pran.f90

prana.o: prana.f90
	$(F90) $(FFLAG90) -c prana.f90

Cal_Sto_Coef.o: Cal_Sto_Coef.f90
	$(F90) $(FFLAG90) -c Cal_Sto_Coef.f90
