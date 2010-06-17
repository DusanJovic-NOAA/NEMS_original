include ../conf/configure.nems

MAKEFILE = makefile

LIBRARY  = ENS_Cpl.a

OPTS     = $(FFLAGS_GFS) $(R8) -O3 $(UTILINCS) $(GOCARTINCS)

FFLAG90  = $(OPTS) $(FREE)
FFLAGS1  = $(OPTS) $(FIXED)

OBJS =  Lib_ESMFStateAddGetMod.o            \
	machine.o                           \
	ENS_Cpl_InternalState_ESMFMod.o     \
	ENS_bcst_global.o                   \
	ENS_Sto_Per_Scheme_Step2.o          \
	DistributeForStep1.o                \
	ENS_CplState_ESMFMod.o              \
	ENS_Cpl_Run_ESMFMod.o               \
	ENS_Cpl_ESMFMod.o                   \
	ENS_GetParameterFromStateMod.o      \
	ENS_CplComp_ESMFMod.o               \
	ENS_Sto_Per_Scheme_Step1.o         

OBJS_90=kinds.o                             \
	peuc.o                              \
	pran.o                              \
	prana.o                             \
	Cal_Sto_Coef.o

MODULES_STUB = ENS_CplComp_ESMFMod_STUB.o

DEPEND_FILES = ${OBJS:.o=.f} ${OBJS_90:.o=.f90}


all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY):  $(OBJS) $(OBJS_90)
	$(AR) $(ARFLAGS) $@ $?

stub: $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(TARGET) $(MODULES_STUB)

.f.o: %.o: %.f
	$(FC) $(FFLAGS1) -c $*.f

.f90.o: %.o: %.f90
	$(FC) $(FFLAGS1) -c $*.f90

$(MODULES_STUB): %.o: %.f
	$(FC) $(FFLAGS) -c $*.f


clean:
	$(RM) -f $(LIBRARY) *.o *.mod *.lst lm map depend

#
# *****************************************************************
#


Lib_ESMFStateAddGetMod.o: Lib_ESMFStateAddGetMod.f
	$(FC) $(FFLAG90) -c Lib_ESMFStateAddGetMod.f

machine.o: machine.f
	$(FC) $(FFLAGS1) -c machine.f

ENS_Cpl_InternalState_ESMFMod.o: ENS_Cpl_InternalState_ESMFMod.f
	$(FC) $(FFLAG90) -c ENS_Cpl_InternalState_ESMFMod.f

DistributeForStep1.o: DistributeForStep1.f
	$(FC) $(FFLAG90) -c DistributeForStep1.f

ENS_CplState_ESMFMod.o: ENS_CplState_ESMFMod.f
	$(FC) $(FFLAG90) -c ENS_CplState_ESMFMod.f

ENS_Sto_Per_Scheme_Step1.o: ENS_Sto_Per_Scheme_Step1.f
	$(FC) $(FFLAG90) -c ENS_Sto_Per_Scheme_Step1.f

ENS_Sto_Per_Scheme_Step2.o: ENS_Sto_Per_Scheme_Step2.f
	$(FC) $(FFLAG90) -c ENS_Sto_Per_Scheme_Step2.f

ENS_bcst_global.o: ENS_bcst_global.f
	$(FC) $(FFLAG90) -c ENS_bcst_global.f

ENS_Cpl_Run_ESMFMod.o:	ENS_Cpl_Run_ESMFMod.f
	$(FC) $(FFLAG90) -c ENS_Cpl_Run_ESMFMod.f

ENS_Cpl_ESMFMod.o: ENS_Cpl_ESMFMod.f
	$(FC) $(FFLAG90) -c ENS_Cpl_ESMFMod.f

ENS_GetParameterFromStateMod.o: ENS_GetParameterFromStateMod.f
	$(FC) $(FFLAG90) -c ENS_GetParameterFromStateMod.f

ENS_CplComp_ESMFMod.o:	ENS_CplComp_ESMFMod.f
	$(FC) $(FFLAG90) -c ENS_CplComp_ESMFMod.f

kinds.o: kinds.f90
	$(FC) $(FFLAG90) -c kinds.f90

peuc.o: peuc.f90
	$(FC) $(FFLAG90) -c peuc.f90

pran.o: pran.f90
	$(FC) $(FFLAG90) -c pran.f90

prana.o: prana.f90
	$(FC) $(FFLAG90) -c prana.f90

Cal_Sto_Coef.o: Cal_Sto_Coef.f90
	$(FC) $(FFLAG90) -c Cal_Sto_Coef.f90



#
#
#

include ../conf/make.rules

include depend
