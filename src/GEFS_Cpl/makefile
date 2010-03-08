##  This is an elementary grid component or say independent grid component ##
include ../../ush/compile_settings
MAKEFILE=makefile
TARGET = GEFS_Cpl.a

 ARFLAGS = -r

.SUFFIXES: .o .f .F .h
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
	GEFS_Sto_Per_Scheme_Step1.o         \
	kinds.o                             \
	peuc.o                              \
	pran.o                              \
	prana.o                             \
	Cal_Sto_Coef.o

#
#
# *****************************************************************
#
default:
	@gmake -f $(MAKEFILE) all

all: $(LIBDIR)/lib/$(TARGET)

$(LIBDIR)/lib/$(TARGET):  $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	echo " $(TARGET) created!"

clean:
	rm -f *.o *.mod *.a

.F.o:
	$(F77) $(FFLAGS) -c -d $< 
	#$(F77) $(FFLAGS) -c -d -WF,-DCLR:${RASV} $< 
.f.o:
	$(F77) $(FFLAGS) -c $< 

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

kinds.o: GEFS_Cpl_Cal_Sto_Coef.fd/kinds.f90
	$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/kinds.f90

peuc.o: GEFS_Cpl_Cal_Sto_Coef.fd/peuc.f90
	$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/peuc.f90

pran.o: GEFS_Cpl_Cal_Sto_Coef.fd/pran.f90
	$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/pran.f90

prana.o: GEFS_Cpl_Cal_Sto_Coef.fd/prana.f90
	$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/prana.f90

Cal_Sto_Coef.o: GEFS_Cpl_Cal_Sto_Coef.fd/Cal_Sto_Coef.f90
	$(F90) $(FFLAG90) -c GEFS_Cpl_Cal_Sto_Coef.fd/Cal_Sto_Coef.f90
