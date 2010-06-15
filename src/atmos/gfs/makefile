include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCS = -I../share -Ilibutil -Iio

GOCARTINCS = -I../../chem/gocart/AIX/include/Chem_Base \
             -I../../chem/gocart/AIX/include/GMAO_mpeu \
             -I../../chem/gocart/AIX/include/GEOSchem_GridComp


 DIR_LIBUTIL = libutil
 DIR_IO = io
 DIR_DYN = dyn
 DIR_PHY = phys
 GCINCS  = -I$(DIR_DYN) -I$(DIR_PHY)

 LIBRARY=libgfs.a

 DEBUG    = -qcheck -qsource -qlist

#--------------
# For Non-Debug.
#---------------
FFLAGM_DEBUG  = -O3

 OPTS     = -qsuffix=cpp=f -qrealsize=8 -qstrict -qxlf77=leadzero -qmaxmem=-1 -qnolm -qsmp=noauto  -qnosave $(FFLAGM_DEBUG)

 FFLAG90  = $(FFLAGS) $(OPTS)  -qfree -NS2048 $(UTILINCS) $(GCINCS) $(GOCARTINCS)

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

OBJS=module_GFS_CORE_SETUP.o module_GFS_WRITE.o module_GOCART_ROUTINES.o \
	module_GFS_INTERNAL_STATE.o module_GFS_INTEGRATE.o module_GFS_GRID_COMP.o

MODULES_STUB  = module_GFS_GRID_COMP_stub.o
#
#
# *****************************************************************
#
default:
	cd $(DIR_LIBUTIL) && gmake && cd ..
	cd $(DIR_IO) && gmake && cd ..
	cd $(DIR_DYN) && gmake && cd ..
	cd $(DIR_PHY) && gmake && cd ..
	@gmake -f makefile $(LIBRARY)

$(LIBRARY): $(OBJS) 
	$(AR) $(ARFLAGS) $@ $?

stub: $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_STUB)

clean:
	cd $(DIR_LIBUTIL) && gmake clean && cd ..
	cd $(DIR_IO) && gmake clean && cd ..
	cd $(DIR_DYN) && gmake clean && cd ..
	cd $(DIR_PHY) && gmake clean && cd ..
	rm -f *.f90 *.o *.mod *.lst *.a

$(OBJS): %.o: %.f90
	$(FC) $(FFLAG90) -c $*.f90

$(MODULES_STUB): %.o: %.F90
	$(FC) $(FFLAGS_GFS) $(UTILINCS) -c $*.F90
