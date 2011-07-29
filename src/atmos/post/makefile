include ../../conf/configure.nems

MAKEFILE = makefile

UTILINCSPOST = -I../share -I../nmm/ -I../gfs/io -I/climate/save/wx20wa/esmf/nems/unipost_r14561/incmod/post_4
UTILINCSPOSTGFS = -I../share -I../gfs/io -I/climate/save/wx20wa/esmf/nems/unipost_r14561/incmod/post_4
UTILINCSPOSTNMM = -I../share -I../nmm -I/climate/save/wx20wa/esmf/nems/unipost_r14561/incmod/post_4
UTILINCSPOSTSTUB = -I../share

LIBRARY  = libpost.a

MODULE   = post_gfs.o post_nmm.o post_nems_routines.o

MODULES_STUB  = post_gfs_stub.o post_nmm_stub.o post_nems_routines_stub.o
MODULES_GFS_STUB  = post_nmm.o post_gfs_stub.o post_nems_routines.o
MODULES_NMM_STUB  = post_nmm_stub.o post_gfs.o post_nems_routines.o

DEPEND_FILES = ${MODULE_POST:.o=.f}

all default: depend
	@gmake -f $(MAKEFILE) $(LIBRARY)

$(LIBRARY): $(MODULE) 
	$(AR) $(ARFLAGS) $@ $(MODULE)

stub: $(MODULES_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_STUB)
	
gfs_stub: $(MODULES_GFS_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_GFS_STUB)
	
nmm_stub: $(MODULES_NMM_STUB)
	$(AR) $(ARFLAGS) $(LIBRARY) $(MODULES_NMM_STUB)
	
$(MODULE): %.o: %.f
	$(FC) $(FFLAGS) $(TRAPS) $(UTILINCSPOST) -c $*.f

$(MODULE_GFS_STUB): %.o: %.f
	$(FC) $(FFLAGS) $(TRAPS) $(UTILINCSPOSTNMM) -c $*.f

$(MODULE_NMM_STUB): %.o: %.f
	$(FC) $(FFLAGS) $(TRAPS) $(UTILINCSPOSTGFS) -c $*.f

$(MODULES_STUB): %.o: %.f
	$(FC) $(FFLAGS) $(UTILINCSPOSTSTUB) -c $*.f

clean:
	$(RM) -f $(LIBRARY) *.o *.mod *.lst lm map depend

include ../../conf/make.rules

include depend
