#revision history:
#   Aug,2011  J. Wang  compile GOCART only when GOCART_MODE is full

TARGET = ../exe/NEMS.x

include ./conf/configure.nems

INCS = -Iatmos -Iocean -Ichem -Iatmos/share -IENS_Cpl -Iatmos/nmm -Iatmos/fim -Iatmos/gfs -Iatmos/phys -Iatmos/gen

MAIN = MAIN_NEMS.o

OBJS = module_EARTH_INTERNAL_STATE.o module_EARTH_GRID_COMP.o \
       module_NEMS_INTERNAL_STATE.o module_NEMS_GRID_COMP.o

LIBS = atmos/share/libshare.a           \
       atmos/phys/libphys.a             \
       atmos/nmm/libnmm.a               \
       atmos/gfs/libutil/gfs_libutil.a  \
       atmos/gfs/io/gfs_io.a            \
       atmos/gfs/phys/gfs_physics.a     \
       atmos/gfs/dyn/gfs_dynamics.a     \
       atmos/gfs/libgfs.a               \
       atmos/fim/libfimtop.a            \
       atmos/gen/libgen.a               \
       atmos/libatmos.a                 \
       atmos/post/libpost.a             \
       ENS_Cpl/ENS_Cpl.a

LIBS_NMM = atmos/libatmos.a             \
           atmos/nmm/libnmm.a           \
           atmos/gfs/libgfs.a           \
           atmos/fim/libfimtop.a        \
           atmos/gen/libgen.a           \
           atmos/phys/libphys.a         \
           atmos/share/libshare.a       \
           atmos/post/libpost.a         \
           ENS_Cpl/ENS_Cpl.a

LIBS_GEN = atmos/libatmos.a             \
           atmos/nmm/libnmm.a           \
           atmos/gfs/libgfs.a           \
           atmos/fim/libfimtop.a        \
           atmos/gen/libgen.a           \
           atmos/phys/libphys.a         \
           atmos/share/libshare.a       \
           atmos/post/libpost.a         \
           ENS_Cpl/ENS_Cpl.a

#TODO:  DRY these settings, they are redundant with atmos/fim/makefile and 
#TODO:  atmos/fim/FIMsrc ...  
SMS      = /gpfs/t3/global/save/wx20tbh/SMS/SMS_r237
SMSLIB   = $(SMS)/lib/libsms.a
# libraries created during build inside atmos/fim/FIM_src_devccs
FIMSRCLIBS = atmos/fim/FIMsrc_devccs/lib/libfim.a      \
             atmos/fim/FIMsrc_devccs/lib/libwrfio.a    \
             atmos/fim/FIMsrc_devccs/lib/libgribio.a   \
             atmos/fim/FIMsrc_devccs/lib/libvlint.a    \
             atmos/fim/FIMsrc_devccs/lib/libsp_4.a     \
             atmos/fim/FIMsrc_devccs/lib/libslint.a    \
             atmos/fim/FIMsrc_devccs/lib/libsigio_4.a  \
             atmos/fim/FIMsrc_devccs/lib/libsfcio_4.a  \
             atmos/fim/FIMsrc_devccs/lib/libw3_4.a     \
             atmos/fim/FIMsrc_devccs/lib/libbacio_4.a  \
             atmos/fim/FIMsrc_devccs/lib/libcntl.a     \
             atmos/fim/FIMsrc_devccs/lib/libsysshare.a

LIBS_FIM = atmos/libatmos.a                     \
           atmos/nmm/libnmm.a                   \
           atmos/gfs/libgfs.a                   \
           atmos/fim/libfimtop.a                \
           $(SMSLIB) $(FIMSRCLIBS)              \
           atmos/gen/libgen.a                   \
           atmos/share/libshare.a               \
           ENS_Cpl/ENS_Cpl.a                    \
           -L/nwprod/lib -lnetcdf_64

MAKEFILE = makefile

#
# GOCART specific compilation variables/flags
#
GOCART_MODE=stub
ESMADIR=chem/gocart
include $(ESMADIR)/src/Config/NCEP_base.mk

#
# log info
#
WHOAMI  := $(shell whoami)
CMP_DATE := $(shell date)
SVN_REV := $(shell svnversion -n .)
SVN_URL := $(shell svn info | grep URL | cut -f 2 -d ' ')
CMP_YEAR := $(shell date +"%Y" )
CMP_JD := $(shell date +"%j" )

CPPFLAGS += -D'SVN_INFO="($(WHOAMI)) $(CMP_DATE) r$(SVN_REV) $(SVN_URL)"'

ifdef CMP_YEAR
  CPPFLAGS += -D'CMP_YEAR=$(CMP_YEAR)'
endif
ifdef CMP_JD
  CPPFLAGS += -D'CMP_JD=$(CMP_JD)'
endif

#for gocart:
ifeq ($(GOCART_MODE),stub)
	LIBGOCART=
else
	LIBGOCART=$(LIB_GOCART)
endif
export GOCART_MODE

################################################################################
#
# Ask user which target he/she wants to build
#
################################################################################
all:
	@echo "gmake nmm_gfs --- for both NMM and GFS core - without the GOCART"
	@echo " or"
	@echo "gmake gen --- for GEN core only"
	@echo " or"
	@echo "gmake nmm_gfs_gen_post GOCART_MODE=full --- for GEN, NMM and GFS cores - with the GOCART and post"
	@echo " or"
	@echo "gmake nmm_gfs GOCART_MODE=full --- for both NMM and GFS core - with the GOCART"
	@echo " or"
	@echo "gmake nmm --- for NMM core only"
	@echo " or"
	@echo "gmake nmm_post --- for NMM core only - with post"
	@echo " or"
	@echo "gmake gfs --- for GFS core only - without the GOCART"
	@echo " or"
	@echo "gmake gfs_post --- for GFS core only - with post without the GOCART"
	@echo " or"
	@echo "gmake gfs GOCART_MODE=full --- for GFS core only - with the GOCART"
	@echo " or"
	@echo "gmake fim --- for FIM core only"
#	@echo " or"
#	@echo "gmake fim_gfs --- for both FIM and GFS core only - without the GOCART"
#	@echo " or"
#	@echo "gmake fim_gfs GOCART_MODE=full --- for both FIM and GFS core only - with the GOCART"


################################################################################
#
# NMM
#
################################################################################

nmm: nmm_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_NMM) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS_NMM) $(EXTLIBS) 
	@echo "$(TARGET) is created for NMM core."

nmm_libs:
	cd atmos/share         && gmake
	cd atmos/phys          && gmake gocart_stub
	cd atmos/nmm           && gmake
	cd atmos/gfs           && gmake stub
	cd atmos/gen           && gmake stub
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake stub

################################################################################
#
# NMM with post
#
################################################################################

nmm_post: nmm_post_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_NMM) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS_NMM) $(EXTLIBS_POST) 
	@echo "$(TARGET) is created for NMM core."

nmm_post_libs:
	cd atmos/share         && gmake
	cd atmos/phys          && gmake gocart_stub
	cd atmos/nmm           && gmake
	cd atmos/gfs           && gmake stub
	cd atmos/gen           && gmake stub
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake gfs_stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake stub

################################################################################
#
# GFS
#
################################################################################

gfs: gfs_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS) $(EXTLIBS) $(LIBGOCART)
	@echo "nems.x is created for GFS core."

gfs_libs:
ifeq ($(GOCART_MODE),full)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
endif
	cd atmos/share         && gmake
	cd atmos/phys          && gmake gocart_$(GOCART_MODE)
	cd atmos/nmm           && gmake stub
	cd atmos/gen           && gmake stub
ifeq ($(GOCART_MODE),full)
	cd atmos/gfs           && gmake gocart_full
else
	cd atmos/gfs           && gmake gocart_stub
endif
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# GFS with post
#
################################################################################

gfs_post: gfs_post_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS) $(EXTLIBS_POST) $(LIBGOCART)
	@echo "nems.x is created for GFS core."

gfs_post_libs:
ifeq ($(GOCART_MODE),full)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
endif
	cd atmos/share         && gmake
	cd atmos/phys          && gmake gocart_$(GOCART_MODE)
	cd atmos/nmm           && gmake stub
	cd atmos/gen           && gmake stub
ifeq ($(GOCART_MODE),full)
	cd atmos/gfs           && gmake gocart_full
else
	cd atmos/gfs           && gmake gocart_stub
endif
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake nmm_stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# GEN
#
################################################################################

gen: gen_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS_GEN) $(EXTLIBS)
	@echo "nems.x is created for GEN core."

gen_libs:
	cd atmos/share         && gmake
	cd atmos/phys          && gmake
	cd atmos/nmm           && gmake stub
	cd atmos/gfs           && gmake stub
	cd atmos/fim           && gmake stub
	cd atmos/gen           && gmake
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake


################################################################################
#
# GEN, NMM and GFS
#
################################################################################

nmm_gfs_gen: nmm_gfs_gen_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS) $(EXTLIBS) $(LIBGOCART)
	@echo "nems.x is created for GEN, NMM and GFS core."

nmm_gfs_gen_libs:
ifeq ($(GOCART_MODE),full)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
endif
	cd atmos/share         && gmake
	cd atmos/phys          && gmake gocart_$(GOCART_MODE)
	cd atmos/nmm           && gmake
ifeq ($(GOCART_MODE),full)
	cd atmos/gfs           && gmake gocart_full
else
	cd atmos/gfs           && gmake gocart_stub
endif
	cd atmos/gen           && gmake 
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# GEN, NMM and GFS with post
#
################################################################################

nmm_gfs_gen_post: nmm_gfs_gen_post_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS) $(EXTLIBS_POST) $(LIBGOCART)
	@echo "nems.x is created for GEN, NMM and GFS core."

nmm_gfs_gen_post_libs:
ifeq ($(GOCART_MODE),full)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
endif
	cd atmos/share         && gmake
	cd atmos/phys          && gmake gocart_$(GOCART_MODE)
	cd atmos/nmm           && gmake
ifeq ($(GOCART_MODE),full)
	cd atmos/gfs           && gmake gocart_full
else
	cd atmos/gfs           && gmake gocart_stub
endif
	cd atmos/gen           && gmake
	cd atmos/fim           && gmake stub
	cd atmos/post          && gmake 
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# NMM and GFS
#
################################################################################

nmm_gfs: nmm_gfs_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS) $(EXTLIBS) $(LIBGOCART)
	@echo "nems.x is created for NMM and GFS core."

nmm_gfs_libs:
ifeq ($(GOCART_MODE),full)
	cd chem/gocart/src     && gmake install GOCART_MODE=$(GOCART_MODE)
endif
	cd atmos/share         && gmake
	cd atmos/phys          && gmake gocart_$(GOCART_MODE)
	cd atmos/nmm           && gmake
ifeq ($(GOCART_MODE),full)
	cd atmos/gfs           && gmake gocart_full
else
	cd atmos/gfs           && gmake gocart_stub
endif
	cd atmos/fim           && gmake stub
	cd atmos/gen           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# FIM
#
################################################################################

# TODO:  Change to re-use NEMS standard $(EXTLIBS) instead of FIM versions 
# TODO:  currently in $(LIBS_FIM) .  
fim: fim_libs $(OBJS) $(MAIN)
	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS_FIM) $(ESMF_LIBS)
	@echo "nems.x is created for FIM core."

fim_libs:
	cd atmos/share         && gmake
	cd atmos/nmm           && gmake stub
	cd atmos/gfs           && gmake stub 
	cd atmos/fim           && gmake
	cd atmos/gen           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# FIM and GFS
#TBH:  not yet integrated with FIM build
#
################################################################################

fim_gfs: fim_gfs_libs $(OBJS) $(MAIN)
#	$(FC) $(FFLAGS_GFS) $(SMP) -o $(TARGET) $(MAIN) $(OBJS) $(LIBS) $(EXTLIBS)
#	@echo "nems.x is created for FIM and GFS core."
	@echo "ERROR:  Cannot build FIM with other cores yet.  nems.x is NOT created for FIM and GFS core."
	@exit -1

fim_gfs_libs:
	cd atmos/share         && gmake
	cd atmos/nmm           && gmake stub
	cd atmos/gfs           && gmake
	cd atmos/fim           && gmake
	cd atmos/gen           && gmake stub
	cd atmos/post          && gmake stub
	cd atmos               && gmake
	cd ENS_Cpl             && gmake

################################################################################
#
# Clean everything
#
################################################################################

clean:
	cd atmos               && gmake clean
	cd ENS_Cpl             && gmake clean
	cd chem                && gmake clean GOCART_MODE=full
	$(RM) -f *.f90 *.lst *.o *.mod lm map $(TARGET)

.SUFFIXES:
.SUFFIXES: .F90 .f90 .o

.F90.f90:
	$(CPP) $(CPPFLAGS) $< > $*.f90

$(OBJS): %.o: %.f90
	$(FC) $(FFLAGS) $(INCS) -c $*.f90

MAIN_NEMS.o: %.o: %.f90
	$(FC) $(FFLAGS) $(INCS) -c $*.f90
	$(RM) $*.f90
