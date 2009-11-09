#
# Makefile for ESMA components.
#

# Make sure ESMADIR is defined
# ----------------------------
ifndef ESMADIR
       ESMADIR := $(PWD)/..
endif


# Compilation rules, flags, etc
# -----------------------------
  include $(ESMADIR)/Config/ESMA_base.mk  # Generic stuff
  include $(ESMADIR)/Config/ESMA_arch.mk  # System dependencies
  include $(ESMADIR)/Config/GMAO_base.mk  # GMAO specific macros

#                  ---------------------
#                  Standard ESMA Targets
#                  ---------------------

esma_help :
	@echo "Standard ESMA targets:"
	@echo "% make esma_install    (builds and install under ESMADIR)"
	@echo "% make esma_clean      (removes deliverables: *.[aox], etc)"
	@echo "% make esma_distclean  (leaves in the same state as cvs co)"
	@echo "% make esma_doc        (generates PDF, installs under ESMADIR)"
	@echo "% make esma_help       (this message)"
	@echo "Environment:"
	@echo "      ESMADIR = $(ESMADIR)"
	@echo "      BASEDIR = $(BASEDIR)"
	@echo "         ARCH = $(ARCH)"
	@echo "         SITE = $(SITE) "


THIS := $(shell basename `pwd`)
LIB   = $(LIB_THIS)


#                  --------------------------------
#                   Recurse Make in Sub-directories
#                  --------------------------------

# Handle stubs for NCEP implementarion of GOCART Stand alone
# ----------------------------------------------------------
  ifeq ($(GOCART_MODE),stub)
     ALLDIRS = 
     SRCS := GOCART_StubCompMod.F90
  else
     ALLDIRS = O3_GridComp CO_GridComp CO2_GridComp CFC_GridComp CARMA_GridComp \
               DU_GridComp SS_GridComp BC_GridComp OC_GridComp SU_GridComp Rn_GridComp 
     SRCS := GOCART_GridCompMod.F90 Aero_GridCompMod.F90  
  endif

SUBDIRS = $(wildcard $(ALLDIRS))

TARGETS = esma_install esma_clean esma_distclean esma_doc \
          install clean distclean doc 

DEBUG = # $(D)DEBUG

ifeq ($(ESMA_REAL),$(FREAL8))
	FREAL = $(FREAL8)
else
	FREAL = $(FREAL4)
endif

export ESMADIR BASEDIR ARCH SITE FREAL DEBUG

$(TARGETS): 
	@ t=$@; argv="$(SUBDIRS)" ;\
	  for d in $$argv; do			 \
	    ( cd $$d				;\
	      echo ""; echo Making $$t in `pwd`          ;\
	      $(MAKE) -e $$t ) \
	  done
	$(MAKE) local_$@

local_esma_install local_install: $(LIB)
	$(MKDIR) -p $(INC_CHEM)
	$(CP) -p *.mod          $(INC_CHEM)

local_esma_clean local_clean:
	-$(RM) *~ *.[aox] *.[Mm][Oo][Dd]

local_esma_distclean local_distclean:
	-$(RM) *~ *.[aoxd] *.[Mm][Oo][Dd] *___.[fF]*

local_esma_doc local_doc:
	@echo "Target $@ not implemented yet in `pwd`"


#                  --------------------
#                  User Defined Targets
#                  --------------------

OBJS := $(addsuffix .o, $(basename $(SRCS))) 
DEPS := $(addsuffix .d, $(basename $(SRCS))) 

INC_DIRS = . $(INC_MAPL_BASE) $(INC_CHEM)

MOD_DIRS = . $(INC_ESMF) $(INC_MAPL_BASE) $(INC_GEOS_SHARED) \
             $(INC_CHEM_BASE) $(INC_CHEM_SHARED) $(INC_CHEM) $(INC_MPEU)

USER_FINCS  = $(foreach dir,$(INC_DIRS),$(I)$(dir)) 
USER_FMODS  = $(foreach dir,$(MOD_DIRS),$(M)$(dir)) 
USER_FFLAGS = $(BIG_ENDIAN) $(DEBUG)

vpath % $(MOD_DIRS)

$(LIB) lib : $(OBJS)
	$(AR) $(AR_FLAGS) $(LIB) $(OBJS)

# Hack to prevent remaking dep files during cleaning
# --------------------------------------------------
  ifneq ($(findstring clean,$(MAKECMDGOALS)),clean)
    -include $(DEPS)
  endif

#.

  -include $(ESMADIR)/Config/ESMA_post.mk  # ESMA additional targets, macros

