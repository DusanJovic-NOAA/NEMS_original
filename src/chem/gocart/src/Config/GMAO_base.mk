#
# Earth System Modeling Applications (ESMA) base makefile fragment.
# This fragment defines MAO specific macros associated with basic, shared 
# packages. This file must be included after ESMA_base.mk and ESMA_arch.mk.
#
# REVISION HISTORY:
#
# 13oct04  da Silva  First Crack
# 04apr05  Todling   Commented out REAL4 definition
#
#--------------------------------------------------------------------------

#                       ----------------
#                        Shared Packages
#                       ----------------

# GEOS-5 Build environment
# ------------------------

DOING_GEOS5 = TRUE

THIS_MPEU = GMAO_mpeu
INC_MPEU = $(ESMAINC)/GMAO_mpeu
LIB_MPEU = $(ESMALIB)/libGMAO_mpeu.a
LIB_EU   = $(ESMALIB)/libGMAO_eu.a

INC_ODS = $(ESMAINC)/GMAO_ods
LIB_ODS = $(ESMALIB)/libGMAO_ods.a

ifeq ($(ESMA_SDF),hdf)
   INC_MFHDF3 = $(ESMAINC)/GMAO_mfhdf3
   LIB_MFHDF3 = $(ESMALIB)/libGMAO_mfhdf3.a
else
   INC_MFHDF3 = $(INC_NETCDF)
   LIB_MFHDF3 = 
endif

THIS_GFIO = GMAO_gfio
INC_GFIO = $(ESMAINC)/$(THIS_GFIO)
LIB_GFIO = $(ESMALIB)/lib$(THIS_GFIO).a

THIS_CFIO = MAPL_cfio
INC_CFIO = $(ESMAINC)/$(THIS_CFIO)
LIB_CFIO = $(ESMALIB)/lib$(THIS_CFIO).a

THIS_CFIO_EOS = MAPL_cfio_eos
INC_CFIO_EOS = $(ESMAINC)/$(THIS_CFIO_EOS)
LIB_CFIO_EOS = $(ESMALIB)/lib$(THIS_CFIO_EOS).a

INC_GFIOEOS = $(ESMAINC)/GMAO_gfioeos
LIB_GFIOEOS = $(ESMALIB)/libGMAO_gfioeos.a $(LIB_EOS)

THIS_HERMES = GMAO_hermes
INC_HERMES = $(ESMAINC)/GMAO_hermes
LIB_HERMES = $(ESMALIB)/libGMAO_hermes.a

INC_CHEM_BASE = $(ESMAINC)/Chem_Base
LIB_CHEM_BASE = $(ESMALIB)/libChem_Base.a

INC_CHEM_SHARED = $(ESMAINC)/Chem_Shared
LIB_CHEM_SHARED = $(ESMALIB)/libChem_Shared.a

INC_CHEM = $(ESMAINC)/GEOSchem_GridComp
LIB_CHEM = $(ESMALIB)/libGEOSchem_GridComp.a

INC_GEOS_BASE = $(ESMAINC)/GEOS_Base
LIB_GEOS_BASE = $(ESMALIB)/libGEOS_Base.a

# This should be renamed MAPL_BASE
INC_MAPL_BASE = $(ESMAINC)/MAPL_Base
LIB_MAPL_BASE = $(ESMALIB)/libMAPL_Base.a

INC_GEOS_SHARED = $(ESMAINC)/GEOS_Shared
LIB_GEOS_SHARED = $(ESMALIB)/libGEOS_Shared.a

INC_PILGRIM = $(ESMAINC)/GMAO_pilgrim
LIB_PILGRIM = $(ESMALIB)/libGMAO_pilgrim.a

INC_TRANSF = $(ESMAINC)/GMAO_transf
LIB_TRANSF = $(ESMALIB)/libGMAO_transf.a

INC_PSAS = $(ESMAINC)/GMAO_psas
LIB_PSAS = $(ESMALIB)/libGMAO_psas.a

THIS_GSI = GSI_GridComp
INC_GSI = $(ESMAINC)/$(THIS_GSI)
LIB_GSI = $(ESMALIB)/lib$(THIS_GSI).a

THIS_GEOSGSI_COUPLER = GEOSgsi_Coupler
INC_GEOSGSI_COUPLER = $(ESMAINC)/$(THIS_GEOSGSI_COUPLER)
LIB_GEOSGSI_COUPLER = $(ESMALIB)/lib$(THIS_GEOSGSI_COUPLER).a

THIS_BACIO = NCEP_bacio
INC_BACIO  = $(ESMAINC)/$(THIS_BACIO)
LIB_BACIO  = $(ESMALIB)/lib$(THIS_BACIO).a

THIS_BUFR = NCEP_bufr
INC_BUFR  = $(ESMAINC)/$(THIS_BUFR)
LIB_BUFR  = $(ESMALIB)/lib$(THIS_BUFR).a

THIS_CRTM = NCEP_crtm
INC_CRTM  = $(ESMAINC)/$(THIS_CRTM)
LIB_CRTM  = $(ESMALIB)/lib$(THIS_CRTM).a

THIS_IRSSE = NCEP_irsse
INC_IRSSE  = $(ESMAINC)/$(THIS_IRSSE)
LIB_IRSSE  = $(ESMALIB)/lib$(THIS_IRSSE).a

THIS_SFCIO = NCEP_sfcio
INC_SFCIO  = $(ESMAINC)/$(THIS_SFCIO)
LIB_SFCIO  = $(ESMALIB)/lib$(THIS_SFCIO).a

THIS_SIGIO = NCEP_sigio
INC_SIGIO  = $(ESMAINC)/$(THIS_SIGIO)
LIB_SIGIO  = $(ESMALIB)/lib$(THIS_SIGIO).a

THIS_SP = NCEP_sp
INC_SP  = $(ESMAINC)/$(THIS_SP)
LIB_SP  = $(ESMALIB)/lib$(THIS_SP).a

THIS_W3 = NCEP_w3
INC_W3  = $(ESMAINC)/$(THIS_W3)
LIB_W3  = $(ESMALIB)/lib$(THIS_W3).a

INC_AIRS = $(ESMAINC)/GMAO_iret--airs
LIB_AIRS = $(ESMALIB)/libGMAO_iret--airs.a

INC_FASTEM = $(ESMAINC)/GMAO_iret--fastem
LIB_FASTEM = $(ESMALIB)/libGMAO_iret--fastem.a

INC_GLATOVS = $(ESMAINC)/GMAO_iret--glatovs
LIB_GLATOVS = $(ESMALIB)/libGMAO_iret--glatovs.a

INC_HFFP = $(ESMAINC)/GMAO_iret--hffp
LIB_HFFP = $(ESMALIB)/libGMAO_iret--hffp.a

INC_IDL = $(ESMAINC)/GMAO_iret--idl
LIB_IDL = $(ESMALIB)/libGMAO_iret--idl.a

INC_IRUTIL = $(ESMAINC)/GMAO_iret--irutil
LIB_IRUTIL = $(ESMALIB)/libGMAO_iret--irutil.a

INC_MIT = $(ESMAINC)/GMAO_iret--mit
LIB_MIT = $(ESMALIB)/libGMAO_iret--mit.a

INC_OPTRAN = $(ESMAINC)/GMAO_iret--optran
LIB_OPTRAN = $(ESMALIB)/libGMAO_iret--optran.a

INC_RADTRANS = $(ESMAINC)/GMAO_iret--radtrans
LIB_RADTRANS = $(ESMALIB)/libGMAO_iret--radtrans.a

INC_SARTA = $(ESMAINC)/GMAO_iret--sarta
LIB_SARTA = $(ESMALIB)/libGMAO_iret--sarta.a


#
# The modules below are probably not needed.
#

INC_OBS = $(ESMAINC)/obs
LIB_OBS = $(ESMALIB)/libobs.a

INC_ANA = $(ESMAINC)/ana
LIB_ANA = $(ESMALIB)/libana.a

INC_FVGCM = $(ESMAINC)/fvgcm
LIB_FVGCM = $(ESMALIB)/libfvgcm.a


