#if 0
//
// Make this header file available as ESMFVersionDefine.h in order to build
// NEMS against a public ESMF 5 installation (i.e. 5.2.0rp1 and up).
//
#endif

#undef ESMF_3

#ifndef ESMF_MAJOR_VERSION
#define ESMF_MAJOR_VERSION 5
#define ESMF_MINOR_VERSION 2
#endif

#include "./ESMFVersionLogic.h"
