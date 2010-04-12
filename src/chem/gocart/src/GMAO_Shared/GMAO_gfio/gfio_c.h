#include "cfortran.h"
#include <netcdf.h>

PROTOCCALLSFSUB4(GFIO_OPEN,gfio_open,STRING,INT,PINT,PINT)
#define Gfio_Open(A1,A2,A3,A4) CCALLSFSUB4(GFIO_OPEN,gfio_open,STRING,INT,PINT,PINT,\
                                        A1,A2,A3,A4)

PROTOCCALLSFSUB25(GFIO_CREATE,gfio_create,\
                   STRING,STRING,STRING,STRING,FLOAT,\
                   INT,INT,INT,FLOATV,FLOATV,FLOATV,STRING,\
                   INT,INT,INT,\
                   INT,STRINGV,STRINGV,STRINGV,INTV,\
                   FLOATVV,FLOATVV,INT,\
                   PINT,PINT)
#define Gfio_Create(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,\
                    A16,A17,A18,A19,A20,A21,A22,A23,A24,A25) \
                   CCALLSFSUB25(GFIO_CREATE,gfio_create,\
                   STRING,STRING,STRING,STRING,FLOAT,\
                   INT,INT,INT,FLOATV,FLOATV,FLOATV,STRING,\
                   INT,INT,INT,\
                   INT,STRINGV,STRINGV,STRINGV,INTV,\
                   FLOATVV,FLOATVV,INT,\
                   PINT,PINT,\
                   A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,\
                   A16,A17,A18,A19,A20,A21,A22,A23,A24,A25)

PROTOCCALLSFSUB2(GFIO_CLOSE,gfio_close,INT,PINT)
#define Gfio_Close(A1,A2) CCALLSFSUB2(GFIO_CLOSE,gfio_close,INT,PINT,A1,A2)

PROTOCCALLSFSUB10(GFIO_PUTVAR,gfio_putvar,INT,STRING,INT,INT,\
                  INT,INT,INT,INT,FLOATVV,PINT)
#define Gfio_PutVar(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10) CCALLSFSUB10(GFIO_PUTVAR,\
                  gfio_putvar,INT,STRING,INT,INT,\
                  INT,INT,INT,INT,FLOATVV,PINT,A1,A2,A3,A4,A5,A6,\
                  A7,A8,A9,A10)

PROTOCCALLSFSUB10(GFIO_GETVAR,gfio_getvar,INT,STRING,INT,INT,\
                  INT,INT,INT,INT,FLOATVV,PINT)
#define Gfio_GetVar(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10) CCALLSFSUB10(GFIO_GETVAR,\
                  gfio_getvar,INT,STRING,INT,INT,\
                  INT,INT,INT,INT,FLOATVV,PINT,A1,A2,A3,A4,A5,A6,\
                  A7,A8,A9,A10)

PROTOCCALLSFSUB8(GFIO_DIMINQUIRE,gfio_diminquire,INT,PINT,PINT,PINT,PINT,\
                 PINT,PINT,PINT)
#define Gfio_DimInquire(A1,A2,A3,A4,A5,A6,A7,A8) CCALLSFSUB8(GFIO_DIMINQUIRE,\
                 gfio_diminquire,INT,PINT,PINT,PINT,PINT,PINT,PINT,PINT,\
                 A1,A2,A3,A4,A5,A6,A7,A8)

PROTOCCALLSFSUB24(GFIO_INQUIRE,gfio_inquire,INT,PINT,PINT,PINT,PINT,PINT,\
                  PSTRING,PSTRING,PSTRING,PFLOAT,PVOID,PVOID,PVOID,\
                  PSTRING,PVOID,PVOID,PINT,PSTRINGV,PSTRINGV,PSTRINGV,PVOID,\
                  PVOID,PVOID,PINT)
#define Gfio_Inquire(A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,\
                     A16,A17,A18,A19,A20,A21,A22,A23,A24) \
        CCALLSFSUB24(GFIO_INQUIRE,gfio_inquire,INT,PINT,PINT,PINT,PINT,PINT,\
                  PSTRING,PSTRING,PSTRING,PFLOAT,\
                  PVOID,PVOID,PVOID,PSTRING,\
                  PVOID,PVOID,PINT,\
                  PSTRINGV,PSTRINGV,PSTRINGV,PVOID,\
                  PVOID,PVOID,PINT,A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,\
                  A13,A14,A15,A16,A17,A18,A19,A20,A21,A22,A23,A24)

PROTOCCALLSFSUB6(GFIO_PUTINTATT,gfio_putintatt,INT,STRING,INT,INTV,INT,PINT)
#define Gfio_PutIntAtt(A1,A2,A3,A4,A5,A6) CCALLSFSUB6(GFIO_PUTINTATT,\
        gfio_putintatt,INT,STRING,INT,INTV,INT,PINT,A1,A2,A3,A4,A5,A6)

PROTOCCALLSFSUB6(GFIO_PUTREALATT,gfio_putrealatt,INT,STRING,INT,FLOATV,INT,PINT)
#define Gfio_PutRealAtt(A1,A2,A3,A4,A5,A6) CCALLSFSUB6(GFIO_PUTREALATT,\
        gfio_putrealatt,INT,STRING,INT,FLOATV,INT,PINT,A1,A2,A3,A4,A5,A6)

PROTOCCALLSFSUB5(GFIO_PUTCHARATT,gfio_putcharatt,INT,STRING,INT,STRING,PINT)
#define Gfio_PutCharAtt(A1,A2,A3,A4,A5) CCALLSFSUB5(GFIO_PUTCHARATT,\
        gfio_putcharatt,INT,STRING,INT,STRING,PINT,A1,A2,A3,A4,A5)

PROTOCCALLSFSUB4(GFIO_GETATTNAMES,gfio_getattnames,INT,PINT,PSTRINGV,PINT)
#define Gfio_GetAttNames(A1,A2,A3,A4) CCALLSFSUB4(GFIO_GETATTNAMES,\
        gfio_getattnames,INT,PINT,PSTRINGV,PINT,A1,A2,A3,A4)

PROTOCCALLSFSUB5(GFIO_ATTINQUIRE,gfio_attinquire,INT,STRING,PINT,PINT,PINT)
#define Gfio_AttInquire(A1,A2,A3,A4,A5) CCALLSFSUB5(GFIO_ATTINQUIRE,\
        gfio_attinquire,INT,STRING,PINT,PINT,PINT,A1,A2,A3,A4,A5)

PROTOCCALLSFSUB5(GFIO_GETINTATT,gfio_getintatt,INT,STRING,PINT,PVOID,PINT)
#define Gfio_GetIntAtt(A1,A2,A3,A4,A5) CCALLSFSUB5(GFIO_GETINTATT,\
        gfio_getintatt,INT,STRING,PINT,PVOID,PINT,A1,A2,A3,A4,A5)

PROTOCCALLSFSUB5(GFIO_GETREALATT,gfio_getrealatt,INT,STRING,PINT,PVOID,PINT)
#define Gfio_GetRealAtt(A1,A2,A3,A4,A5) CCALLSFSUB5(GFIO_GETREALATT,\
        gfio_getrealatt,INT,STRING,PINT,PVOID,PINT,A1,A2,A3,A4,A5)

PROTOCCALLSFSUB5(GFIO_GETCHARATT,gfio_getcharatt,INT,STRING,PINT,PSTRING,PINT)
#define Gfio_GetCharAtt(A1,A2,A3,A4,A5) CCALLSFSUB5(GFIO_GETCHARATT,\
        gfio_getcharatt,INT,STRING,PINT,PSTRING,PINT,A1,A2,A3,A4,A5)
