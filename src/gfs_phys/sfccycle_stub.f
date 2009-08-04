      SUBROUTINE SFCCYCLE(LUGB,LEN,LSOIL,SIG1T,DELTSFC &
     ,                   IY,IM,ID,IH,FH &
     ,                   RLA, RLO, SLMASK,OROG &
     ,                   SIHFCS,SICFCS,SITFCS  &               
     ,                   SWDFCS,SLCFCS   &   
     ,                   VMNFCS,VMXFCS,SLPFCS,ABSFCS &
     ,                   TSFFCS,SNOFCS,ZORFCS,ALBFCS,TG3FCS &
     ,                   CNPFCS,SMCFCS,STCFCS,SLIFCS,AISFCS,F10M &
     ,                   VEGFCS,VETFCS,SOTFCS,ALFFCS &
     ,                   CVFCS,CVBFCS,CVTFCS,me,NLUNIT,IALB)
!
      implicit none

      INTEGER imsk,jmsk,ifp,irtscv,irtacn,irtais,irtsno,irtzor &
     ,        irtalb,irtsot,irtalf,j,irtvet,irtsmc,irtstc,irtveg &
     ,        irtwet,k,iprnt,kk,irttsf,iret,i,igrdbg,iy,im,id &
     ,        icalbl,icalbs,icalfl,ictsfs,lugb,len,lsoil,ih &
     ,        ictsfl,iczors,icplrl,icplrs,iczorl,icalfs,icsnol &
     ,        icsnos,irttg3,me,KQCM, NLUNIT,IALB &
     ,       irtvmn, irtvmx, irtslp, irtabs
      LOGICAL LGCHEK
      DATA LGCHEK/.TRUE./
      REAL (KIND=8) SLMASK(LEN),OROG(LEN)
      REAL (KIND=8) RLA(LEN), RLO(LEN)
      CHARACTER*500 FNGLAC,FNMXIC
      real (kind=8), allocatable :: GLACIR(:),AMXICE(:),TSFCL0(:)
      CHARACTER*500 FNTSFC,FNWETC,FNSNOC,FNZORC,FNALBC,FNAISC &
     ,              FNPLRC,FNTG3C,FNSCVC,FNSMCC,FNSTCC,FNACNC &
     ,              FNVEGC,fnvetc,fnsotc &
     ,             FNVMNC,FNVMXC,FNSLPC,FNABSC, FNALBC2 
      REAL (KIND=8) TSFCLM(LEN), WETCLM(LEN),   SNOCLM(LEN) &
     ,     ZORCLM(LEN), ALBCLM(LEN,4), AISCLM(LEN) &
     ,     TG3CLM(LEN), ACNCLM(LEN),   CNPCLM(LEN) &
     ,     CVCLM (LEN), CVBCLM(LEN),   CVTCLM(LEN) &
     ,     SCVCLM(LEN), TSFCL2(LEN),   VEGCLM(LEN) &
     ,     vetclm(LEN), sotclm(LEN),   ALFCLM(LEN,2), SLICLM(LEN) &
     ,     SMCCLM(LEN,LSOIL), STCCLM(LEN,LSOIL) &
     ,    SIHCLM(LEN), SICCLM(LEN) &
     ,    VMNCLM(LEN), VMXCLM(LEN), SLPCLM(LEN), ABSCLM(LEN)
      CHARACTER*500 FNTSFA,FNWETA,FNSNOA,FNZORA,FNALBA,FNAISA &
     ,             FNPLRA,FNTG3A,FNSCVA,FNSMCA,FNSTCA,FNACNA &
     ,             FNVEGA,fnveta,fnsota &
     ,            FNVMNA,FNVMXA,FNSLPA,FNABSA       
!
      REAL (KIND=8) TSFANL(LEN), WETANL(LEN),   SNOANL(LEN) &
     ,     ZORANL(LEN), ALBANL(LEN,4), AISANL(LEN) &
     ,     TG3ANL(LEN), ACNANL(LEN),   CNPANL(LEN) &
     ,     CVANL (LEN), CVBANL(LEN),   CVTANL(LEN) &
     ,     SCVANL(LEN), TSFAN2(LEN),   VEGANL(LEN) &
     ,     vetanl(LEN), sotanl(LEN),   ALFANL(LEN,2), SLIANL(LEN) &
     ,     SMCANL(LEN,LSOIL), STCANL(LEN,LSOIL) &
     ,    SIHANL(LEN), SICANL(LEN) &
     ,    VMNANL(LEN), VMXANL(LEN), SLPANL(LEN), ABSANL(LEN)
      REAL (KIND=8) TSFAN0(LEN) !  Sea surface temperature analysis at FT=0.
      REAL (KIND=8) TSFFCS(LEN), WETFCS(LEN),   SNOFCS(LEN) &
     ,     ZORFCS(LEN), ALBFCS(LEN,4), AISFCS(LEN) &
     ,     TG3FCS(LEN), ACNFCS(LEN),   CNPFCS(LEN) &
     ,     CVFCS (LEN), CVBFCS(LEN),   CVTFCS(LEN) &
     ,     SLIFCS(LEN), VEGFCS(LEN) &
     ,     vetfcs(LEN), sotfcs(LEN),   alffcs(LEN,2) &
     ,     SMCFCS(LEN,LSOIL), STCFCS(LEN,LSOIL) &
     ,    SIHFCS(LEN), SICFCS(LEN), SITFCS(LEN) &
     ,    VMNFCS(LEN), VMXFCS(LEN), SLPFCS(LEN), ABSFCS(LEN) &
     ,    SWDFCS(LEN), SLCFCS(LEN,LSOIL)
!
! Ratio of sigma level 1 wind and 10m wind (diagnozed by model and not touched
! in this program).
!
      REAL (KIND=8) F10M  (LEN)
      REAL (KIND=8) FSMCL(25),FSMCS(25),FSTCL(25),FSTCS(25)
      REAL (KIND=8) FCSMCL(25),FCSMCS(25),FCSTCL(25),FCSTCS(25)

      REAL (KIND=8) SWRATIO(LEN,LSOIL)
      REAL (KIND=8) deltsfc, fh
      LOGICAL FIXRATIO(LSOIL)
!
      INTEGER ICSMCL(25), ICSMCS(25), ICSTCL(25), ICSTCS(25)
!
      Integer kpd9
!
      logical icefl1(len), icefl2(len)
!
!  Input and output SURFACE FIELDS (BGES) file names
!
!
!  Sigma level 1 temperature for dead start
!
      REAL (KIND=8) SIG1T(LEN), WRK(LEN)
!
      CHARACTER*32 LABEL
!
!  = 1 ==> FORECAST IS USED
!  = 0 ==> ANALYSIS (OR CLIMATOLOGY) IS USED
!
!     OUTPUT FILE  ... PRIMARY SURFACE FILE FOR RADIATION AND FORECAST
!
!       REC.  1    LABEL
!       REC.  2    DATE RECORD
!       REC.  3    TSF
!       REC.  4    SOILM(TWO LAYERS)              ----> 4 layers
!       REC.  5    SNOW
!       REC.  6    SOILT(TWO LAYERS)              ----> 4 layers
!       REC.  7    TG3
!       REC.  8    ZOR
!       REC.  9    CV
!       REC. 10    CVB
!       REC. 11    CVT
!       REC. 12    ALBEDO (four types)
!       REC. 13    SLIMSK
!       REC. 14    vegetation cover
!       REC. 14    PLANTR                         -----> skip this record
!       REC. 15    F10M                           -----> CANOPY
!       REC. 16    CANOPY WATER CONTENT (CNPANL)  -----> F10M
!       REC. 17    vegetation type
!       REC. 18    soil type
!       REC. 19    zeneith angle dependent vegetation fraction (two types)
!       REC. 20    UUSTAR
!       REC. 21    FFMM
!       REC. 22    FFHH
!Cwu add SIH & SIC
!       REC. 23    SIH(one category only)
!       REC. 24    SIC
!Clu [+8L] add PRCP, FLAG, SWD, SLC, VMN, VMX, SLP, ABS
!       REC. 25    TPRCP
!       REC. 26    SRFLAG
!       REC. 27    SWD
!       REC. 28    SLC (4 LAYERS)
!       REC. 29    VMN
!       REC. 30    VMX
!       REC. 31    SLP
!       REC. 32    ABS

!
!  Debug only
!   LDEBUG=.TRUE. creates BGES files for climatology and analysis
!   LQCBGS=.TRUE. Quality controls input BGES file before merging (should have been
!              QCed in the forecast program)
!
      LOGICAL LDEBUG,LQCBGS
      logical lprnt
!
!  Debug only
!
      CHARACTER*500 FNDCLM,FNDANL
!
      LOGICAL LANOM
!
!
      RETURN
      END SUBROUTINE SFCCYCLE
