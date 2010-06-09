!----------------------------------------------------------------------------
module module_nemsio_mpi
!$$$ module document block
!
! module:   nemsio_module      API for NEMS input/output 
!
! Abstract: This module handles NEMS input/output
!
! Program history log
!    2006-11-10    Jun Wang  for gfsio
!    2008-02-29    Jun Wang
!
! Public Variables
! Public Defined Types
!   nemsio_gfile
!     private
!        gtype:   character(nemsio_charkind*2)  NEMSIO file identifier
!        gdatatype:character(nemsio_charkind) data format
!        modelname:character(nemsio_charkind) modelname
!        version: integer(nemsio_intkind)   verion number
!        nmeta:   integer(nemsio_intkind)   number of metadata rec
!        lmeta:   integer(nemsio_intkind)   length of metadata rec 2 for model paramodels
!        nrec:    integer(nemsio_intkind)   number of data rec
!        idate(1:7):integer(nemsio_intkind) initial date (yyyy/mm/dd/hh/mm/ssn/ssd)
!        nfday:   integer(nemsio_intkind)   forecast day
!        nfhour:  integer(nemsio_intkind)   forecast hour
!        nfminute:integer(nemsio_intkind)   forecast minutes
!        nfsecondn:integer(nemsio_intkind)  numerator of forecast second fraction
!        nfsecondd:integer(nemsio_intkind)  denominator of forecast second fraction
!        dimy:    integer(nemsio_intkind)   dimension in latitude
!        dimx:    integer(nemsio_intkind)   dimension in Longitude
!        dimz:    integer(nemsio_intkind)   number of levels
!        nframe:  integer(nemsio_intkind)   dimension of halo
!        nsoil:    integer(nemsio_intkind)  number of soil layers
!        ntrac:    integer(nemsio_intkind)  number of tracers
!        jcap:    integer(nemsio_intkind)   spectral truncation
!        ncldt:   integer(nemsio_intkind)   number of cloud types
!        idsl:    integer(nemsio_intkind)   semi-lagrangian id
!        idvc:    integer(nemsio_intkind)   vertical coordinate id
!        idvm:    integer(nemsio_intkind)   mass variable id
!        idrt:    integer(nemsio_intkind)   grid identifier
!                 (idrt=4 for gaussian grid,
!                  idrt=0 for equally-spaced grid including poles,
!                  idrt=256 for equally-spaced grid excluding poles)
!        rlon_min:real(nemsio_realkind)     minimal longtitude of regional domain (global:set to 0)
!        rlon_max:real(nemsio_realkind)     maximal longtitude of regional domain (global:set to 360.)
!        rlat_min:real(nemsio_realkind)     minimal longtitude of regional domain (global:set to -90)
!        rlat_max:real(nemsio_realkind)     maximal longtitude of regional domain (global:set to 90)
!        extrameta:logical(nemsio_logickind)extra meta data flag 
!        nmetavari:integer(nemsio_intkind)  number of extra meta data integer variables
!        nmetavarr:integer(nemsio_intkind)  number of extra meta data real variables
!        nmetavarl:integer(nemsio_intkind)  number of extra meta data logical variables
!        nmetavarc:integer(nemsio_intkind)  number of extra meta data character variables
!        nmetaaryi:integer(nemsio_intkind)  number of extra meta data integer arrays
!        nmetaaryr:integer(nemsio_intkind)  number of extra meta data real arrays
!        nmetaaryl:integer(nemsio_intkind)  number of extra meta data logical arrays
!        nmetaaryc:integer(nemsio_intkind)  number of extra meta data character arrays
!
!        recname: character(nemsio_charkind),allocatable    recname(:)
!        reclevtyp: character(nemsio_charkind*2),allocatable    reclevtyp(:)
!        reclev:  integer(nemsio_intkind),allocatable       reclev(:)
!        vcoord:  real(nemsio_realkind),allocatable         vcoord(:,:,:)
!        lat:  real(nemsio_realkind),allocatable         lat(:) lat for mess point
!        lon:  real(nemsio_realkind),allocatable         lon(:) lon for mess point
!        gvlat1d: real(nemsio_realkind),allocatable         gvlat1d(:) lat for wind point
!        gvlon1d: real(nemsio_realkind),allocatable         gvlon1d(:) lon for wind point
!        Cpi:     real(nemsio_realkind),allocatable         cpi(:)
!        Ri:      real(nemsio_realkind),allocatable         ri(:)
!
!        variname:character(nemsio_charkind)  names of extra meta data integer variables
!        varrname:character(nemsio_charkind)  names of extra meta data real variables
!        varlname:character(nemsio_charkind)  names of extra meta data logical variables
!        varcname:character(nemsio_charkind)  names of extra meta data character variables
!        varival: integer(nemsio_intkind)     values of extra meta data integer variables
!        varrval: real(nemsio_realkind)       values of extra meta data integer variables
!        varlval: logical(nemsio_logickind)   values of extra meta data integer variables
!        varcval: character(nemsio_charkind)  values of extra meta data integer variables
!        aryiname:character(nemsio_charkind)  names of extra meta data integer arrays
!        aryrname:character(nemsio_charkind)  names of extra meta data real arrays
!        arylname:character(nemsio_charkind)  names of extra meta data logical arrays
!        arycname:character(nemsio_charkind)  names of extra meta data character arrays
!        aryilen: integer(nemsio_intkind)     lengths of extra meta data integer arrays
!        aryilen: integer(nemsio_intkind)     number of extra meta data integer arrays
!        aryilen: integer(nemsio_intkind)     number of extra meta data integer arrays

!!--- file handler
!        gfname:  character(255)  file name
!        gaction: character(nemsio_charkind)  read/write
!        flunit:  integer(nemsio_intkind)  unit number  
!
! Public method
!   nemsio_init
!   nemsio_finalize
!   nemsio_open
!   nemsio_writerec
!   nemsio_readirec
!   nemsio_writerecv
!   nemsio_readirecv
!   nemsio_writerecw34
!   nemsio_readirecw34
!   nemsio_writerecvw34
!   nemsio_readirecvw34
!   nemsio_close
!   nemsio_getfilehead
! Possible return code
!          0   Successful call
!         -1   Open or close I/O error
!         -2   array size
!         -3   Meta data I/O error (possible EOF)
!         -4   GETGB/PUTGB error
!         -5   Search record and set GRIB message info error
!         -6   allocate/deallocate error
!         -7   set grib table
!         -8   file meta data initialization (default:1152*576)
!         -9   NOT nemsio type file
!         -10  get/close file unit
!         -11  read/write bin data
!         -12  read/write NMM B grid lat lon
!         -13  read/write NMM sfc var
!         -15  read/write gsi 
!         -17  get var from file header
!
!$$$ end module document block
!
  use mpi
!
  implicit none
  private
!------------------------------------------------------------------------------
! private variables and type needed by nemsio_gfile
  integer,parameter:: nemsio_lmeta1=48,nemsio_lmeta3=32
  integer,parameter:: nemsio_intkind=4,nemsio_intkind8=8
  integer,parameter:: nemsio_realkind=4,nemsio_dblekind=8
  integer,parameter:: nemsio_charkind=16,nemsio_charkind8=8, nemsio_charkind4=4
  integer,parameter:: nemsio_logickind=4
  integer,parameter:: nemsio_maxint=2147483647
  real(nemsio_intkind),parameter     :: nemsio_intfill=-9999_nemsio_intkind
  logical(nemsio_logickind),parameter:: nemsio_logicfill=.false.
  real(nemsio_intkind),parameter     :: nemsio_kpds_intfill=-1_nemsio_intkind
  real(nemsio_realkind),parameter    :: nemsio_realfill=-9999._nemsio_realkind
  real(nemsio_dblekind),parameter    :: nemsio_dblefill=-9999._nemsio_dblekind
!
!------------------------------------------------------------------------------
!---  public types
  type,public :: nemsio_gfile
    private
    character(nemsio_charkind8) :: gtype=' '
    integer(nemsio_intkind):: version=nemsio_intfill
    character(nemsio_charkind8):: gdatatype=' '
    character(nemsio_charkind8):: modelname=' '
    integer(nemsio_intkind):: nmeta=nemsio_intfill
    integer(nemsio_intkind):: lmeta=nemsio_intfill
    integer(nemsio_intkind):: nrec=nemsio_intfill
!
    integer(nemsio_intkind):: idate(7)=nemsio_intfill
    integer(nemsio_intkind):: nfday=nemsio_intfill
    integer(nemsio_intkind):: nfhour=nemsio_intfill
    integer(nemsio_intkind):: nfminute=nemsio_intfill
    integer(nemsio_intkind):: nfsecondn=nemsio_intfill
    integer(nemsio_intkind):: nfsecondd=nemsio_intfill
!    integer(nemsio_intkind):: ifdate(7)=nemsio_intfill
!
    integer(nemsio_intkind):: dimx=nemsio_intfill
    integer(nemsio_intkind):: dimy=nemsio_intfill
    integer(nemsio_intkind):: dimz=nemsio_intfill
    integer(nemsio_intkind):: nframe=nemsio_intfill
    integer(nemsio_intkind):: nsoil=nemsio_intfill
    integer(nemsio_intkind):: ntrac=nemsio_intfill
!
    integer(nemsio_intkind) :: jcap=nemsio_intfill
    integer(nemsio_intkind) :: ncldt=nemsio_intfill
    integer(nemsio_intkind) :: idvc=nemsio_intfill
    integer(nemsio_intkind) :: idsl=nemsio_intfill
    integer(nemsio_intkind) :: idvm=nemsio_intfill
    integer(nemsio_intkind) :: idrt=nemsio_intfill
    real(nemsio_realkind) :: rlon_min=nemsio_realfill
    real(nemsio_realkind) :: rlon_max=nemsio_realfill
    real(nemsio_realkind) :: rlat_min=nemsio_realfill
    real(nemsio_realkind) :: rlat_max=nemsio_realfill
    logical(nemsio_logickind) :: extrameta=nemsio_logicfill
!
    integer(nemsio_intkind):: nmetavari=nemsio_intfill
    integer(nemsio_intkind):: nmetavarr=nemsio_intfill
    integer(nemsio_intkind):: nmetavarl=nemsio_intfill
    integer(nemsio_intkind):: nmetavarc=nemsio_intfill
    integer(nemsio_intkind):: nmetaaryi=nemsio_intfill
    integer(nemsio_intkind):: nmetaaryr=nemsio_intfill
    integer(nemsio_intkind):: nmetaaryl=nemsio_intfill
    integer(nemsio_intkind):: nmetaaryc=nemsio_intfill
!
    character(nemsio_charkind),allocatable :: recname(:)
    character(nemsio_charkind),allocatable :: reclevtyp(:)
    integer(nemsio_intkind),allocatable    :: reclev(:)
!
    real(nemsio_realkind),allocatable      :: vcoord(:,:,:)
    real(nemsio_realkind),allocatable      :: lat(:)
    real(nemsio_realkind),allocatable      :: lon(:)
    real(nemsio_realkind),allocatable      :: dx(:)
    real(nemsio_realkind),allocatable      :: dy(:)
!
    real(nemsio_realkind),allocatable      :: Cpi(:)
    real(nemsio_realkind),allocatable      :: Ri(:)
!
    character(nemsio_charkind),allocatable :: variname(:)
    integer(nemsio_intkind),allocatable    :: varival(:)
    character(nemsio_charkind),allocatable :: varrname(:)
    real(nemsio_realkind),allocatable      :: varrval(:)
    character(nemsio_charkind),allocatable :: varlname(:)
    logical(nemsio_logickind),allocatable  :: varlval(:)
    character(nemsio_charkind),allocatable :: varcname(:)
    character(nemsio_charkind),allocatable :: varcval(:)
!
    character(nemsio_charkind),allocatable :: aryiname(:)
    integer(nemsio_intkind),allocatable    :: aryilen(:)
    integer(nemsio_intkind),allocatable    :: aryival(:,:)
    character(nemsio_charkind),allocatable :: aryrname(:)
    integer(nemsio_intkind),allocatable    :: aryrlen(:)
    real(nemsio_realkind),allocatable      :: aryrval(:,:)
    character(nemsio_charkind),allocatable :: arylname(:)
    integer(nemsio_intkind),allocatable    :: aryllen(:)
    logical(nemsio_logickind),allocatable  :: arylval(:,:)
    character(nemsio_charkind),allocatable :: arycname(:)
    integer(nemsio_intkind),allocatable    :: aryclen(:)
    character(nemsio_charkind),allocatable :: arycval(:,:)
!  
    character(255) :: gfname
    character(nemsio_charkind8) :: gaction
    integer(nemsio_intkind8)    :: tlmeta=nemsio_intfill
    integer(nemsio_intkind)    :: fieldsize=nemsio_intfill
    integer(nemsio_intkind)    :: flunit=nemsio_intfill
    integer(nemsio_intkind)    :: headvarinum=nemsio_intfill
    integer(nemsio_intkind)    :: headvarrnum=nemsio_intfill
    integer(nemsio_intkind)    :: headvarcnum=nemsio_intfill
    integer(nemsio_intkind)    :: headvarlnum=nemsio_intfill
    integer(nemsio_intkind)    :: headaryinum=nemsio_intfill
    integer(nemsio_intkind)    :: headaryrnum=nemsio_intfill
    integer(nemsio_intkind)    :: headarycnum=nemsio_intfill
    character(nemsio_charkind),allocatable :: headvarcname(:)
    character(nemsio_charkind),allocatable :: headvariname(:)
    character(nemsio_charkind),allocatable :: headvarrname(:)
    character(nemsio_charkind),allocatable :: headvarlname(:)
    character(nemsio_charkind),allocatable :: headaryiname(:)
    character(nemsio_charkind),allocatable :: headaryrname(:)
    character(nemsio_charkind),allocatable :: headarycname(:)
    integer(nemsio_intkind),allocatable    :: headvarival(:)
    real(nemsio_realkind),allocatable      :: headvarrval(:)
    character(nemsio_charkind),allocatable :: headvarcval(:)
    logical(nemsio_logickind),allocatable  :: headvarlval(:)
    integer(nemsio_intkind),allocatable    :: headaryival(:,:)
    real(nemsio_realkind),allocatable      :: headaryrval(:,:)
    character(nemsio_charkind),allocatable :: headarycval(:,:)
    character,allocatable      :: cbuf(:)
    integer(nemsio_intkind):: mbuf=0,nlen,nnum,mnum
!-- for MPI I/O
    integer(nemsio_intkind)     :: mpi_comm=nemsio_intfill
    integer(nemsio_intkind)     :: lead_task=nemsio_intfill
    integer(nemsio_intkind)     :: mype=nemsio_intfill
    integer(nemsio_intkind)     :: npes=nemsio_intfill
    integer(nemsio_intkind)     :: fh=nemsio_intfill
    real(nemsio_realkind)       :: fieldsize_real4=nemsio_realfill
    real(nemsio_dblekind)       :: fieldsize_real8=nemsio_realfill
  end type nemsio_gfile
!
!------------------------------------------------------------------------------
!--- private types
!
  type :: nemsio_meta1
    sequence
     character(nemsio_charkind8) :: gtype
     character(nemsio_charkind8) :: modelname
     character(nemsio_charkind8) :: gdatatype
     integer(nemsio_intkind) :: version,nmeta,lmeta
     integer(nemsio_intkind) :: reserve(3)
  end type nemsio_meta1
!
  type :: nemsio_meta2
    sequence
    integer(nemsio_intkind) :: nrec 
    integer(nemsio_intkind) :: idate(1:7),nfday,nfhour,nfminute,nfsecondn, &
                               nfsecondd,dimx,dimy,dimz,nframe,nsoil,ntrac,&
                               jcap,ncldt,idvc,idsl,idvm,idrt
    real(nemsio_realkind)   :: rlon_min,rlon_max,rlat_min,rlat_max 
    logical(nemsio_logickind) :: extrameta
  end type nemsio_meta2
!
  type :: nemsio_meta3
    integer(nemsio_intkind) :: nmetavari,nmetavarr,nmetavarl,nmetavarc, &
                               nmetaaryi,nmetaaryr,nmetaaryl,nmetaaryc
  end type nemsio_meta3
!
!*** for mpi
  integer(nemsio_intkind)   :: itypemeta1,itypemeta2,itypemeta3
!
  type  :: nemsio_grbmeta
    integer(nemsio_intkind)   :: jf=nemsio_intfill
    integer(nemsio_intkind)   :: j=nemsio_kpds_intfill
    logical*1,allocatable     :: lbms(:)
    integer(nemsio_intkind)   :: jpds(200)=nemsio_kpds_intfill
    integer(nemsio_intkind)   :: jgds(200)=nemsio_kpds_intfill
  end type nemsio_grbmeta
!
!----- interface
  interface nemsio_getheadvar
    module procedure nemsio_getfheadvari
    module procedure nemsio_getfheadvarr
    module procedure nemsio_getfheadvarl
    module procedure nemsio_getfheadvarc
    module procedure nemsio_getfheadaryi
    module procedure nemsio_getfheadaryr
    module procedure nemsio_getfheadaryl
    module procedure nemsio_getfheadaryc
  end interface nemsio_getheadvar
!
  interface nemsio_readrec
    module procedure nemsio_readrecbin4
    module procedure nemsio_readrecbin8
  end interface nemsio_readrec
!
  interface nemsio_readrecv
    module procedure nemsio_readrecvbin4
    module procedure nemsio_readrecvbin8
  end interface nemsio_readrecv
!
  interface nemsio_writerec
    module procedure nemsio_writerecbin4
    module procedure nemsio_writerecbin8
  end interface nemsio_writerec
!
  interface nemsio_writerecv
    module procedure nemsio_writerecvbin4
    module procedure nemsio_writerecvbin8
  end interface nemsio_writerecv
!
  interface nemsio_denseread
    module procedure nemsio_denseread4
    module procedure nemsio_denseread8
  end interface nemsio_denseread
!
  interface nemsio_densewrite
    module procedure nemsio_densewrite4
    module procedure nemsio_densewrite8
  end interface nemsio_densewrite
!
!--- file unit for putgb/getgb ----
  integer(nemsio_intkind),save   :: fileunit(600:699)=0
!------------------------------------------------------------------------------
!public mehtods
  public nemsio_intkind,nemsio_intkind8,nemsio_realkind,nemsio_dblekind
  public nemsio_charkind,nemsio_charkind8,nemsio_logickind
  public nemsio_init,nemsio_finalize,nemsio_open,nemsio_close
  public nemsio_readrec,nemsio_writerec,nemsio_readrecv,nemsio_writerecv
  public nemsio_denseread,nemsio_densewrite
  public nemsio_getfilehead,nemsio_getheadvar,nemsio_getrechead
!
contains
!-------------------------------------------------------------------------------
  subroutine nemsio_init(iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! initialization
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    implicit none
    integer(nemsio_intkind),optional,intent(out):: iret
!-- local vars
    integer :: meta1_type(2),meta1_block(2),meta1_disp(2)
    integer :: meta2_type(3),meta2_block(3),meta2_disp(3)
    integer :: ios
!------------------------------------------------------------
! MPI set meta data type
!------------------------------------------------------------
! 1. meta1
    meta1_type(1)=MPI_CHARACTER
    meta1_type(2)=MPI_INTEGER
    meta1_block(1)=24
    meta1_block(2)=6
    meta1_disp(1)=0
    meta1_disp(2)=meta1_disp(1)+meta1_block(1)*1
    call mpi_type_struct(2,meta1_block,meta1_disp,meta1_type,            &
           itypemeta1,ios)
    call mpi_type_commit(itypemeta1,ios)
    if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop('nemsio, stop at mpi_type_struct for meta1')
       endif
    endif
!
! 2. meta2
    meta2_type(1)=MPI_INTEGER
    meta2_type(2)=MPI_REAL
    meta2_type(3)=MPI_LOGICAL
    meta2_block(1)=25
    meta2_block(2)=4
    meta2_block(3)=1
    meta2_disp(1)=0
    meta2_disp(2)=meta2_block(1)*4+meta2_disp(1)
    meta2_disp(3)=meta2_block(2)*4+meta2_disp(2)
    call mpi_type_struct(3,meta2_block,meta2_disp,meta2_type,            &
         itypemeta2,ios)
    call mpi_type_commit(itypemeta2,ios)
    if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop('nemsio, stop at mpi_type_struct for meta2')
       endif
    endif
!
! 3. meta3
    call mpi_type_contiguous(8,MPI_INTEGER,itypemeta3,ios)
    call mpi_type_commit(itypemeta3,ios)
    if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop('nemsio, stop at mpi_type_struct for meta3')
       endif
    endif
!
  end subroutine nemsio_init
!------------------------------------------------------------------------------
  subroutine nemsio_finalize()
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! abstract: finalization
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    implicit none
!--
  end subroutine nemsio_finalize
!------------------------------------------------------------------------------
  subroutine nemsio_open(gfile,gfname,gaction,mpi_comm,                     &
             iret,gdatatype,version,mype,npes,                              &
      nmeta,lmeta,modelname,nrec,idate,nfday,nfhour,nfminute,nfsecondn,     &
      nfsecondd, &
      dimx,dimy,dimz,nframe,nsoil,ntrac,jcap,ncldt,idvc,idsl,idvm,idrt,     &
      rlon_min,rlon_max,rlat_min,rlat_max,extrameta,           &
      nmetavari,nmetavarr,nmetavarl,nmetavarc,                              &
      nmetaaryi,nmetaaryr,nmetaaryl,nmetaaryc,                              &
      recname,reclevtyp,reclev,vcoord,lat,lon,dx,dy,cpi,ri,                 &
      variname,varival,varrname,varrval,varlname,varlval,varcname,varcval,  &
      aryiname,aryilen,aryival,aryrname,aryrlen,aryrval,                    &
      arylname,aryllen,arylval,arycname,aryclen,arycval  )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! abstract: open nemsio file, and read/write the meta data
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)    :: gfile
    character*(*),intent(in)            :: gfname
    character*(*),intent(in)            :: gaction
    integer,intent(in)                  :: mpi_comm
!-------------------------------------------------------------------------------
! optional variables
!-------------------------------------------------------------------------------
    integer(nemsio_intkind),optional,intent(out) :: iret
    character*(*),optional,intent(in)            :: gdatatype,modelname
    integer(nemsio_intkind),optional,intent(in)  :: version,nmeta,lmeta,nrec
    integer,optional,intent(in)                  :: mype,npes
    integer(nemsio_intkind),optional,intent(in)  :: idate(7),nfday,nfhour,    &
            nfminute, nfsecondn,nfsecondd
    integer(nemsio_logickind),optional,intent(in):: dimx,dimy,dimz,nframe,    &
            nsoil,ntrac
    integer(nemsio_logickind),optional,intent(in):: jcap,ncldt,idvc,idsl,     &
            idvm,idrt
    real(nemsio_realkind),optional,intent(in)    :: rlat_min,rlat_max,   &
             rlon_min,rlon_max
    logical(nemsio_logickind),optional,intent(in):: extrameta
    integer(nemsio_intkind),optional,intent(in)  :: nmetavari,nmetavarr, &   
            nmetavarl,nmetavarc,nmetaaryi,nmetaaryr,nmetaaryl,nmetaaryc
!
    character*(*),optional,intent(in)            :: recname(:),reclevtyp(:)
    integer(nemsio_intkind),optional,intent(in)  :: reclev(:)
    real(nemsio_realkind),optional,intent(in)    :: vcoord(:,:,:)
    real(nemsio_realkind),optional,intent(in)    :: lat(:),lon(:)
    real(nemsio_realkind),optional,intent(in)    :: dx(:),dy(:)
    real(nemsio_realkind),optional,intent(in)    :: Cpi(:),Ri(:)
!
    character*(*),optional,intent(in)            :: variname(:),varrname(:),&
          varlname(:),varcname(:),aryiname(:),aryrname(:),arylname(:),arycname(:)
    integer(nemsio_intkind),optional,intent(in)  :: aryilen(:),aryrlen(:),  &
          aryllen(:),aryclen(:)
    integer(nemsio_intkind),optional,intent(in)  :: varival(:),aryival(:,:)
    real(nemsio_realkind),optional,intent(in)    :: varrval(:),aryrval(:,:)
    logical(nemsio_logickind),optional,intent(in):: varlval(:),arylval(:,:)
    character(*),optional,intent(in)             :: varcval(:),arycval(:,:)
!
    integer :: ios
!------------------------------------------------------------
!### for MPI IO, just need this part for read header ###### 
!    assign a unit number 
!------------------------------------------------------------
    if (present(iret)) iret=-1
!
    gfile%gfname=gfname
    gfile%gaction=gaction
    gfile%mpi_comm=mpi_comm
    gfile%lead_task=0
!
    call nemsio_getlu(gfile,ios)
    if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop
       endif
    endif
    if(present(mype)) then
      gfile%mype=mype
    else
      call mpi_comm_rank(mpi_comm,gfile%mype,ios)
      if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop
       endif
      endif
    endif
    if(present(npes)) then
      gfile%mype=npes
    else
      call mpi_comm_size(mpi_comm,gfile%npes,ios)
      if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop
       endif
      endif
    endif
!
!------------------------------------------------------------
! open and read meta data for READ
!------------------------------------------------------------
    if ( gaction .eq. "read" .or. gaction .eq. "READ") then
!
!-read  meta data for gfile, use non-mpi read for header
!
       call nemsio_rcreate(gfile,ios)
!       write(0,*)'after nemsio_rcreate'
       if ( ios.ne.0) then
        if ( present(iret))  then
          iret=ios
          return
        else
          call nemsio_stop
        endif
       endif
!
!-read 2D field using MPI I/O
!
       call mpi_file_open(mpi_comm,gfname,MPI_MODE_RDONLY,MPI_INFO_NULL,gfile%fh,ios)
       if ( ios.ne.0) then
        if ( present(iret))  then
          return
        else
          call nemsio_stop
        endif
       endif
!------------------------------------------------------------
! open and write meta data for WRITE
!------------------------------------------------------------
    elseif (gaction .eq. "write" .or. gaction .eq. "WRITE") then
!
!-write  meta data for gfile, use non-mpi write for header
!
      call nemsio_wcreate(gfile,ios,gdatatype=gdatatype, &
        version=version, nmeta=nmeta,lmeta=lmeta,modelname=modelname,  &
        nrec=nrec,idate=idate,nfday=nfday,nfhour=nfhour,nfminute=nfminute,&
        nfsecondn=nfsecondn, nfsecondd=nfsecondd, &
        dimx=dimx,dimy=dimy,dimz=dimz,nframe=nframe,nsoil=nsoil,   &
        ntrac=ntrac,jcap=jcap,ncldt=ncldt,idvc=idvc,idsl=idsl,    &
        idvm=idvm,idrt=idrt,                          &
        rlon_min=rlon_min,rlon_max=rlon_max,rlat_min=rlat_min, &
        rlat_max=rlat_max,extrameta=extrameta, &
        nmetavari=nmetavari,nmetavarr=nmetavarr,  &
        nmetavarl=nmetavarl,nmetaaryi=nmetaaryi,nmetaaryr=nmetaaryr,&
        nmetaaryl=nmetaaryl,recname=recname,reclevtyp=reclevtyp,    &
        reclev=reclev,vcoord=vcoord,lat=lat,lon=lon,dx=dx,dy=dy,    &
        cpi=cpi,ri=ri,variname=variname,varival=varival,varrname=varrname,&
        varrval=varrval,varlname=varlname,varlval=varlval, &
        varcname=varcname,varcval=varcval, &
        aryiname=aryiname,aryilen=aryilen,aryival=aryival, &
        aryrname=aryrname,aryrlen=aryrlen,aryrval=aryrval, &
        arylname=arylname,aryllen=aryllen,arylval=arylval, &
        arycname=arycname,aryclen=aryclen,arycval=arycval  )
      if ( ios.ne.0) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop
       endif
     endif
!
!-write 2D field using MPI I/O
!
      call mpi_file_open(mpi_comm,gfname,MPI_MODE_CREATE+MPI_MODE_WRONLY,    &
           MPI_INFO_NULL,gfile%fh,ios)
      if ( ios.ne.0) then
       if ( present(iret))  then
         return
       else
         call nemsio_stop
       endif
      endif
!
!------------------------------------------------------------
! if gaction is wrong
!------------------------------------------------------------
    else
       if ( present(iret))  then
         return
       else
         call nemsio_stop
       endif
    endif
!------------------------------------------------------------
! set default header
!------------------------------------------------------------
    if(.not.allocated(gfile%headvariname).or. &
       .not.allocated(gfile%headvarrname).or. &
       .not.allocated(gfile%headvarcname).or. &
       .not.allocated(gfile%headvarlname).or. &
       .not.allocated(gfile%headaryiname).or. &
       .not.allocated(gfile%headaryrname) ) then
      call nemsio_setfhead(gfile,ios)
      if ( present(iret)) iret=ios
      if ( ios.ne.0) then
        if (present(iret)) return
        call nemsio_stop
      endif
    endif
!
    iret=0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine nemsio_open
!-------------------------------------------------------------------------------
  subroutine nemsio_close(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! abstract: close gfile including closing the file, returning unit number, 
!           setting file meta data empty
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)     :: gfile
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer(nemsio_intkind)      :: ios
!------------------------------------------------------------
! close the file
!------------------------------------------------------------
    if ( present(iret) ) iret=-1
    call mpi_file_close(gfile%fh,ios)
    if ( ios.ne.0) then
       if ( present(iret))  then
         return
       else
         call nemsio_stop
       endif
    endif
!------------------------------------------------------------
! empty gfile meta data
!------------------------------------------------------------
    call nemsio_axmeta(gfile,ios)
    if ( ios.ne.0) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop
       endif
    endif
    if ( present(iret)) iret=0
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  end subroutine nemsio_close
!------------------------------------------------------------------------------
  subroutine nemsio_rcreate(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio meta data
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)     :: gfile
    integer(nemsio_intkind),intent(out)  :: iret
!local variables
    integer(nemsio_intkind)      :: ios,nmeta,tlmeta4
    integer(nemsio_intkind8)     :: iskip,iread,nread
    type(nemsio_meta1)           :: meta1
    type(nemsio_meta2)           :: meta2
    type(nemsio_meta3)           :: meta3
    integer(nemsio_intkind) :: i
    character(nemsio_charkind8),allocatable :: char8var(:)
!------------------------------------------------------------
! open gfile for read header
!------------------------------------------------------------
    iret=-3
!    print *,'in rcreate ',gfile%gfname,gfile%flunit,gfile%mype,gfile%lead_task
    if(gfile%mype.eq.gfile%lead_task) then
      call baopenr(gfile%flunit,gfile%gfname,ios)
      if ( ios.ne.0) return
!------------------------------------------------------------
! read first meta data record
!------------------------------------------------------------
      iskip=0
      iread=nemsio_lmeta1
      call bafrreadl(gfile%flunit,iskip,iread,nread,meta1)
      if(nread.lt.iread) return
      gfile%tlmeta=nread
!    print *,'lead_task,after meta1,gtype=',meta1%gtype,meta1%gdatatype,  &
!      meta1%modelname,meta1%version,meta1%nmeta,meta1%lmeta
    endif
!
    call MPI_BCAST(meta1,1,itypemeta1,gfile%lead_task,     &
       gfile%mpi_comm,ios)
    gfile%gtype=meta1%gtype
    gfile%gdatatype=meta1%gdatatype
    gfile%modelname=meta1%modelname
    gfile%version=meta1%version
    gfile%nmeta=meta1%nmeta
    gfile%lmeta=meta1%lmeta
!    print *,'after meta1,gtype=',meta1%gtype,meta1%gdatatype,  &
!      meta1%modelname,meta1%version,meta1%nmeta,meta1%lmeta,ios
    if ( trim(gfile%gdatatype).ne."bin4" .and. trim(gfile%gdatatype).ne."bin8" &
         .and. trim(gfile%gdatatype).ne."grib" ) then
      gfile%gdatatype="grib"
    endif
    if ( gfile%gtype(1:6) .ne. 'NEMSIO' ) then
      iret=-9
      return
    endif
    if ( gfile%nmeta .ne. 12 ) then
      print*,'WARNING: Not standard meta data, may not be ingested into GSI!!!'
      iret=-9
      return
    endif
!------------------------------------------------------------
! read second meta data record
!------------------------------------------------------------
    if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=gfile%lmeta
      call bafrreadl(gfile%flunit,iskip,iread,nread,meta2)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'2 tlmeta=',gfile%tlmeta
    endif
!
    call MPI_BCAST(meta2,1,itypemeta2,gfile%lead_task,       &
       gfile%mpi_comm,ios)
    gfile%nrec=meta2%nrec
    gfile%idate(1:7)=meta2%idate(1:7)
    gfile%nfday=meta2%nfday
    gfile%nfhour=meta2%nfhour
    gfile%nfminute=meta2%nfminute
    gfile%nfsecondn=meta2%nfsecondn
    gfile%nfsecondd=meta2%nfsecondd
    gfile%dimx=meta2%dimx
    gfile%dimy=meta2%dimy
    gfile%dimz=meta2%dimz
    gfile%nframe=meta2%nframe
    gfile%nsoil=meta2%nsoil
    gfile%ntrac=meta2%ntrac
    gfile%jcap=meta2%jcap
    gfile%ncldt=meta2%ncldt
    gfile%idvc=meta2%idvc
    gfile%idsl=meta2%idsl
    gfile%idvm=meta2%idvm
    gfile%idrt=meta2%idrt
    gfile%rlon_min=meta2%rlon_min
    gfile%rlon_max=meta2%rlon_max
    gfile%rlat_min=meta2%rlat_min
    gfile%rlat_max=meta2%rlat_max
    gfile%extrameta=meta2%extrameta
    gfile%fieldsize=(gfile%dimx+2*gfile%nframe)*(gfile%dimy+2*gfile%nframe)
!    print *,'meta2,nrec=',gfile%nrec,gfile%idate(1:7),gfile%nfday,  &
!      gfile%nfhour,gfile%nfminute,gfile%nfsecondn,gfile%nfsecondd,  &
!      gfile%dimx,gfile%dimy,gfile%dimz,gfile%nframe,gfile%nsoil,    &
!      gfile%ntrac,gfile%jcap,gfile%ncldt,gfile%idvc,gfile%idsl,     &
!      gfile%idvm,gfile%idrt,gfile%rlon_min,gfile%rlon_max,          &
!      gfile%rlat_min,gfile%rlat_max,gfile%extrameta

    nmeta=gfile%nmeta-2
!------------------------------------------------------------
! set up gfile required meata arrays
!------------------------------------------------------------
    call nemsio_almeta(gfile,ios)
    if ( ios .ne. 0 ) then
      iret=ios
      return
    endif
!------------------------------------------------------------
! read gfile meta data array (meta rec 3:13)
!------------------------------------------------------------
!meta3:recname
    if(gfile%mype.eq.gfile%lead_task) then
      if ( nmeta.le.3 ) then
        print *,'WRONG: Please provide names,level type and  &
     &   levs for the fields in the nemsio file'
        return
      endif
      if(gfile%nmeta-2>0) then
        iskip=iskip+nread
        iread=len(gfile%recname)*size(gfile%recname)
        call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%recname)
        if(nread.lt.iread)  then
           iread=nemsio_charkind8*size(gfile%recname)
           allocate(char8var(size(gfile%recname)))
           call bafrreadl(gfile%flunit,iskip,iread,nread,char8var)
           gfile%recname=char8var
           deallocate(char8var)
           if (nread.lt.iread) return
        endif
        nmeta=nmeta-1
        gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetarecname =',gfile%tlmeta,'nread=',nread
     endif
    endif
    call MPI_BCAST(gfile%recname,gfile%nrec*nemsio_charkind,   &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
!
!meta4:reclevtyp
    if(gfile%mype.eq.gfile%lead_task) then
      if (gfile%nmeta-3>0 ) then
        iskip=iskip+nread
        iread=len(gfile%reclevtyp)*size(gfile%reclevtyp)
        call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%reclevtyp)
        if(nread.lt.iread) return
        nmeta=nmeta-1
        gfile%tlmeta=gfile%tlmeta+nread
!        print *,'tlmetareclwvtyp =',gfile%tlmeta,'nread=',nread
      endif
    endif
    call MPI_BCAST(gfile%reclevtyp,gfile%nrec*nemsio_charkind,   &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
!
!meta5:reclev
    if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=nemsio_intkind*size(gfile%reclev)
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%reclev)
      if(nread.lt.iread) return
      nmeta=nmeta-1
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetareclev =',gfile%tlmeta,'nread=',nread
    endif
    call MPI_BCAST(gfile%reclev,size(gfile%reclev),   &
       MPI_INTEGER,gfile%lead_task,gfile%mpi_comm,ios)
!meta6:vcoord
    if(gfile%mype.eq.gfile%lead_task) then
     if ( nmeta.gt.0 ) then
      iskip=iskip+nread
      iread=nemsio_realkind*size(gfile%vcoord)
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%vcoord)
      if(nread.lt.iread) return
      nmeta=nmeta-1
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetavcoord =',gfile%tlmeta,'nread=',nread
     endif
    endif
    call MPI_BCAST(gfile%vcoord,size(gfile%vcoord),   &
       MPI_REAL,gfile%lead_task,gfile%mpi_comm,ios)
!meta7:lat
    if(gfile%mype.eq.gfile%lead_task) then
     if ( nmeta.gt.0 ) then
      iskip=iskip+nread
      iread=nemsio_realkind*size(gfile%lat)
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%lat)
      if(nread.lt.iread) return
      nmeta=nmeta-1
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetareclat =',gfile%tlmeta,                  &
!         maxval(gfile%lat),minval(gfile%lat)
     endif
    endif
    call MPI_BCAST(gfile%lat,size(gfile%lat),   &
       MPI_REAL,gfile%lead_task,gfile%mpi_comm,ios)
!
!meta8:lon
    if(gfile%mype.eq.gfile%lead_task) then
     if ( nmeta.gt.0 ) then
      iskip=iskip+nread
      iread=nemsio_realkind*size(gfile%lon)
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%lon)
      if(nread.lt.iread) return
      nmeta=nmeta-1
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetareclon =',gfile%tlmeta,                 &
!         maxval(gfile%lon),minval(gfile%lon)
     endif
    endif
    call MPI_BCAST(gfile%lon,size(gfile%lon),   &
       MPI_REAL,gfile%lead_task,gfile%mpi_comm,ios)
!
!meta9:dx
    if(gfile%mype.eq.gfile%lead_task) then
     if ( nmeta.gt.0 ) then
      iskip=iskip+nread
      iread=nemsio_realkind*size(gfile%dx)
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%dx)
      if(nread.lt.iread) return
      nmeta=nmeta-1
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetarecdx =',gfile%tlmeta,                 &
!        maxval(gfile%dx),minval(gfile%dx)
     endif
    endif
    call MPI_BCAST(gfile%dx,size(gfile%dx),               &
       MPI_REAL,gfile%lead_task,gfile%mpi_comm,ios)
!meta10:dy
    if(gfile%mype.eq.gfile%lead_task) then
     if ( nmeta.gt.0 ) then
      iskip=iskip+nread
      iread=nemsio_realkind*size(gfile%dy)
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%dy)
      if(nread.lt.iread) return
      nmeta=nmeta-1
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetarecdy =',gfile%tlmeta,'nread=',nread, &
!         maxval(gfile%dy),maxval(gfile%dy)
     endif
    endif
    call MPI_BCAST(gfile%dy,size(gfile%dy),               &
       MPI_REAL,gfile%lead_task,gfile%mpi_comm,ios)
!meta11:cpi
    if(gfile%mype.eq.gfile%lead_task) then
     if ( nmeta .gt.0 ) then
      iskip=iskip+nread
      iread=nemsio_realkind*size(gfile%Cpi)
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%Cpi)
      if(nread.lt.iread) return
      nmeta=nmeta-1
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetacpi =',gfile%tlmeta,'cpi=',maxval(gfile%Cpi), &
!        minval(gfile%cpi)
     endif
    endif
    call MPI_BCAST(gfile%Cpi,size(gfile%Cpi),               &
       MPI_REAL,gfile%lead_task,gfile%mpi_comm,ios)
!Ri
    if(gfile%mype.eq.gfile%lead_task) then
     if ( nmeta.gt.0 ) then
      iskip=iskip+nread
      iread=nemsio_realkind*size(gfile%Ri)
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%Ri)
      if(nread.lt.iread) return
      nmeta=nmeta-1
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetri =',gfile%tlmeta,'ri=',maxval(gfile%ri), &
!        minval(gfile%ri)
     endif
    endif
    call MPI_BCAST(gfile%Ri,size(gfile%Ri),               &
       MPI_REAL,gfile%lead_task,gfile%mpi_comm,ios)
!
!    if ( nmeta.gt.0 ) then
!      print *,'nmeta=',nmeta,' WARNING:there are more meta to be read!'
!    endif
     
    if(gfile%extrameta) then
!------------------------------------------------------------
! read out extra meta data
!------------------------------------------------------------
    if(gfile%mype.eq.gfile%lead_task) then
     iskip=iskip+nread
     iread=nemsio_lmeta3
     call bafrreadl(gfile%flunit,iskip,iread,nread,meta3)
     if(nread.lt.iread) return
     gfile%tlmeta=gfile%tlmeta+nread
!     print *,'tlmetameta3 =',gfile%tlmeta,'nread=',nread
      
    endif
    call MPI_BCAST(meta3,1,itypemeta3,gfile%lead_task,        &
       gfile%mpi_comm,ios)
    gfile%nmetavari=meta3%nmetavari
    gfile%nmetavarr=meta3%nmetavarr
    gfile%nmetavarl=meta3%nmetavarl
    gfile%nmetavarc=meta3%nmetavarc
    gfile%nmetaaryi=meta3%nmetaaryi
    gfile%nmetaaryr=meta3%nmetaaryr
    gfile%nmetaaryl=meta3%nmetaaryl
    gfile%nmetaaryc=meta3%nmetaaryc
!      print *,'after meta3,nread=',nread, &
!     'nmetavari=',gfile%nmetavari,'nvarr=',gfile%nmetavarr, &
!     'varl=',gfile%nmetavarl,'varc=',gfile%nmetavarc, &
!     gfile%nmetaaryi,gfile%nmetaaryr,gfile%nmetaaryl,gfile%nmetaaryc
   
    call nemsio_alextrameta(gfile,ios)
    if ( ios .ne. 0 ) then
      iret=ios
      return
    endif

!meta var integer
    if (gfile%nmetavari.gt.0) then
     if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=nemsio_charkind*gfile%nmetavari
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%variname)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetavari =',gfile%tlmeta,'nread=',nread
!
      iskip=iskip+nread
      iread=nemsio_intkind*gfile%nmetavari
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%varival)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetavarival =',gfile%tlmeta,'nread=',nread
     endif
     call MPI_BCAST(gfile%variname,gfile%nmetavari*nemsio_charkind,  &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
     call MPI_BCAST(gfile%varival,gfile%nmetavari,               &
       MPI_INTEGER,gfile%lead_task,gfile%mpi_comm,ios)
!     print *,'in rcreate,after bcast ios=',ios, gfile%variname,gfile%varival
    endif
!
!meta var real
    if (gfile%nmetavarr.gt.0) then
     if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=nemsio_charkind*gfile%nmetavarr
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%varrname)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetavarr =',gfile%tlmeta,'nread=',nread

      iskip=iskip+nread
      iread=nemsio_realkind*gfile%nmetavarr
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%varrval)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetavarrval =',gfile%tlmeta,'nread=',nread
     endif
     call MPI_BCAST(gfile%varrname,gfile%nmetavarr*nemsio_charkind,   &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
     call MPI_BCAST(gfile%varrval,gfile%nmetavarr,               &
       MPI_REAL,gfile%lead_task,gfile%mpi_comm,ios)
    endif
!
!meta var logical
    if (gfile%nmetavarl.gt.0) then
     if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=nemsio_charkind*gfile%nmetavarl
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%varlname)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetavarl =',gfile%tlmeta,'nread=',nread

      iskip=iskip+nread
      iread=nemsio_logickind*gfile%nmetavarl
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%varlval)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetavarlval =',gfile%tlmeta,'nread=',nread
     endif
     call MPI_BCAST(gfile%varlname,gfile%nmetavarl*nemsio_charkind,   &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
     call MPI_BCAST(gfile%varlval,gfile%nmetavarl,               &
       MPI_LOGICAL,gfile%lead_task,gfile%mpi_comm,ios)
    endif
!
!meta var character
    if (gfile%nmetavarc.gt.0) then
     if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=nemsio_charkind*gfile%nmetavarc
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%varcname)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetavarc =',gfile%tlmeta,'nread=',nread
!
      iskip=iskip+nread
      iread=nemsio_charkind*gfile%nmetavarc
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%varcval)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetavarcval =',gfile%tlmeta,'nread=',nread
     endif
     call MPI_BCAST(gfile%varcname,gfile%nmetavarc*nemsio_charkind,  &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
     call MPI_BCAST(gfile%varcval,gfile%nmetavarc*nemsio_charkind,   &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
    endif
!
!meta arr integer
    if (gfile%nmetaaryi.gt.0) then
     if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=nemsio_charkind*gfile%nmetaaryi
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%aryiname)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetaaryinam =',gfile%tlmeta,'nread=',nread
!
      iskip=iskip+nread
      iread=nemsio_intkind*gfile%nmetaaryi
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%aryilen)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetaaryilen =',gfile%tlmeta,'nread=',nread
     endif
     call MPI_BCAST(gfile%aryiname,gfile%nmetaaryi*nemsio_charkind,  &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
     call MPI_BCAST(gfile%aryilen,gfile%nmetaaryi,   &
       MPI_INTEGER,gfile%lead_task,gfile%mpi_comm,ios)
     allocate(gfile%aryival(maxval(gfile%aryilen),gfile%nmetaaryi))
     do i=1,gfile%nmetaaryi
      if(gfile%mype.eq.gfile%lead_task) then
        iskip=iskip+nread
        iread=nemsio_intkind*gfile%aryilen(i)
        call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%aryival(:,i))
        if(nread.lt.iread) return
        gfile%tlmeta=gfile%tlmeta+nread
!        print *,'tlmetaaryival =',gfile%tlmeta,'nread=',nread
      endif
      call MPI_BCAST(gfile%aryival(:,i),gfile%aryilen(i),   &
       MPI_INTEGER,gfile%lead_task,gfile%mpi_comm,ios)
     enddo
    endif
!meta arr real
    if (gfile%nmetaaryr.gt.0) then
     if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=nemsio_charkind*gfile%nmetaaryr
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%aryrname)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetaaryrnam =',gfile%tlmeta,'nread=',nread

      iskip=iskip+nread
      iread=nemsio_intkind*gfile%nmetaaryr
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%aryrlen)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!      print *,'tlmetaaryrlen =',gfile%tlmeta,'nread=',nread,'nmetaaryr=',gfile%nmetaaryr
     endif
     call MPI_BCAST(gfile%aryrname,gfile%nmetaaryr*nemsio_charkind,  &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
     call MPI_BCAST(gfile%aryrlen,gfile%nmetaaryr,   &
       MPI_INTEGER,gfile%lead_task,gfile%mpi_comm,ios)
     allocate(gfile%aryrval(maxval(gfile%aryrlen),gfile%nmetaaryr) )
     do i=1,gfile%nmetaaryr
       if(gfile%mype.eq.gfile%lead_task) then
        iskip=iskip+nread
        iread=nemsio_realkind*gfile%aryrlen(i)
        call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%aryrval(:,i))
        if(nread.lt.iread) return
        gfile%tlmeta=gfile%tlmeta+nread
!        print *,'tlmetaaryrval =',gfile%tlmeta,'nread=',nread
       endif
       call MPI_BCAST(gfile%aryrval(:,i),gfile%aryrlen(i),   &
        MPI_REAL,gfile%lead_task,gfile%mpi_comm,ios)
     enddo
    endif
!meta arr logical
    if (gfile%nmetaaryl.gt.0) then
     if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=nemsio_charkind*gfile%nmetaaryl
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%arylname)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!
      iskip=iskip+nread
      iread=nemsio_intkind*gfile%nmetaaryl
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%aryllen)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
     endif
     call MPI_BCAST(gfile%arylname,gfile%nmetaaryl*nemsio_charkind,  &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
     call MPI_BCAST(gfile%aryllen,gfile%nmetaaryl,   &
       MPI_INTEGER,gfile%lead_task,gfile%mpi_comm,ios)
     allocate(gfile%arylval(maxval(gfile%aryllen),gfile%nmetaaryl) )
     do i=1,gfile%nmetaaryl
       if(gfile%mype.eq.gfile%lead_task) then
        iskip=iskip+nread
        iread=nemsio_logickind*gfile%aryllen(i)
        call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%arylval(:,i))
        if(nread.lt.iread) return
        gfile%tlmeta=gfile%tlmeta+nread
       endif
       call MPI_BCAST(gfile%arylval(:,i),gfile%aryllen(i),   &
        MPI_LOGICAL,gfile%lead_task,gfile%mpi_comm,ios)
     enddo
    endif
!meta arr char
    if (gfile%nmetaaryc.gt.0) then
     if(gfile%mype.eq.gfile%lead_task) then
      iskip=iskip+nread
      iread=nemsio_charkind*gfile%nmetaaryc
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%arycname)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
!
      iskip=iskip+nread
      iread=nemsio_intkind*gfile%nmetaaryc
      call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%aryclen)
      if(nread.lt.iread) return
      gfile%tlmeta=gfile%tlmeta+nread
     endif
     call MPI_BCAST(gfile%arycname,gfile%nmetaaryc*nemsio_charkind,  &
       MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
     call MPI_BCAST(gfile%aryclen,gfile%nmetaaryc,   &
       MPI_INTEGER,gfile%lead_task,gfile%mpi_comm,ios)
     allocate(gfile%arycval(maxval(gfile%aryclen),gfile%nmetaaryc) )
     do i=1,gfile%nmetaaryc
       if(gfile%mype.eq.gfile%lead_task) then
        iskip=iskip+nread
        iread=2*nemsio_charkind*gfile%aryclen(i)
        call bafrreadl(gfile%flunit,iskip,iread,nread,gfile%arycval(:,i))
        if(nread.lt.iread) return
        gfile%tlmeta=gfile%tlmeta+nread
       endif
       call MPI_BCAST(gfile%arycval(:,i),gfile%aryclen(i)*nemsio_charkind,   &
        MPI_CHARACTER,gfile%lead_task,gfile%mpi_comm,ios)
     enddo
    endif
!
!end if extrameta
   endif
!
!bcst tlmeta
    tlmeta4=gfile%tlmeta
    call MPI_BCAST(tlmeta4,1,MPI_INTEGER,   &
       gfile%lead_task,gfile%mpi_comm,ios)
    gfile%tlmeta=tlmeta4
!------------------------------------------------------------
! close the file
!------------------------------------------------------------
    if(gfile%mype.eq.gfile%lead_task) then
      call baclose(gfile%flunit,ios)
      if ( ios.ne.0) return
    endif
!
    call MPI_Barrier(gfile%mpi_comm, ios)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
   iret=0
  end subroutine nemsio_rcreate
!------------------------------------------------------------------------------
  subroutine nemsio_wcreate(gfile,iret,gdatatype,version,  &
      nmeta,lmeta,modelname,nrec,idate,nfday,nfhour,nfminute,nfsecondn,      &
      nfsecondd, &
      dimx,dimy,dimz,nframe,nsoil,ntrac,jcap,ncldt,idvc,idsl,idvm,idrt,     &
      rlon_min,rlon_max,rlat_min,rlat_max,extrameta,                        &
      nmetavari,nmetavarr,nmetavarl,nmetavarc,                              &
      nmetaaryi,nmetaaryr,nmetaaryl,nmetaaryc,                              &
      recname,reclevtyp,reclev,vcoord,lat,lon,dx,dy,cpi,ri,                 &
      variname,varival,varrname,varrval,varlname,varlval,varcname,varcval,  &
      aryiname,aryilen,aryival,aryrname,aryrlen,aryrval,                    &
      arylname,aryllen,arylval,arycname,aryclen,arycval  )
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: write nemsio meta data
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)             :: gfile
    integer(nemsio_intkind),intent(out)          :: iret
!optional variables
    character*(*),optional,intent(in)            :: gdatatype,modelname
    integer(nemsio_intkind),optional,intent(in)  :: version,nmeta,lmeta,nrec
    integer(nemsio_intkind),optional,intent(in)  :: idate(7),nfday,nfhour,  &
            nfminute,nfsecondn,nfsecondd
    integer(nemsio_logickind),optional,intent(in):: dimx,dimy,dimz,nframe,    &
            nsoil,ntrac
    integer(nemsio_logickind),optional,intent(in):: jcap,ncldt,idvc,idsl,     &
            idvm,idrt
    real(nemsio_realkind),optional,intent(in)    :: rlat_min,rlat_max,   &
             rlon_min,rlon_max
    logical(nemsio_logickind),optional,intent(in):: extrameta
    integer(nemsio_intkind),optional,intent(in)  :: nmetavari,nmetavarr, &
            nmetavarl,nmetavarc,nmetaaryi,nmetaaryr,nmetaaryl,nmetaaryc
!
    character*(*),optional,intent(in)            :: recname(:),reclevtyp(:)
    integer(nemsio_intkind),optional,intent(in)  :: reclev(:)
    real(nemsio_realkind),optional,intent(in)    :: vcoord(:,:,:)
    real(nemsio_realkind),optional,intent(in)    :: lat(:),lon(:)
    real(nemsio_realkind),optional,intent(in)    :: dx(:),dy(:)
    real(nemsio_realkind),optional,intent(in)    :: Cpi(:),Ri(:)
!
    character*(*),optional,intent(in)            :: variname(:),varrname(:),&
          varlname(:),varcname(:),aryiname(:),aryrname(:),arylname(:),arycname(:)
    integer(nemsio_intkind),optional,intent(in)  :: aryilen(:),aryrlen(:),  &
          aryllen(:),aryclen(:)
    integer(nemsio_intkind),optional,intent(in)  :: varival(:),aryival(:,:)
    real(nemsio_realkind),optional,intent(in)    :: varrval(:),aryrval(:,:)
    logical(nemsio_logickind),optional,intent(in):: varlval(:),arylval(:,:)
    character(*),optional,intent(in)             :: varcval(:),arycval(:,:)
!
!---  local variables
!
    real(nemsio_realkind) :: radi
    integer(nemsio_intkind8) :: iskip,iwrite,nwrite
    type(nemsio_meta1)      :: meta1
    type(nemsio_meta2)      :: meta2
    type(nemsio_meta3)      :: meta3
    integer(nemsio_intkind) :: i,n,ios,nummeta
    integer :: status(MPI_STATUS_SIZE) 
    logical :: linit
!------------------------------------------------------------
! set gfile meta data to operational model (default) if it's empty
!------------------------------------------------------------
    iret=-3
    gfile%gtype="NEMSIO"
    if(present(gdatatype)) then
      if ( trim(gdatatype).ne.'grib'.and.trim(gdatatype).ne.'bin4'.and. &
           trim(gdatatype).ne.'bin8' ) return
      gfile%gdatatype=gdatatype
    else
      gfile%gdatatype='grib'
    endif
    if(present(modelname)) then 
      gfile%modelname=modelname
    else
      gfile%modelname="GFS"
    endif
!
!    print *,'NEMSIO file,datatype,model is ',gfile%gtype, &
!        gfile%gdatatype,gfile%modelname,idate(1:7)
    if(present(version)) gfile%version=version
    if(present(dimx)) gfile%dimx=dimx
    if(present(dimy)) gfile%dimy=dimy
    if(present(dimz)) gfile%dimz=dimz
    if(present(nrec)) gfile%nrec=nrec
    if(present(nmeta)) gfile%nmeta=nmeta
    if(gfile%nmeta==nemsio_intfill) gfile%nmeta=12
    if(present(lmeta)) gfile%lmeta=lmeta
    if(gfile%lmeta==nemsio_intfill)   &
      gfile%lmeta=25*nemsio_intkind+4*nemsio_realkind+nemsio_logickind
    if(present(nsoil)) gfile%nsoil=nsoil
    if(gfile%nsoil.eq.nemsio_intfill) gfile%nsoil=4
    if(present(nframe)) gfile%nframe=nframe
    if(gfile%nframe.eq.nemsio_intfill) gfile%nframe=0
    if(trim(gfile%modelname)=='GFS')gfile%nframe=0
    if(present(idate)) gfile%idate=idate
    if ( gfile%idate(1) .lt. 50) then
        gfile%idate(1)=2000+gfile%idate(1)
    else if (gfile%idate(1) .lt. 100) then
        gfile%idate(1)=1999+gfile%idate(1)
    endif
    if ( gfile%idate(1).eq.nemsio_intfill) then
      print *,'idate=',gfile%idate,' WRONG: please provide idate(1:7)(yyyy/mm/dd/hh/min/secn/secd)!!!'
      call nemsio_stop()
    endif
!
    linit= gfile%dimx .eq. nemsio_intfill .or. gfile%dimy .eq. nemsio_intfill &
      .or. gfile%dimz .eq. nemsio_intfill .or. gfile%nrec .eq. nemsio_intfill &
      .or. gfile%nmeta .eq. 12 
!    
     
    if ( gfile%gtype(1:6).eq."NEMSIO" .and. linit ) then
      call nemsio_gfinit(gfile,ios,recname=recname,reclevtyp=reclevtyp,reclev=reclev)
      if (ios .ne.0 ) then
        iret=ios
        return
      endif
    endif
!
!------------------------------------------------------------
! set up basic gfile meta data variables from outsides to 
! define meta data array
!------------------------------------------------------------
    if(present(nfday)) gfile%nfday=nfday
    if(present(nfhour)) gfile%nfhour=nfhour
    if(present(nfminute)) gfile%nfminute=nfminute
    if(present(nfsecondn)) gfile%nfsecondn=nfsecondn
    if(present(nfsecondd)) gfile%nfsecondd=nfsecondd
    if(present(ntrac)) gfile%ntrac=ntrac
    if(gfile%ntrac.eq.nemsio_intfill) gfile%ntrac=0
    if(present(ncldt)) gfile%ncldt=ncldt
    if(present(jcap)) gfile%jcap=jcap
    if(present(idvc)) gfile%idvc=idvc
    if(present(idsl)) gfile%idsl=idsl
    if(present(idvm)) gfile%idvm=idvm
    if(present(idrt)) gfile%idrt=idrt
    if(present(rlon_min)) gfile%rlon_min=rlon_min
    if(present(rlon_max)) gfile%rlon_max=rlon_max
    if(present(rlat_min)) gfile%rlat_min=rlat_min
    if(present(rlat_max)) gfile%rlat_max=rlat_max
    if(present(extrameta)) gfile%extrameta=extrameta
    if(gfile%fieldsize.eq.nemsio_intfill) &
       gfile%fieldsize=(gfile%dimx+2*gfile%nframe)*(gfile%dimy+2*gfile%nframe)
    if(gfile%mype.eq.gfile%lead_task) then
     if(gfile%gdatatype.eq.'bin4') then
      call mpi_send(gfile%fieldsize*nemsio_realkind,1,MPI_integer,0,99,gfile%mpi_comm,ios)
      call mpi_recv(gfile%fieldsize_real4,1,MPI_real4,0,99,gfile%mpi_comm,status,ios)
     elseif(gfile%gdatatype.eq.'bin8') then
      call mpi_send(gfile%fieldsize*nemsio_dblekind,1,MPI_integer,0,99,gfile%mpi_comm,ios)
      call mpi_recv(gfile%fieldsize_real8,1,MPI_real8,0,99,gfile%mpi_comm,status,ios)
     endif
    endif
!
    if( gfile%extrameta )then
      if(present(nmetavari).and.nmetavari.gt.0.and.present(variname) &
         .and.size(variname).eq.nmetavari .and. &
         present(varival).and.size(varival).eq.nmetavari) then
           gfile%nmetavari=nmetavari
           if(allocated(gfile%variname)) deallocate(gfile%variname)
           if(allocated(gfile%varival)) deallocate(gfile%varival)
           allocate(gfile%variname(nmetavari),gfile%varival(nmetavari))
           gfile%variname=variname
           gfile%varival=varival
      endif
      if(present(nmetavarr).and.nmetavarr.gt.0.and.present(varrname) &
         .and.size(varrname).eq.nmetavarr .and. &
         present(varrval).and.size(varrval).eq.nmetavarr) then
           gfile%nmetavarr=nmetavarr
           if(allocated(gfile%varrname)) deallocate(gfile%varrname)
           if(allocated(gfile%varrval)) deallocate(gfile%varrval)
           allocate(gfile%varrname(nmetavarr),gfile%varrval(nmetavarr))
           gfile%varrname=varrname
           gfile%varrval=varrval
      endif
      if(present(nmetavarl).and.nmetavarl.gt.0.and.present(varlname) &
         .and.size(varlname).eq.nmetavarl .and. &
         present(varlval).and.size(varlval).eq.nmetavarl) then
           gfile%nmetavarl=nmetavarl
           if(allocated(gfile%varlname)) deallocate(gfile%varlname)
           if(allocated(gfile%varlval)) deallocate(gfile%varlval)
           allocate(gfile%varlname(nmetavarl),gfile%varlval(nmetavarl))
           gfile%varlname=varlname
           gfile%varlval=varlval
      endif
      if(present(nmetavarc).and.nmetavarc.gt.0.and.present(varcname) &
         .and.size(varcname).eq.nmetavarc .and. &
         present(varcval).and.size(varcval).eq.nmetavarc) then
           gfile%nmetavarc=nmetavarc
           if(allocated(gfile%varcname)) deallocate(gfile%varcname)
           if(allocated(gfile%varcval)) deallocate(gfile%varcval)
           allocate(gfile%varcname(nmetavarc),gfile%varcval(nmetavarc))
           gfile%varcname=varcname
           gfile%varcval=varcval
      endif
      if(present(nmetaaryi).and.nmetaaryi.gt.0.and.present(aryiname) &
         .and.size(aryiname).eq.nmetaaryi .and. &
         present(aryilen).and.size(aryilen).eq.nmetaaryi) then
           gfile%nmetaaryi=nmetaaryi
           if(allocated(gfile%aryiname)) deallocate(gfile%aryiname)
           if(allocated(gfile%aryilen)) deallocate(gfile%aryilen)
           allocate(gfile%aryiname(nmetaaryi),gfile%aryilen(nmetaaryi))
           gfile%aryiname=aryiname
           gfile%aryilen=aryilen
           if(present(aryival).and.size(aryival).eq.nmetaaryi*maxval(gfile%aryilen) ) then
             if(allocated(gfile%aryival)) deallocate(gfile%aryival)
             allocate(gfile%aryival(maxval(gfile%aryilen),nmetaaryi))
             gfile%aryival=aryival
           endif
      endif
      if(present(nmetaaryr).and.nmetaaryr.gt.0.and.present(aryrname) &
         .and.size(aryrname).eq.nmetaaryr .and. &
         present(aryrlen).and.size(aryrlen).eq.nmetaaryr) then
           gfile%nmetaaryr=nmetaaryr
           if(allocated(gfile%aryrname)) deallocate(gfile%aryrname)
           if(allocated(gfile%aryrlen)) deallocate(gfile%aryrlen)
           allocate(gfile%aryrname(nmetaaryr),gfile%aryrlen(nmetaaryr))
           gfile%aryrname=aryrname
           gfile%aryrlen=aryrlen
!           print *,'in wcreate,gfile%aryrname=',gfile%aryrname
!           print *,'in wcreate,gfile%aryrlen=',gfile%aryrlen
           if(present(aryrval).and.size(aryrval).eq.nmetaaryr*maxval(gfile%aryrlen)) then
             if(allocated(gfile%aryrval)) deallocate(gfile%aryrval)
             allocate(gfile%aryrval(maxval(gfile%aryrlen),nmetaaryr))
             gfile%aryrval=aryrval
           endif
      endif
      if(present(nmetaaryl).and.nmetaaryl.gt.0.and.present(arylname) &
          .and.size(arylname).eq.nmetaaryl .and. &
           present(aryllen).and.size(aryllen).eq.nmetaaryl) then
           gfile%nmetaaryl=nmetaaryl
           if(allocated(gfile%arylname)) deallocate(gfile%arylname)
           if(allocated(gfile%aryllen)) deallocate(gfile%aryllen)
           allocate(gfile%arylname(nmetaaryl),gfile%aryllen(nmetaaryl))
           gfile%arylname=arylname
           gfile%aryllen=aryllen
           if(present(arylval).and.size(arylval).eq.nmetaaryl*maxval(gfile%aryllen)) then
             if(allocated(gfile%arylval)) deallocate(gfile%arylval)
             allocate(gfile%arylval(maxval(gfile%aryllen),nmetaaryl))
             gfile%arylval=arylval
           endif
      endif
      if(present(nmetaaryc).and.nmetaaryc.gt.0.and.present(arycname) &
          .and.size(arycname).eq.nmetaaryc .and. &
          present(aryclen).and.size(aryclen).eq.nmetaaryc) then
           gfile%nmetaaryc=nmetaaryc
           if(allocated(gfile%arycname)) deallocate(gfile%arycname)
           if(allocated(gfile%aryclen)) deallocate(gfile%aryclen)
           allocate(gfile%arycname(nmetaaryc),gfile%aryclen(nmetaaryc))
           gfile%arycname=arycname
           gfile%aryclen=aryclen
           if(present(arycval).and.size(arycval).eq.nmetaaryc*maxval(gfile%aryclen)) then
             if(allocated(gfile%arycval)) deallocate(gfile%arycval)
             allocate(gfile%arycval(maxval(gfile%aryclen),nmetaaryc))
             gfile%arycval=arycval
           endif
      endif
      if (gfile%nmetavari+gfile%nmetavarr+gfile%nmetavarl+gfile%nmetavarc+ &
           gfile%nmetaaryi+gfile%nmetaaryr+gfile%nmetaaryl+gfile%nmetaaryc.lt.0 )then
           print *,'WRONG: gfile%extrameta is not compatiable with input extra meta!'
           return
      endif
    endif 
!------------------------------------------------------------
! check gfile meta data array size
!------------------------------------------------------------
!    print *,'before chkgfary'
    call nemsio_chkgfary(gfile,ios)
    if (ios.ne. 0) then
      iret=ios
      return
    endif
!------------------------------------------------------------
! continue to set gfile meta data variables tnd arrays
!------------------------------------------------------------
!set gfile data type to bin/grb, default set to grb
!recname
    if(present(recname) ) then
       if (gfile%nrec.ne.size(recname)) return
       gfile%recname=recname
    endif
!reclevtyp
    if(present(reclevtyp)) then
       if (gfile%nrec.ne.size(reclevtyp)) return
       gfile%reclevtyp=reclevtyp
    endif
!reclev
    if(present(reclev) ) then
       if (gfile%nrec.ne.size(reclev)) return
       gfile%reclev=reclev
    endif
!vcoord vcoord(levs+1
    if(present(vcoord) ) then
       if ((gfile%dimz+1)*3*2.ne.size(vcoord)) return
       gfile%vcoord=vcoord
    endif
!lat
    if(present(lat) ) then
       if (gfile%fieldsize.ne.size(lat)) return
       gfile%lat=lat
    endif
!lon
    if(present(lon) ) then
       if (gfile%fieldsize.ne.size(lon)) return
       gfile%lon=lon
    endif
!dx
    if(present(dx) ) then
       if (gfile%fieldsize.ne.size(dx)) return
       gfile%dx=dx
    endif
!dy
    if(present(dy) ) then
       if (gfile%fieldsize.ne.size(dy)) return
       gfile%dy=dy
    endif
!Cpi
    if( present(Cpi) ) then
       if (gfile%ntrac+1.ne.size(gfile%Cpi)) return
       gfile%Cpi = Cpi
    endif
!Ri
    if( present(Ri) ) then
       if (gfile%ntrac+1.ne.size(gfile%Ri)) return
       gfile%Ri = Ri
    endif
!------------------------------------------------------------
! write out the header by lead_task
!------------------------------------------------------------
    if(gfile%mype.eq.gfile%lead_task) then
      call baopenwt(gfile%flunit,gfile%gfname,ios)
      if ( ios.ne.0) return
!------------------------------------------------------------
! write out first meta data record
!------------------------------------------------------------
      meta1%gtype=gfile%gtype
      meta1%gdatatype=gfile%gdatatype
      meta1%modelname=gfile%modelname
      meta1%version=gfile%version
      meta1%nmeta=gfile%nmeta
      meta1%lmeta=gfile%lmeta
      meta1%reserve=0
      iskip=0
      iwrite=nemsio_lmeta1
      call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,meta1)
      if(nwrite.lt.iwrite) return
      gfile%tlmeta=nwrite
!      print *,'tlmet1 =',gfile%tlmeta,'nwrite=',nwrite,meta1%gdatatype
!------------------------------------------------------------
! write out second meta data record
!------------------------------------------------------------
      meta2%nrec=gfile%nrec
      meta2%idate(1:7)=gfile%idate(1:7)
      meta2%nfday=gfile%nfday
      meta2%nfhour=gfile%nfhour
      meta2%nfminute=gfile%nfminute
      meta2%nfsecondn=gfile%nfsecondn
      meta2%nfsecondd=gfile%nfsecondd
      meta2%dimx=gfile%dimx
      meta2%dimy=gfile%dimy
      meta2%dimz=gfile%dimz
      meta2%nframe=gfile%nframe
      meta2%nsoil=gfile%nsoil
      meta2%ntrac=gfile%ntrac
      meta2%jcap=gfile%jcap
      meta2%ncldt=gfile%ncldt
      meta2%idvc=gfile%idvc
      meta2%idsl=gfile%idsl
      meta2%idvm=gfile%idvm
      meta2%idrt=gfile%idrt
      meta2%rlon_min=gfile%rlon_min
      meta2%rlon_max=gfile%rlon_max
      meta2%rlat_min=gfile%rlat_min
      meta2%rlat_max=gfile%rlat_max
      meta2%extrameta=gfile%extrameta
     iskip=iskip+nwrite
     iwrite=gfile%lmeta
     call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,meta2)
     if(nwrite.lt.iwrite) return
     gfile%tlmeta=gfile%tlmeta+nwrite
!     print *,'tlmet2 =',gfile%tlmeta,'nwrite=',nwrite,'meta2=', &
!      meta2%dimx,meta2%dimy,meta2%dimz,meta2%nframe,meta2%nsoil, &
!      meta2%ntrac,meta2%jcap,meta2%ncldt,meta2%idvc,meta2%idsl,  &
!      meta2%idvm,meta2%idrt,meta2%rlon_min,meta2%rlon_max,       &
!      meta2%rlat_min,meta2%rlat_max,meta2%extrameta
!------------------------------------------------------------
! write out 3rd-13th meta data record (arrays)
!------------------------------------------------------------
!recname
     iskip=iskip+nwrite
     iwrite=nemsio_charkind*size(gfile%recname)
     call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%recname)
     if(nwrite.lt.iwrite) return
     gfile%tlmeta=gfile%tlmeta+nwrite
!     print *,'tlmetrecname =',gfile%tlmeta,'nwrite=',nwrite
!reclevtyp
     iskip=iskip+nwrite
     iwrite=nemsio_charkind*size(gfile%reclevtyp)
     call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%reclevtyp)
     if(nwrite.lt.iwrite) return
     gfile%tlmeta=gfile%tlmeta+nwrite
!     print *,'tlmetreclevty=',gfile%tlmeta,'nwrite=',nwrite
!reclev
     iskip=iskip+nwrite
     iwrite=nemsio_intkind*size(gfile%reclev)
     call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%reclev)
     if(nwrite.lt.iwrite) return
     gfile%tlmeta=gfile%tlmeta+nwrite
!     print *,'tlmetreclev=',gfile%tlmeta,'nwrite=',nwrite
!vcoord
    nummeta=gfile%nmeta-5
    if ( nummeta.gt.0 ) then
      iskip=iskip+nwrite
      iwrite=nemsio_realkind*size(gfile%vcoord)
      call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%vcoord)
      if(nwrite.lt.iwrite) return
      gfile%tlmeta=gfile%tlmeta+nwrite
!      print *,'tlmetavcoord=',gfile%tlmeta,'nwrite=',nwrite
      nummeta=nummeta-1
    endif
!lat
    if ( nummeta.gt.0 ) then
      iskip=iskip+nwrite
      iwrite=nemsio_realkind*size(gfile%lat)
      call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%lat)
      if(nwrite.lt.iwrite) return
      gfile%tlmeta=gfile%tlmeta+nwrite
!      print *,'tlmetreclat=',gfile%tlmeta,'nwrite=',nwrite
      nummeta=nummeta-1
    endif
!lon
    if ( nummeta.gt.0 ) then
      iskip=iskip+nwrite
      iwrite=nemsio_realkind*size(gfile%lon)
      call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%lon)
      if(nwrite.lt.iwrite) return
      gfile%tlmeta=gfile%tlmeta+nwrite
!      print *,'tlmetreclon=',gfile%tlmeta,'nwrite=',nwrite
      nummeta=nummeta-1
    endif
!dx
    if ( nummeta.gt.0 ) then
      iskip=iskip+nwrite
      iwrite=nemsio_realkind*size(gfile%dx)
      call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%dx)
      if(nwrite.lt.iwrite) return
      gfile%tlmeta=gfile%tlmeta+nwrite
!      print *,'tlmetrecdx=',gfile%tlmeta,'nwrite=',nwrite,  &
!        maxval(gfile%dx),minval(gfile%dx),maxval(gfile%dy),maxval(gfile%dy)
      nummeta=nummeta-1
    endif
!dy
    if ( nummeta.gt.0 ) then
      iskip=iskip+nwrite
      iwrite=nemsio_realkind*size(gfile%dy)
      call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%dy)
      if(nwrite.lt.iwrite) return
      gfile%tlmeta=gfile%tlmeta+nwrite
!      print *,'tlmetrecdy=',gfile%tlmeta,'nwrite=',nwrite
      nummeta=nummeta-1
    endif
!Cpi
    if ( nummeta.gt.0 ) then
      iskip=iskip+nwrite
      iwrite=nemsio_realkind*size(gfile%Cpi)
      call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%Cpi)
      if(nwrite.lt.iwrite) return
      gfile%tlmeta=gfile%tlmeta+nwrite
!      print *,'tlmetreccpi=',gfile%tlmeta,'nwrite=',nwrite
      nummeta=nummeta-1
    endif
!Ri
    if ( nummeta.gt.0 ) then
      iskip=iskip+nwrite
      iwrite=nemsio_realkind*size(gfile%Ri)
      call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%Ri)
      if(nwrite.lt.iwrite) return
      gfile%tlmeta=gfile%tlmeta+nwrite
!      print *,'tlmetrecri=',gfile%tlmeta,'nwrite=',nwrite
      nummeta=nummeta-1
    endif
!------------------------------------------------------------
! write out extra meta data record 
!------------------------------------------------------------
    if(gfile%extrameta) then
      meta3%nmetavari=gfile%nmetavari
      meta3%nmetavarr=gfile%nmetavarr
      meta3%nmetavarl=gfile%nmetavarl
      meta3%nmetavarc=gfile%nmetavarc
      meta3%nmetaaryi=gfile%nmetaaryi
      meta3%nmetaaryr=gfile%nmetaaryr
      meta3%nmetaaryl=gfile%nmetaaryl
      meta3%nmetaaryc=gfile%nmetaaryc
      iskip=iskip+nwrite
      iwrite=nemsio_lmeta3
      call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,meta3)
      if(nwrite.lt.iwrite) return
      gfile%tlmeta=gfile%tlmeta+nwrite
!      print *,'tlmetameta3=',gfile%tlmeta
!
!-- write meta var integer
      if (gfile%nmetavari.gt.0) then
        iskip=iskip+nwrite
        iwrite=nemsio_charkind*gfile%nmetavari
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%variname)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite

        iskip=iskip+nwrite
        iwrite=nemsio_intkind*gfile%nmetavari
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%varival)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
!        print *,'rlmetavari=',gfile%tlmeta
      endif
      if (gfile%nmetavarr.gt.0) then
        iskip=iskip+nwrite
        iwrite=nemsio_charkind*gfile%nmetavarr
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%varrname)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite

        iskip=iskip+nwrite
        iwrite=nemsio_realkind*gfile%nmetavarr
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%varrval)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
      endif
      if (gfile%nmetavarl.gt.0) then
        iskip=iskip+nwrite
        iwrite=nemsio_charkind*gfile%nmetavarl
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%varlname)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
!        print *,'tlmetavarl =',gfile%tlmeta,'nwrite=',nwrite

        iskip=iskip+nwrite
        iwrite=nemsio_logickind*gfile%nmetavarl
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%varlval)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
!        print *,'tlmetavarlval =',gfile%tlmeta,'nwrite=',nwrite
      endif
      if (gfile%nmetavarc.gt.0) then
        iskip=iskip+nwrite
        iwrite=nemsio_charkind*gfile%nmetavarc
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%varcname)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
!        print *,'tlmetaaryinam =',gfile%tlmeta,'write=',nwrite

        iskip=iskip+nwrite
        iwrite=nemsio_charkind*gfile%nmetavarc
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%varcval)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
!        print *,'tlmetaaryilen =',gfile%tlmeta,'nwrite=',nwrite
      endif
!meta arr integer
      if (gfile%nmetaaryi.gt.0) then
        iskip=iskip+nwrite
        iwrite=nemsio_charkind*gfile%nmetaaryi
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%aryiname)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite

        iskip=iskip+nwrite
        iwrite=nemsio_intkind*gfile%nmetaaryi
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%aryilen)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
        do i=1,gfile%nmetaaryi
          iskip=iskip+nwrite
          iwrite=nemsio_intkind*gfile%aryilen(i)
          call bafrwritel(gfile%flunit,iskip,iwrite,nwrite, &
                         gfile%aryival(1:gfile%aryilen(i),i))
          if(nwrite.lt.iwrite) return
          gfile%tlmeta=gfile%tlmeta+nwrite
!          print *,'tlmetaryint=',i,gfile%tlmeta,'nwrite=',nwrite
        enddo
!        print *,'after tlmetaryi ',gfile%nmetaaryr,gfile%nmetaaryl,gfile%nmetaaryc
      endif
!meta arr real
      if (gfile%nmetaaryr.gt.0) then
        iskip=iskip+nwrite
        iwrite=nemsio_charkind*gfile%nmetaaryr
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%aryrname)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
!!          print *,'before tlmetaryr 1'

        iskip=iskip+nwrite
        iwrite=nemsio_intkind*gfile%nmetaaryr
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%aryrlen)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
        do i=1,gfile%nmetaaryr
          iskip=iskip+nwrite
          iwrite=nemsio_intkind*gfile%aryrlen(i)
          call bafrwritel(gfile%flunit,iskip,iwrite,nwrite, &
                         gfile%aryrval(1:gfile%aryrlen(i),i))
          if(nwrite.lt.iwrite) return
          gfile%tlmeta=gfile%tlmeta+nwrite
!          print *,'tlmetaryreal=',i,gfile%tlmeta,'nwrite=',nwrite
        enddo
      endif
!meta arr logical
      if (gfile%nmetaaryl.gt.0) then
        iskip=iskip+nwrite
        iwrite=nemsio_charkind*gfile%nmetaaryl
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%arylname)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite

        iskip=iskip+nwrite
        iwrite=nemsio_intkind*gfile%nmetaaryl
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%aryllen)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
        do i=1,gfile%nmetaaryl
          iskip=iskip+nwrite
          iwrite=nemsio_logickind*gfile%aryllen(i)
          call bafrwritel(gfile%flunit,iskip,iwrite,nwrite, &
                         gfile%arylval(1:gfile%aryllen(i),i))
          if(nwrite.lt.iwrite) return
          gfile%tlmeta=gfile%tlmeta+nwrite
!          print *,'tlmetarylogic=',i,gfile%tlmeta,'nwrite=',nwrite
        enddo
      endif
!meta arr char
      if (gfile%nmetaaryc.gt.0) then
        iskip=iskip+nwrite
        iwrite=nemsio_charkind*gfile%nmetaaryc
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%arycname)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite

        iskip=iskip+nwrite
        iwrite=nemsio_intkind*gfile%nmetaaryc
        call bafrwritel(gfile%flunit,iskip,iwrite,nwrite,gfile%aryclen)
        if(nwrite.lt.iwrite) return
        gfile%tlmeta=gfile%tlmeta+nwrite
        do i=1,gfile%nmetaaryc
          iskip=iskip+nwrite
          iwrite=nemsio_charkind*gfile%aryclen(i)
          call bafrwritel(gfile%flunit,iskip,iwrite,nwrite, &
                         gfile%arycval(1:gfile%aryclen(i),i))
          if(nwrite.lt.iwrite) return
          gfile%tlmeta=gfile%tlmeta+nwrite
!          print *,'tlmetarycogic=',i,gfile%tlmeta,'nwrite=',nwrite
        enddo
      endif
    endif     !end of gfile%extrameta
   endif      !end of lead_task
!mpi
    call MPI_Barrier(gfile%mpi_comm, ios)
    call mpi_bcast(gfile%tlmeta,1,MPI_INTEGER,gfile%lead_task,gfile%mpi_comm,ios)
!    print *,'after mpi_bcasttlmeta,',gfile%tlmeta, 'end of wcreate'
!
    iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine nemsio_wcreate
!------------------------------------------------------------------------------
  subroutine nemsio_getfilehead(gfile,iret,gtype,gdatatype,gfname,gaction, &
      modelname,version,nmeta,lmeta,nrec,idate,nfday,nfhour,nfminute, &
      nfsecondn,nfsecondd,dimx,dimy,dimz,nframe,nsoil,ntrac,ncldt,jcap,&
      idvc,idsl,idvm,idrt, rlon_min,rlon_max,rlat_min,rlat_max,extrameta, &
      nmetavari,nmetavarr,nmetavarl,nmetavarc,nmetaaryi,nmetaaryr, &
      nmetaaryl,nmetaaryc,    &
      recname,reclevtyp,reclev,vcoord,lon,lat,dx,dy,cpi,ri, &
      variname,varival,varrname,varrval,varlname,varlval,varcname,varcval, &
      aryiname,aryilen,aryival,aryrname,aryrlen,aryrval,    &
      arylname,aryllen,arylval,arycname,aryclen,arycval    )

!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: get nemsio meta data information from outside
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                :: gfile
    integer(nemsio_intkind),optional,intent(out) :: iret
    character*(*),optional,intent(out)           :: gtype,gdatatype,gfname, &
                                                    gaction,modelname
    integer(nemsio_intkind),optional,intent(out) :: version,nmeta,lmeta
    integer(nemsio_realkind),optional,intent(out):: nrec,idate(7),nfday,nfhour, &
                                                    nfminute,nfsecondn,nfsecondd
    integer(nemsio_realkind),optional,intent(out):: dimx,dimy,dimz,nframe, &
                                                    nsoil,ntrac
    integer(nemsio_realkind),optional,intent(out):: ncldt,jcap,idvc,idsl,idvm,idrt
    real(nemsio_realkind),optional,intent(out)   :: rlon_min,rlon_max,rlat_min, &
                                                    rlat_max
    logical(nemsio_logickind),optional,intent(out):: extrameta
    integer(nemsio_realkind),optional,intent(out):: nmetavari,nmetavarr, &
                                                    nmetavarl,nmetavarc,nmetaaryi, &
                                                    nmetaaryr,nmetaaryl,nmetaaryc
    character(*),optional,intent(out)           :: recname(:)
    character(*),optional,intent(out)           :: reclevtyp(:)
    integer(nemsio_intkind),optional,intent(out) :: reclev(:)
    real(nemsio_realkind),optional,intent(out)   :: vcoord(:,:,:)
    real(nemsio_realkind),optional,intent(out)   :: lat(:),lon(:)
    real(nemsio_realkind),optional,intent(out)   :: dx(:),dy(:)
    real(nemsio_realkind),optional,intent(out)   :: Cpi(:),Ri(:)
    character(*),optional,intent(out)            :: variname(:),varrname(:)
    character(*),optional,intent(out)            :: varlname(:),varcname(:)
    character(*),optional,intent(out)            :: aryiname(:),aryrname(:)
    character(*),optional,intent(out)            :: arylname(:),arycname(:)
    integer(nemsio_intkind),optional,intent(out) :: aryilen(:),aryrlen(:)
    integer(nemsio_intkind),optional,intent(out) :: aryllen(:),aryclen(:)
    integer(nemsio_intkind),optional,intent(out) :: varival(:),aryival(:,:)
    real(nemsio_realkind),optional,intent(out)   :: varrval(:),aryrval(:,:)
    logical(nemsio_logickind),optional,intent(out):: varlval(:),arylval(:,:)
    character(*),optional,intent(out):: varcval(:)
    character(*),optional,intent(out):: arycval(:)
    integer ierr
!------------------------------------------------------------
    if (present(iret)) iret=-3
    if(present(gtype)) gtype=gfile%gtype
    if(present(gdatatype)) gdatatype=gfile%gdatatype
    if(present(gfname)) gfname=trim(gfile%gfname)
    if(present(gaction)) gaction=gfile%gaction
    if(present(modelname)) modelname=gfile%modelname
    if(present(version)) version=gfile%version
    if(present(nmeta)) nmeta=gfile%nmeta
    if(present(lmeta)) lmeta=gfile%lmeta
    if(present(nrec)) nrec=gfile%nrec
    if(present(nfday)) nfday=gfile%nfday
    if(present(nfhour)) nfhour=gfile%nfhour
    if(present(nfminute)) nfminute=gfile%nfminute
    if(present(nfsecondn)) nfsecondn=gfile%nfsecondn
    if(present(nfsecondd)) nfsecondd=gfile%nfsecondd
    if(present(idate)) idate=gfile%idate
    if(present(dimx)) dimx=gfile%dimx
    if(present(dimy)) dimy=gfile%dimy
    if(present(dimz)) dimz=gfile%dimz
    if(present(nframe)) nframe=gfile%nframe
    if(present(nsoil)) nsoil=gfile%nsoil
    if(present(ntrac)) ntrac=gfile%ntrac
    if(present(jcap)) jcap=gfile%jcap
    if(present(ncldt)) ncldt=gfile%ncldt
    if(present(idvc)) idvc=gfile%idvc
    if(present(idsl)) idsl=gfile%idsl
    if(present(idvm)) idvm=gfile%idvm
    if(present(idrt)) idrt=gfile%idrt
    if(present(rlon_min)) rlon_min=gfile%rlon_min
    if(present(rlon_max)) rlon_max=gfile%rlon_max
    if(present(rlat_min)) rlat_min=gfile%rlat_min
    if(present(rlat_max)) rlat_max=gfile%rlat_max
    if(present(rlat_max)) rlat_max=gfile%rlat_max
    if(present(extrameta)) extrameta=gfile%extrameta
!
!    print *,'in getfilehead, 1extrameta=',gfile%extrameta,        &
!     'nrec=',gfile%nrec,'size(recname)=',size(recname),           &
!     size(reclevtyp),size(reclev)
!--- rec
    if(present(recname) ) then
       if (gfile%nrec.ne.size(recname)) then
         if ( present(iret)) return
         call nemsio_stop
       else
         recname=gfile%recname
       endif
    endif
    if(present(reclevtyp)) then
       if (gfile%nrec.ne.size(reclevtyp)) then
         if ( present(iret)) return
         call nemsio_stop
       else
         reclevtyp=gfile%reclevtyp
       endif
    endif
    if(present(reclev) ) then
       if (gfile%nrec.ne.size(reclev)) then
         if ( present(iret)) return
         call nemsio_stop
       else
         reclev=gfile%reclev
       endif
    endif
!--- vcoord
    if(present(vcoord)) then
       if (size(vcoord) .ne. (gfile%dimz+1)*2*3 ) then
         if ( present(iret))  return
         call nemsio_stop
       else
         vcoord=gfile%vcoord
       endif
    endif
!--- lat
    if(present(lat) ) then
       if (size(lat).ne.gfile%fieldsize) then
         if ( present(iret))  return
         call nemsio_stop
       else
         lat=gfile%lat
       endif
    endif
    if(present(lon) ) then
       if (size(lon).ne.gfile%fieldsize) then
         if ( present(iret)) return
         call nemsio_stop
       else
         lon=gfile%lon
       endif
    endif
!--- dx
    if(present(dx) ) then
!       print *,'getfilehead, size(dx)=',size(dx),gfile%fieldsize,  &
!          maxval(gfile%dx),minval(gfile%dx)
       if (size(dx).ne.gfile%fieldsize) then
         if ( present(iret))  return
         call nemsio_stop
       else
         dx=gfile%dx
       endif
    endif
    if(present(dy) ) then
       if (size(dy).ne.gfile%fieldsize) then
         if ( present(iret)) return
         call nemsio_stop
       else
         dy=gfile%dy
       endif
    endif
!--- Cpi
    if(present(Cpi) ) then
       if (gfile%ntrac+1.ne.size(Cpi)) then
         if ( present(iret)) return
         call nemsio_stop
       else
         Cpi=gfile%Cpi
       endif
    endif
!Ri
    if(present(Ri) ) then 
       if (gfile%ntrac+1.ne.size(Ri)) then
         if ( present(iret)) return
         call nemsio_stop
       else
         Ri=gfile%Ri
       endif
    endif
!------------------------------------------------------------------------------
!*** for extra meta field
!------------------------------------------------------------------------------
!extrameta
    if(present(extrameta) ) extrameta=gfile%extrameta
    if(gfile%extrameta) then
      if (present(nmetavari) ) nmetavari=gfile%nmetavari
      if (present(nmetavarr) ) nmetavarr=gfile%nmetavarr
      if (present(nmetavarl) ) nmetavarl=gfile%nmetavarl
      if (present(nmetaaryi) ) nmetaaryi=gfile%nmetaaryi
      if (present(nmetaaryr) ) nmetaaryr=gfile%nmetaaryr
      if (present(nmetaaryl) ) nmetaaryl=gfile%nmetaaryl
      if ( gfile%nmetavari.gt.0 ) then
         if (present(variname).and.size(variname).eq.nmetavari) &
             variname=gfile%variname
         if (present(varival).and.size(varival).eq.nmetavari) &
             varival=gfile%varival
      endif
      if ( gfile%nmetavarr.gt.0 ) then
         if (present(varrname).and.size(varrname).eq.nmetavarr) &
             varrname=gfile%varrname
         if (present(varrval).and.size(varrval).eq.nmetavarr) &
             varrval=gfile%varrval
      endif
      if ( gfile%nmetavarl.gt.0 ) then
         if (present(varlname).and.size(varlname).eq.nmetavarl) &
             varlname=gfile%varlname
         if (present(varlval).and.size(varlval).eq.nmetavarl) &
             varlval=gfile%varlval
      endif
      if ( gfile%nmetaaryi.gt.0 ) then
         if (present(aryiname).and.size(aryiname).eq.nmetaaryi) &
             aryiname=gfile%aryiname
         if (present(aryilen).and.size(aryilen).eq.nmetaaryi) &
             aryilen=gfile%aryilen
         if (present(aryival).and.size(aryival).eq.nmetaaryi*maxval(gfile%aryilen) ) &
             aryival=gfile%aryival
      endif
      if ( gfile%nmetaaryr.gt.0 ) then
         if (present(aryrname).and.size(aryrname).eq.nmetaaryr) &
             aryiname=gfile%aryiname
         if (present(aryrlen).and.size(aryrlen).eq.nmetaaryr) &
             aryrlen=gfile%aryrlen
         if (present(aryrval).and.size(aryrval).eq.nmetaaryr*maxval(gfile%aryrlen) ) &
             aryrval=gfile%aryrval
      endif
      if ( gfile%nmetaaryl.gt.0 ) then
         if (present(arylname).and.size(arylname).eq.nmetaaryl) &
             arylname=gfile%arylname
         if (present(aryllen).and.size(aryllen).eq.nmetaaryl) &
             aryllen=gfile%aryllen
         if (present(arylval).and.size(arylval).eq.nmetaaryl*maxval(gfile%aryllen) ) &
             arylval=gfile%arylval
      endif
    endif

    call mpi_barrier(gfile%mpi_comm,ierr)
    if ( ierr.ne.0) return
    if ( present(iret)) iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine nemsio_getfilehead
!------------------------------------------------------------------------------
   subroutine nemsio_getfheadvari(gfile,varname,varval,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: get meta data var value from file header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                 :: gfile
    character(*),  intent(in)                     :: varname
    integer(nemsio_intkind),intent(out)           :: varval
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer i,j
!---
    if(present(iret) ) iret=-17
    do i=1,gfile%headvarinum
      if(equal_str_nocase(trim(varname),trim(gfile%headvariname(i))) ) then
           varval=gfile%headvarival(i)
           if(present(iret) ) iret=0
           return
      endif
    enddo
!---
    if(gfile%nmetavari.gt.0) then
      do i=1,gfile%nmetavari
        if(equal_str_nocase(trim(varname),trim(gfile%variname(i))) ) then
           varval=gfile%varival(i)
           if(present(iret) ) iret=0
           return
        endif
      enddo
    endif
!---    
    if(.not.present(iret) ) call nemsio_stop
    return
  end subroutine nemsio_getfheadvari
!------------------------------------------------------------------------------
   subroutine nemsio_getfheadvarr(gfile,varname,varval,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: get meta data var value from file header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                 :: gfile
    character(*),  intent(in)                     :: varname
    real(nemsio_realkind),intent(out)             :: varval
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer i,j
!---
    if(present(iret) ) iret=-17
    do i=1,gfile%headvarrnum
      if(equal_str_nocase(trim(varname),trim(gfile%headvarrname(i))) ) then
           varval=gfile%headvarrval(i)
           if(present(iret) ) iret=0
           return
      endif
    enddo
!---
    if(gfile%nmetavarr.gt.0) then
      do i=1,gfile%nmetavarr
        if(equal_str_nocase(trim(varname),trim(gfile%varrname(i))) ) then
           varval=gfile%varrval(i)
           if(present(iret) ) iret=0
           return
        endif
      enddo
    endif
!---
    if(.not.present(iret) ) call nemsio_stop
    return
  end subroutine nemsio_getfheadvarr
!------------------------------------------------------------------------------
   subroutine nemsio_getfheadvarl(gfile,varname,varval,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: get meta data var value from file header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                 :: gfile
    character(*),  intent(in)                     :: varname
    logical(nemsio_logickind),intent(out)         :: varval
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer i,j
!---
    if(present(iret) ) iret=-17
    if(gfile%nmetavarl.gt.0) then
      do i=1,gfile%nmetavarl
        if(equal_str_nocase(trim(varname),trim(gfile%varlname(i))) ) then
           varval=gfile%varlval(i)
           if(present(iret) ) iret=0
           return
        endif
      enddo
    endif
!---
    if(.not.present(iret) ) call nemsio_stop
    return
  end subroutine nemsio_getfheadvarl
!------------------------------------------------------------------------------
   subroutine nemsio_getfheadvarc(gfile,varname,varval,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: get meta data var value from file header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                 :: gfile
    character(*),  intent(in)                     :: varname
    character(*),intent(out)                      :: varval
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer i,j
!---
    if(present(iret) ) iret=-17
    do i=1,gfile%headvarcnum
      if(equal_str_nocase(trim(varname),trim(gfile%headvarcname(i))) ) then
           varval=gfile%headvarcval(i)
           if(present(iret) ) iret=0
           return
      endif
    enddo
!---
    if(gfile%nmetavarc.gt.0) then
      do i=1,gfile%nmetavarc
        if(equal_str_nocase(trim(varname),trim(gfile%varcname(i))) ) then
           varval=gfile%varcval(i)
           if(present(iret) ) iret=0
           return
        endif
      enddo
    endif
!---
    if(.not.present(iret) ) call nemsio_stop
    return
  end subroutine nemsio_getfheadvarc
!------------------------------------------------------------------------------
  subroutine nemsio_getfheadaryi(gfile,varname,varval,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: get meta data var value from file header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                 :: gfile
    character(*),  intent(in)                     :: varname
    integer(nemsio_intkind),intent(out)           :: varval(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer i,j,ierr
!---
    if(present(iret) ) iret=-17
    do i=1,gfile%headaryinum
      if(equal_str_nocase(trim(varname),trim(gfile%headaryiname(i))) ) then
           varval(:)=gfile%headaryival(1:gfile%aryilen(i),i)
           if(present(iret) ) iret=0
           return
      endif
    enddo
!---
    if(gfile%nmetaaryi.gt.0) then
      do i=1,gfile%nmetaaryi
        if(equal_str_nocase(trim(varname),trim(gfile%aryiname(i))) ) then
           varval(:)=gfile%aryival(1:gfile%aryilen(i),i)
           if(present(iret) ) iret=0
           ierr=0
           return
        endif
      enddo
    endif
!---    
    if(.not.present(iret) ) call nemsio_stop
    return
  end subroutine nemsio_getfheadaryi
!------------------------------------------------------------------------------
   subroutine nemsio_getfheadaryr(gfile,varname,varval,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: get meta data var value from file header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                 :: gfile
    character(*),  intent(in)                     :: varname
    real(nemsio_realkind),intent(out)             :: varval(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer i,j,ierr
!---
    if(present(iret) ) iret=-17
    if(gfile%headaryrnum>0) then
     do i=1,gfile%headaryrnum
      if(equal_str_nocase(trim(varname),trim(gfile%headaryrname(i))) ) then
           varval(:)=gfile%headaryrval(1:gfile%aryrlen(i),i)
           if(present(iret) ) iret=0
           return
      endif
     enddo
    endif
!---
    if(gfile%nmetaaryr.gt.0) then
      do i=1,gfile%nmetaaryr
        if(equal_str_nocase(trim(varname),trim(gfile%aryrname(i)))) then
           varval(:)=gfile%aryrval(1:gfile%aryrlen(i),i)
           if(present(iret) ) iret=0
           ierr=0
           return
        endif
      enddo
    endif
!---
    if(.not.present(iret) ) call nemsio_stop
    return
  end subroutine nemsio_getfheadaryr
!------------------------------------------------------------------------------
   subroutine nemsio_getfheadaryl(gfile,varname,varval,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: get meta data var value from file header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                 :: gfile
    character(*),  intent(in)                     :: varname
    logical(nemsio_logickind),intent(out)         :: varval(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer i,j,ierr
!---
    if(present(iret) ) iret=-17
    if(gfile%nmetaaryl.gt.0) then
      do i=1,gfile%nmetaaryl
        if(equal_str_nocase(trim(varname),trim(gfile%arylname(i)))) then
           varval(:)=gfile%arylval(1:gfile%aryllen(i),i)
           if(present(iret) ) iret=0
           ierr=0
           return
        endif
      enddo
    endif
!---
    if(.not.present(iret) ) call nemsio_stop
    return
  end subroutine nemsio_getfheadaryl
!------------------------------------------------------------------------------
   subroutine nemsio_getfheadaryc(gfile,varname,varval,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: get meta data var value from file header
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                 :: gfile
    character(*),  intent(in)                     :: varname
    character(*),intent(out)                      :: varval(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer i,j,ierr
!---
    if(present(iret) ) iret=-17
    if(gfile%nmetaaryc.gt.0) then
      do i=1,gfile%nmetaaryc
       if(equal_str_nocase(trim(varname),trim(gfile%headarycname(i))) ) then
           varval(:)=gfile%headarycval(1:gfile%aryclen(i),i)
           if(present(iret) ) iret=0
           return
       endif
      enddo
    endif
!---
    if(gfile%nmetaaryc.gt.0) then
      do i=1,gfile%nmetaaryc
        if(equal_str_nocase(trim(varname),trim(gfile%arycname(i)))) then
           varval(:)=gfile%arycval(1:gfile%aryclen(i),i)
           if(present(iret) ) iret=0
           ierr=0
           return
        endif
      enddo
    endif
!---
    if(.not.present(iret) ) call nemsio_stop
    return
  end subroutine nemsio_getfheadaryc
!------------------------------------------------------------------------------

!*****************   read bin data set :  ********************************
!
!------------------------------------------------------------------------------
  subroutine nemsio_readrecbin4(gfile,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_realkind),intent(inout)           :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!
    real(nemsio_dblekind),allocatable   :: data8(:)
!--
    if(trim(gfile%gdatatype)=='bin4') then
      call nemsio_readrecbind4(gfile,ista,iend,jsta,jend,jrec,data,iret)
    elseif(trim(gfile%gdatatype)=='bin8') then
      allocate(data8(size(data)))
      call nemsio_readrecbind8(gfile,ista,iend,jsta,jend,jrec,data8,iret)
      data=data8
      deallocate(data8)
    endif
!
    return
  end subroutine nemsio_readrecbin4
!------------------------------------------------------------------------------
  subroutine nemsio_readrecbin8(gfile,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_dblekind),intent(inout)           :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!
    real(nemsio_realkind),allocatable :: data4(:)
!--
    if(trim(gfile%gdatatype)=='bin4') then
      allocate(data4(size(data)))
      call nemsio_readrecbind4(gfile,ista,iend,jsta,jend,jrec,data4,iret)
      data=data4
      deallocate(data4)
    elseif(trim(gfile%gdatatype)=='bin8') then
      call nemsio_readrecbind8(gfile,ista,iend,jsta,jend,jrec,data,iret)
    endif
!
    return
  end subroutine nemsio_readrecbin8
!------------------------------------------------------------------------------
  subroutine nemsio_readrecvbin4(gfile,name,levtyp,lev,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    character(*),intent(in)                       :: name, levtyp
    integer(nemsio_intkind),intent(in)            :: lev
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_realkind),intent(inout)           :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!
    real(nemsio_dblekind),allocatable   :: data8(:)
!--
    if(trim(gfile%gdatatype)=='bin4') then
      call nemsio_readrecvbind4(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
    elseif(trim(gfile%gdatatype)=='bin8') then
      allocate(data8(size(data)))
      call nemsio_readrecvbind8(gfile,name,levtyp,lev,ista,iend,jsta,jend,data8,iret)
      data=data8
      deallocate(data8)
    endif
!
    return
  end subroutine nemsio_readrecvbin4
!------------------------------------------------------------------------------
  subroutine nemsio_readrecvbin8(gfile,name,levtyp,lev,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    character(*),intent(in)                       :: name, levtyp
    integer(nemsio_intkind),intent(in)            :: lev
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_dblekind),intent(inout)           :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!
    real(nemsio_realkind),allocatable :: data4(:)
!--
    if(gfile%gdatatype=='bin4') then
      allocate(data4(size(data)))
      call nemsio_readrecvbind4(gfile,name,levtyp,lev,ista,iend,jsta,jend,data4,iret)
      data=data4
      deallocate(data4)
    elseif(gfile%gdatatype=='bin8') then
      call nemsio_readrecvbind8(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
    endif
!
    return
  end subroutine nemsio_readrecvbin8

!------------------------------------------------------------------------------
  subroutine nemsio_readrecbind4(gfile,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_realkind),intent(inout)           :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!local vars
    integer,allocatable  :: fieldmap(:)
    integer fieldmapsize
    integer ios,i,j
    integer(8) idispstt
!
    iret=-11
!
!--- set file map
     fieldmapsize=(iend-ista+1)*(jend-jsta+1)
     allocate(fieldmap(fieldmapsize) )
     call set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,iret,1)
     if(iret.ne.0) return
     idispstt=int(jrec-1,8)*int(gfile%fieldsize*4+8,8)
!---
     call readmpi4(gfile,fieldmapsize,fieldmap,data,iret,idispstt)
     if(iret.ne.0) return
     deallocate(fieldmap)
!---
     iret=0
!
    return
  end subroutine nemsio_readrecbind4
!------------------------------------------------------------------------------
  subroutine nemsio_readrecvbind4(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    character(*),intent(in)                      :: name, levtyp
    integer(nemsio_intkind),intent(in)            :: lev
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    real(nemsio_realkind),intent(out)             :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!
    integer,allocatable  :: fieldmap(:)
    integer :: fieldmapsize
    integer :: jrec, ierr
    integer(8) idispstt
!
    iret=-12
    call nemsio_searchrecv(gfile,jrec,name,levtyp,lev,ierr)
    if ( ierr .ne. 0) return
!
!--- set file map
    fieldmapsize=(iend-ista+1)*(jend-jsta+1)
    allocate(fieldmap(fieldmapsize) )
    call set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,ierr,1)
    if(ierr.ne.0) return
!---
    idispstt=int(jrec-1,8)*int(gfile%fieldsize*4+8,8)
    call readmpi4(gfile,fieldmapsize,fieldmap,data,ierr,idispstt)
    if(ierr.ne.0) return
    deallocate(fieldmap)
!
    iret=0
!
    return
  end subroutine nemsio_readrecvbind4
!------------------------------------------------------------------------------
  subroutine nemsio_readrecbind8(gfile,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_dblekind),intent(out)             :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!local vars
    integer,allocatable  :: fieldmap(:)
    integer fieldmapsize
    integer ierr,i,j
    integer(8) idispstt
!
    iret=-11
!
!--- set file map
     fieldmapsize=(iend-ista+1)*(jend-jsta+1)
     allocate(fieldmap(fieldmapsize) )
     call set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,ierr,1)
     if(ierr.ne.0) return
!---
     idispstt=int(jrec-1,8)*int(gfile%fieldsize*8+8,8)
     call readmpi8(gfile,fieldmapsize,fieldmap,data,ierr)
     if(ierr.ne.0) return
     deallocate(fieldmap)
!---
     iret=0
!
    return
  end subroutine nemsio_readrecbind8
!------------------------------------------------------------------------------
  subroutine nemsio_readrecvbind8(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    character(*),intent(in)                      :: name, levtyp
    integer(nemsio_intkind),intent(in)            :: lev
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    real(nemsio_dblekind),intent(out)             :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!local vars
    integer,allocatable  :: fieldmap(:)
    integer fieldmapsize
    integer :: jrec, ierr
    integer(8) idispstt

    iret=-11
    call nemsio_searchrecv(gfile,jrec,name,levtyp,lev,ierr)
    if ( ierr .ne. 0) return
!
!--- set file map
    fieldmapsize=(iend-ista+1)*(jend-jsta+1)
    allocate(fieldmap(fieldmapsize) )
    call set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,ierr,1)
    if(ierr.ne.0) return
!---
    idispstt=int(jrec-1,8)*int(gfile%fieldsize*8+8,8)
    call readmpi8(gfile,fieldmapsize,fieldmap,data,ierr,idispstt)
    if(ierr.ne.0) return
    deallocate(fieldmap)
!
    iret=0
!
    return
  end subroutine nemsio_readrecvbind8
!------------------------------------------------------------------------------
  subroutine nemsio_searchrecv(gfile,jrec,name,levtyp,lev,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: search rec number giving rec name, levtyp and lev
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                 :: gfile
    integer(nemsio_intkind),intent(out)           :: jrec
    character(*),intent(in)                      :: name, levtyp
    integer(nemsio_intkind),intent(in)            :: lev
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer i, nsize

    iret=-11
    jrec=0
    do i=1,gfile%nrec
      if ( trim(name) .eq. trim(gfile%recname(i)) .and.  &
        trim(levtyp) .eq. trim(gfile%reclevtyp(i)) .and.  &
        lev .eq. gfile%reclev(i) ) then
           jrec=i
           exit
      endif
    enddo
    if ( jrec .ne.0 ) iret=0
!
    return
  end subroutine nemsio_searchrecv
!------------------------------------------------------------------------------
!
!*****************  no read grb1 data set :  **********************************
!
!------------------------------------------------------------------------------
!##############################################################################
!
!*****************   write data set :  ********************************
!
!##############################################################################
!------------------------------------------------------------------------------

!*****************   write out bin data set :  ********************************

!------------------------------------------------------------------------------
  subroutine nemsio_writerecbin4(gfile,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_realkind),intent(in)              :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!
    real(nemsio_dblekind),allocatable   :: data8(:)
!
    if(trim(gfile%gdatatype)=='bin4') then
      call nemsio_writerecbind4(gfile,ista,iend,jsta,jend,jrec,data,iret)
    elseif( trim(gfile%gdatatype)=='bin8') then
      allocate(data8(size(data)))
      data8=data
      call nemsio_writerecbind8(gfile,ista,iend,jsta,jend,jrec,data8,iret)
      deallocate(data8)
    endif
!
    return
  end subroutine nemsio_writerecbin4
!
!------------------------------------------------------------------------------
  subroutine nemsio_writerecbin8(gfile,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_dblekind),intent(in)              :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!
    real(nemsio_realkind),allocatable  :: data4(:)
!
    if(trim(gfile%gdatatype)=='bin4') then
      allocate(data4(size(data)))
      data4=data
      call nemsio_writerecbind4(gfile,ista,iend,jsta,jend,jrec,data4,iret)
      deallocate(data4)
    elseif( trim(gfile%gdatatype)=='bin8') then
      call nemsio_writerecbind8(gfile,ista,iend,jsta,jend,jrec,data,iret)
    endif
!
    return
  end subroutine nemsio_writerecbin8
!
!------------------------------------------------------------------------------
  subroutine nemsio_writerecvbin4(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: write out data (bin4) by field
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    character(*),intent(in)                      :: name, levtyp
    integer(nemsio_intkind),intent(in)            :: lev
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    real(nemsio_realkind),intent(in)              :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer :: jrec, ierr
!
    real(nemsio_dblekind),allocatable   :: data8(:)
!
    if(trim(gfile%gdatatype)=='bin4') then
      call nemsio_writerecvbind4(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
    elseif( trim(gfile%gdatatype)=='bin8') then
      allocate(data8(size(data)))
      data8=data
      call nemsio_writerecvbind8(gfile,name,levtyp,lev,ista,iend,jsta,jend,data8,iret)
      deallocate(data8)
    endif
!
    return
  end subroutine nemsio_writerecvbin4
!
!------------------------------------------------------------------------------
  subroutine nemsio_writerecvbin8(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    character(*),intent(in)                      :: name, levtyp
    integer(nemsio_intkind),intent(in)            :: lev
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    real(nemsio_dblekind),intent(in)              :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer :: jrec, ierr
!
    real(nemsio_realkind),allocatable  :: data4(:)
!
    if(trim(gfile%gdatatype)=='bin4') then
      allocate(data4(size(data)))
      data4=data
      call nemsio_writerecvbind4(gfile,name,levtyp,lev,ista,iend,jsta,jend,data4,iret)
      deallocate(data4)
    elseif( trim(gfile%gdatatype)=='bin8') then
      call nemsio_writerecvbind8(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
    endif
!
    return
  end subroutine nemsio_writerecvbin8
!
!------------------------------------------------------------------------------
  subroutine nemsio_writerecbind4(gfile,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_realkind),intent(in)              :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!local vars
    integer,allocatable                 :: fieldmap(:)
    integer fieldmapsize,fieldmapsize1
    integer ierr
    real(nemsio_realkind),allocatable   :: datatmp(:)
    integer(8) idispstt
!
    iret=-11
!
!--- check lead_task and last_task
    if(gfile%mype==gfile%lead_task .and. ista/=1 .and. jsta/=1 ) &
       call nemsio_stop("lead_task\'s subdomain must cover the (1,1) corner") 
!
!--- set file map
    if(gfile%mype==gfile%lead_task) then
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)+2
      fieldmapsize1=(iend-ista+1)*(jend-jsta+1)
    else
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)
    endif
    allocate(datatmp(fieldmapsize))
    allocate(fieldmap(fieldmapsize) )
    call set_mpimap_wrt(gfile,ista,iend,jsta,jend,fieldmap,ierr,1)
!    print *,'in writerecbin4, after set_mpimap,ierr=',ierr,fieldmap(1:5),  &
!        fieldmap(fieldmapsize-4:fieldmapsize)
    if(ierr.ne.0) return
!
!--- prepare data
    if(gfile%mype.eq.gfile%lead_task) then
      datatmp(1)=gfile%fieldsize_real4
      datatmp(fieldmapsize)=datatmp(1)
      datatmp(2:fieldmapsize-1)=data(1:fieldmapsize1)
!      print *,'datatmp=',datatmp(1),datatmp(fieldmapsize),maxval(datatmp(2:fieldmapsize-1)), &
!        minval(datatmp(2:fieldmapsize-1))
    else
      datatmp(:)=data(:)
!      print *,'datatmp=',maxval(datatmp), minval(datatmp)
    endif
!
!--- write out data into files
    idispstt=int(jrec-1,8)*int(gfile%fieldsize*4+8,8)
    call writempi4(gfile,fieldmapsize,fieldmap,datatmp,ierr,idispstt)
!    print *,'in writempi4, after writempi4,ierr=',ierr
    if(ierr.ne.0) return
    deallocate(fieldmap)
!
    call mpi_barrier(gfile%mpi_comm,ierr)
    if(ierr.ne.0) return
!    print *,'in writempi4, end, after mpi_barrier,ierr=',ierr

    iret=0

    return
  end subroutine nemsio_writerecbind4
!------------------------------------------------------------------------------
  subroutine nemsio_writerecbind8(gfile,ista,iend,jsta,jend,jrec,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    integer(nemsio_intkind),intent(in)            :: jrec
    real(nemsio_dblekind),intent(in)              :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
!local vars
    integer,allocatable                :: fieldmap(:)
    real(nemsio_dblekind),allocatable  :: datatmp(:)
    integer fieldmapsize
    integer ierr
    integer(8) idispstt
!
    iret=-11
!
!--- check lead_task and last_task
    if(gfile%mype==gfile%lead_task .and. ista/=1 .and. jsta/=1 ) &
       call nemsio_stop("lead_task\'s subdomain must cover the (1,1) corner")
!
!--- set file map
    if(gfile%mype==gfile%lead_task) then
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)+2
    else
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)
    endif
    allocate(datatmp(fieldmapsize))
    allocate(fieldmap(fieldmapsize) )
    call set_mpimap_wrt(gfile,ista,iend,jsta,jend,fieldmap,ierr,1)
!    print *,'in writerecbin4, after set_mpimap,ierr=',ierr
    if(ierr.ne.0) return
!
!--- prepare data
    if(gfile%mype.eq.gfile%lead_task) then
      datatmp(1)=gfile%fieldsize_real8
      datatmp(fieldmapsize)=datatmp(1)
      datatmp(2:fieldmapsize-1)=data(:)
    else
      datatmp(:)=data(:)
    endif
!
!---
    idispstt=int(jrec-1,8)*int(gfile%fieldsize*8+8,8)
    call writempi8(gfile,fieldmapsize,fieldmap,datatmp,ierr,idispstt)
    if(ierr.ne.0) return
    deallocate(fieldmap)
!
    iret=0
!
    return
  end subroutine nemsio_writerecbind8
!------------------------------------------------------------------------------
  subroutine nemsio_writerecvbind4(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    character(*),intent(in)                       :: name, levtyp
    integer(nemsio_intkind),intent(in)            :: lev
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    real(nemsio_realkind),intent(in)              :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer :: jrec, ierr
!local vars
    integer,allocatable                :: fieldmap(:)
    real(nemsio_realkind),allocatable  :: datatmp(:)
    integer fieldmapsize
    integer(8) idispstt
!
    iret=-11
!
    call nemsio_searchrecv(gfile,jrec,name,levtyp,lev,ierr)
    if(ierr.ne.0) return
!
!--- check lead_task and last_task
    if(gfile%mype==gfile%lead_task .and. ista/=1 .and. jsta/=1 ) &
       call nemsio_stop("lead_task\'s subdomain must cover the (1,1) corner")
!
!--- set file map
    if(gfile%mype==gfile%lead_task) then
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)+2
    else
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)
    endif
    allocate(datatmp(fieldmapsize))
    allocate(fieldmap(fieldmapsize) )
    call set_mpimap_wrt(gfile,ista,iend,jsta,jend,fieldmap,ierr,1)
!    print *,'in writerecbin4, after set_mpimap,ierr=',ierr
    if(ierr.ne.0) return
!
!--- prepare data
    if(gfile%mype.eq.gfile%lead_task) then
      datatmp(1)=gfile%fieldsize_real4
      datatmp(fieldmapsize)=datatmp(1)
      datatmp(2:fieldmapsize-1)=data(:)
    else
      datatmp(:)=data(:)
    endif
!---
    idispstt=int(jrec-1,8)*int(gfile%fieldsize*4+8,8)
    call writempi4(gfile,fieldmapsize,fieldmap,datatmp,iret,idispstt=idispstt)
    if(iret.ne.0) return
    deallocate(fieldmap)
!
    iret=0
!
    return
  end subroutine nemsio_writerecvbind4
!
!------------------------------------------------------------------------------
  subroutine nemsio_writerecvbind8(gfile,name,levtyp,lev,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read nemsio data (bin) by record number into a 2D 32 bits array
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)              :: gfile
    character(*),intent(in)                      :: name, levtyp
    integer(nemsio_intkind),intent(in)            :: lev
    integer(nemsio_intkind),intent(in)            :: ista,iend,jsta,jend
    real(nemsio_dblekind),intent(in)              :: data(:)
    integer(nemsio_intkind),optional,intent(out)  :: iret
    integer :: jrec, ierr
!local vars
    integer,allocatable               :: fieldmap(:)
    real(nemsio_dblekind),allocatable :: datatmp(:)
    integer fieldmapsize
    integer(8) idispstt
!
    iret=-11
!
    call nemsio_searchrecv(gfile,jrec,name,levtyp,lev,ierr)
!
!--- check lead_task and last_task
    if(gfile%mype==gfile%lead_task .and. ista/=1 .and. jsta/=1 ) &
       call nemsio_stop("lead_task\'s subdomain must cover the (1,1) corner")
!
!--- set file map
    if(gfile%mype==gfile%lead_task) then
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)+2
    else
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)
    endif
    allocate(datatmp(fieldmapsize))
    allocate(fieldmap(fieldmapsize) )
    call set_mpimap_wrt(gfile,ista,iend,jsta,jend,fieldmap,ierr,1)
!    print *,'in writerecbin4, after set_mpimap,ierr=',ierr
    if(ierr.ne.0) return
!
!--- prepare data
    if(gfile%mype.eq.gfile%lead_task) then
      datatmp(1)=gfile%fieldsize_real8
      datatmp(fieldmapsize)=datatmp(1)
      datatmp(2:fieldmapsize-1)=data(:)
    else
      datatmp(:)=data(:)
    endif
!
!--- set file map
!     fieldmapsize=(iend-ista+1)*(jend-jsta+1)
!     allocate(fieldmap(fieldmapsize) )
!     call set_mpimap_wrt(gfile,ista,iend,jsta,jend,fieldmap,iret,1)
!     if(iret.ne.0) return
!---
     idispstt=int(jrec-1,8)*int(gfile%fieldsize*8+8,8)
     call writempi8(gfile,fieldmapsize,fieldmap,data,iret,idispstt=idispstt)
     if(iret.ne.0) return
     deallocate(fieldmap)
!
    iret=0
!
    return
  end subroutine nemsio_writerecvbind8
!------------------------------------------------------------------------------
!
!***************** no write out grb data set :  ********************************
!
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
  subroutine nemsio_chkgfary(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: check if arrays in gfile is allocated and with right size
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)         :: gfile
    integer(nemsio_intkind),intent(out)   :: iret
    integer   :: ios
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    iret=-2
    if ( gfile%dimx .eq. nemsio_intfill .or. gfile%dimy .eq. nemsio_intfill &
        .or. gfile%dimz .eq. nemsio_intfill .or. gfile%nrec .eq. nemsio_intfill &
        .or. gfile%idate(1) .eq.nemsio_intfill .or. gfile%ntrac .eq.nemsio_intfill ) then
        return
    endif
    if (.not. allocated(gfile%vcoord) .or. size(gfile%vcoord).ne. &
       (gfile%dimz+1)*3*2 ) then
       call nemsio_almeta1(gfile,ios)
       if (ios .ne. 0) return
    endif
    if (.not.allocated(gfile%lat) .or. size(gfile%lat).ne.gfile%fieldsize .or.&
        .not.allocated(gfile%lon) .or. size(gfile%lon).ne.gfile%fieldsize .or.&
        .not.allocated(gfile%dx) .or. size(gfile%dx).ne.gfile%fieldsize .or.&
        .not.allocated(gfile%dy) .or. size(gfile%dy).ne.gfile%fieldsize) then
        call nemsio_almeta2(gfile,ios)
        if (ios .ne. 0) return
    endif
    if (.not.allocated(gfile%Cpi) .or. size(gfile%Cpi).ne.gfile%ntrac+1 .or. &
        .not.allocated(gfile%Ri) .or. size(gfile%Ri).ne.gfile%ntrac+1 ) then
        call nemsio_almeta3(gfile,ios)
        if (ios .ne. 0) return
    endif

    if (allocated(gfile%recname) .and. size(gfile%recname).eq.gfile%nrec)&
    then
        if (allocated(gfile%reclevtyp) .and. size(gfile%reclevtyp) &
        .eq.gfile%nrec) then
           if (allocated(gfile%reclev) .and. size(gfile%reclev).eq. &
             gfile%nrec) then
               iret=0
               return
           endif
         endif
   endif
   call  nemsio_almeta4(gfile,ios)
   if (ios .ne. 0) return
   iret=0
  end subroutine nemsio_chkgfary
!------------------------------------------------------------------------------
  subroutine nemsio_almeta(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: allocate all the arrays in gfile
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)  :: gfile 
    integer(nemsio_intkind),intent(out)  :: iret
    integer ::dimvcoord1,dimvcoord2,dimnmmlev
    integer ::dimrecname,dimreclevtyp,dimreclev
    integer ::dimfield
    integer ::dimcpr
    integer ::iret1,iret2,iret3,iret4
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dimvcoord1=gfile%dimz+1
    dimrecname=gfile%nrec
    dimreclevtyp=gfile%nrec
    dimreclev=gfile%nrec
    dimfield=gfile%fieldsize
    dimcpr=gfile%ntrac+1
    if(allocated(gfile%recname)) deallocate(gfile%recname)
    if(allocated(gfile%reclevtyp)) deallocate(gfile%reclevtyp)
    if(allocated(gfile%reclev)) deallocate(gfile%reclev)
    if(allocated(gfile%vcoord)) deallocate(gfile%vcoord)
    if(allocated(gfile%lat)) deallocate(gfile%lat)
    if(allocated(gfile%lon)) deallocate(gfile%lon)
    if(allocated(gfile%dx)) deallocate(gfile%dx)
    if(allocated(gfile%dy)) deallocate(gfile%dy)
    if(allocated(gfile%Cpi)) deallocate(gfile%Cpi)
    if(allocated(gfile%Ri)) deallocate(gfile%Ri)
    allocate(gfile%recname(dimrecname),  gfile%reclevtyp(dimreclevtyp), &
             gfile%reclev(dimreclev), &
             stat=iret1)
    allocate(gfile%vcoord(dimvcoord1,3,2) ,stat=iret2) 
    allocate(gfile%lat(dimfield), gfile%lon(dimfield), &
             gfile%dx(dimfield), gfile%dy(dimfield) ,stat=iret3)
    allocate(gfile%Cpi(dimcpr), gfile%Ri(dimcpr), stat=iret4)

!    print *,'iret1=',iret1,'iret2=',iret2,'dimx=',gfile%dimx,'dimy=',gfile%dimy,'nframe=',gfile%nframe
    iret=abs(iret1)+abs(iret2)+abs(iret3)+abs(iret4)
    if(iret.eq.0) then
      gfile%reclev=nemsio_intfill
      gfile%recname=' '
      gfile%reclevtyp=' '
      gfile%vcoord=nemsio_realfill
      gfile%lat=nemsio_realfill
      gfile%lon=nemsio_realfill
      gfile%dx=nemsio_realfill
      gfile%dy=nemsio_realfill
      gfile%Cpi=nemsio_realfill
      gfile%Ri=nemsio_realfill
    endif
    if(iret.ne.0) iret=-6
  end subroutine nemsio_almeta
!------------------------------------------------------------------------------
  subroutine nemsio_alextrameta(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: allocate all the arrays in gfile
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)  :: gfile
    integer(nemsio_intkind),intent(out)  :: iret
    integer ::iret1,iret2,iret3,iret4
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    iret=-6
    if(gfile%extrameta) then
!      print *,'nmetavari=',gfile%nmetavari,'nmetavarr=',gfile%nmetavarr, &
!              'nmetavarl=',gfile%nmetavarl,'nmetavarc=',gfile%nmetavarc, &
!              'nmetaaryi=',gfile%nmetaaryi,'nmetaaryr=',gfile%nmetaaryi, &
!              'nmetaaryl=',gfile%nmetaaryl,'nmetaaryc=',gfile%nmetaaryc
      if(gfile%nmetavari.gt.0) then
         if(allocated(gfile%variname)) deallocate(gfile%variname)
         if(allocated(gfile%varival)) deallocate(gfile%varival)
         allocate(gfile%variname(gfile%nmetavari), &
                  gfile%varival(gfile%nmetavari), stat=iret1 )
         if(iret1.ne.0) return
      endif
      if(gfile%nmetavarr.gt.0) then
         if(allocated(gfile%varrname)) deallocate(gfile%varrname)
         if(allocated(gfile%varrval)) deallocate(gfile%varrval)
         allocate(gfile%varrname(gfile%nmetavarr), &
                  gfile%varrval(gfile%nmetavarr), stat=iret1 )
         if(iret1.ne.0) return
      endif
      if(gfile%nmetavarl.gt.0) then
         if(allocated(gfile%varlname)) deallocate(gfile%varlname)
         if(allocated(gfile%varlval)) deallocate(gfile%varlval)
         allocate(gfile%varlname(gfile%nmetavarl), &
                  gfile%varlval(gfile%nmetavarl), stat=iret1 )
         if(iret1.ne.0) return
      endif
      if(gfile%nmetavarc.gt.0) then
         if(allocated(gfile%varcname)) deallocate(gfile%varcname)
         if(allocated(gfile%varcval)) deallocate(gfile%varcval)
         allocate(gfile%varcname(gfile%nmetavarc), &
                  gfile%varcval(gfile%nmetavarc), stat=iret1 )
         if(iret1.ne.0) return
      endif
      if(gfile%nmetaaryi.gt.0) then
         if(allocated(gfile%aryiname)) deallocate(gfile%aryiname)
         if(allocated(gfile%aryilen)) deallocate(gfile%aryilen)
         if(allocated(gfile%aryival)) deallocate(gfile%aryival)
         allocate(gfile%aryiname(gfile%nmetaaryi), &
                  gfile%aryilen(gfile%nmetaaryi), stat=iret1 )
         if(iret1.ne.0) return
      endif
      if(gfile%nmetaaryr.gt.0) then
         if(allocated(gfile%aryrname)) deallocate(gfile%aryrname)
         if(allocated(gfile%aryrlen)) deallocate(gfile%aryrlen)
         if(allocated(gfile%aryrval)) deallocate(gfile%aryrval)
         allocate(gfile%aryrname(gfile%nmetaaryr), &
                  gfile%aryrlen(gfile%nmetaaryr), stat=iret1 )
         if(iret1.ne.0) return
      endif
      if(gfile%nmetaaryl.gt.0) then
         if(allocated(gfile%arylname)) deallocate(gfile%arylname)
         if(allocated(gfile%aryllen)) deallocate(gfile%aryllen)
         if(allocated(gfile%arylval)) deallocate(gfile%arylval)
         allocate(gfile%arylname(gfile%nmetaaryl), &
                  gfile%aryllen(gfile%nmetaaryl), stat=iret1 )
         if(iret1.ne.0) return
      endif
      if(gfile%nmetaaryc.gt.0) then
         if(allocated(gfile%arycname)) deallocate(gfile%arycname)
         if(allocated(gfile%aryclen)) deallocate(gfile%aryclen)
         if(allocated(gfile%arycval)) deallocate(gfile%arycval)
         allocate(gfile%arycname(gfile%nmetaaryc), &
                  gfile%aryclen(gfile%nmetaaryc), stat=iret1 )
         if(iret1.ne.0) return
      endif
    endif

    iret=0
!    print *,'end of alextrameta'
  end subroutine nemsio_alextrameta
!------------------------------------------------------------------------------
  subroutine nemsio_almeta1(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: allocate vcoord in gfile
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)  :: gfile
    integer(nemsio_intkind),intent(out)  :: iret
    integer :: dimvcoord1,dimnmmlev,dimnmmnsoil
    integer :: dimgsilev
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dimvcoord1=gfile%dimz+1
    if(allocated(gfile%vcoord)) deallocate(gfile%vcoord)
    allocate(gfile%vcoord(dimvcoord1,3,2), stat=iret)
    if(iret.eq.0) then
      gfile%vcoord=nemsio_realfill
    endif
    if(iret.ne.0) iret=-6
  end subroutine nemsio_almeta1
!------------------------------------------------------------------------------
  subroutine nemsio_almeta2(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: allocate lat1d in gfile
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)  :: gfile
    integer(nemsio_intkind),intent(out)  :: iret
    integer :: dimlat
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dimlat=gfile%fieldsize
    if(allocated(gfile%lat)) deallocate(gfile%lat)
    if(allocated(gfile%lon)) deallocate(gfile%lon)
    if(allocated(gfile%dx)) deallocate(gfile%dx)
    if(allocated(gfile%dy)) deallocate(gfile%dy)
    allocate(gfile%lat(dimlat),gfile%lon(dimlat), &
             gfile%dx(dimlat),gfile%dy(dimlat), stat=iret)
    if(iret.eq.0) then
      gfile%lat=nemsio_realfill
      gfile%lon=nemsio_realfill
      gfile%dx=nemsio_realfill
      gfile%dy=nemsio_realfill
    endif
    if(iret.ne.0) iret=-6
  end subroutine nemsio_almeta2
!------------------------------------------------------------------------------
  subroutine nemsio_almeta3(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: allocate lon1d in gfile
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)  :: gfile
    integer(nemsio_intkind),intent(out)  :: iret
    integer :: dim1d
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dim1d=gfile%ntrac+1
    if(allocated(gfile%Cpi)) deallocate(gfile%Cpi)
    if(allocated(gfile%Ri)) deallocate(gfile%Ri)
    allocate(gfile%Cpi(dim1d),gfile%Ri(dim1d),stat=iret)
    if(iret.eq.0) then
       gfile%Cpi=nemsio_realfill
       gfile%Ri=nemsio_realfill
    endif
    if(iret.ne.0) iret=-6
  end subroutine nemsio_almeta3
!------------------------------------------------------------------------------
  subroutine nemsio_almeta4(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: allocate recnam, reclvevtyp, and reclev in gfile
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)  :: gfile
    integer(nemsio_intkind),intent(out)  :: iret
    integer :: dimrecname,dimreclevtyp,dimreclev
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    dimrecname=gfile%nrec
    dimreclevtyp=gfile%nrec
    dimreclev=gfile%nrec
    if(allocated(gfile%recname)) deallocate(gfile%recname)
    if(allocated(gfile%reclevtyp)) deallocate(gfile%reclevtyp)
    if(allocated(gfile%reclev)) deallocate(gfile%reclev)
    allocate(gfile%recname(dimrecname),  gfile%reclevtyp(dimreclevtyp), &
             gfile%reclev(dimreclev), stat=iret)
    if(iret.eq.0) then
      gfile%reclev=nemsio_intfill
      gfile%recname=' '
      gfile%reclevtyp=' '
    endif
    if(iret.ne.0) iret=-6
  end subroutine nemsio_almeta4
!------------------------------------------------------------------------------
  subroutine nemsio_axmeta(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: empty gfile variables and decallocate arrays in gfile
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)      :: gfile
    integer(nemsio_intkind),intent(out)  :: iret
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    iret=-6
    gfile%gtype=' '
    gfile%gdatatype=' '
    gfile%modelname=' '
    gfile%version=nemsio_intfill
    gfile%nmeta=nemsio_intfill
    gfile%lmeta=nemsio_intfill
    gfile%nrec=nemsio_intfill
    gfile%idate(1:7)=nemsio_intfill
    gfile%nfday=nemsio_intfill
    gfile%nfhour=nemsio_intfill
    gfile%nfminute=nemsio_intfill
    gfile%nfsecondn=nemsio_intfill
    gfile%nfsecondd=nemsio_intfill
    gfile%dimx=nemsio_intfill
    gfile%dimy=nemsio_intfill
    gfile%dimz=nemsio_intfill
    gfile%nframe=nemsio_intfill
    gfile%nsoil=nemsio_intfill
    gfile%ntrac=nemsio_intfill
    gfile%jcap=nemsio_intfill
    gfile%ncldt=nemsio_intfill
    gfile%idvc=nemsio_intfill
    gfile%idsl=nemsio_intfill
    gfile%idvm=nemsio_intfill
    gfile%idrt=nemsio_intfill
    gfile%rlon_min=nemsio_realfill
    gfile%rlon_max=nemsio_realfill
    gfile%rlat_min=nemsio_realfill
    gfile%rlat_max=nemsio_realfill
    gfile%extrameta=nemsio_logicfill
    gfile%nmetavari=nemsio_intfill
    gfile%nmetavarr=nemsio_intfill
    gfile%nmetavarl=nemsio_intfill
    gfile%nmetavarc=nemsio_intfill
    gfile%nmetaaryi=nemsio_intfill
    gfile%nmetaaryr=nemsio_intfill
    gfile%nmetaaryl=nemsio_intfill
    gfile%nmetaaryc=nemsio_intfill

    if(allocated(gfile%recname)) deallocate(gfile%recname)
    if(allocated(gfile%reclevtyp)) deallocate(gfile%reclevtyp)
    if(allocated(gfile%reclev)) deallocate(gfile%reclev)
    if(allocated(gfile%vcoord)) deallocate(gfile%vcoord)
    if(allocated(gfile%lat)) deallocate(gfile%lat)
    if(allocated(gfile%lon)) deallocate(gfile%lon)
    if(allocated(gfile%dx)) deallocate(gfile%dx)
    if(allocated(gfile%dy)) deallocate(gfile%dy)
    if(allocated(gfile%Cpi)) deallocate(gfile%Cpi)
    if(allocated(gfile%Ri)) deallocate(gfile%Ri)
!
    gfile%mbuf=0
    gfile%nnum=0
    gfile%nlen=0
    gfile%mnum=0
    if(allocated(gfile%cbuf)) deallocate(gfile%cbuf)
    if(allocated(gfile%headvariname)) deallocate(gfile%headvariname)
    if(allocated(gfile%headvarrname)) deallocate(gfile%headvarrname)
    if(allocated(gfile%headvarlname)) deallocate(gfile%headvarlname)
    if(allocated(gfile%headvarcname)) deallocate(gfile%headvarcname)
    if(allocated(gfile%headvarival)) deallocate(gfile%headvarival)
    if(allocated(gfile%headvarrval)) deallocate(gfile%headvarrval)
    if(allocated(gfile%headvarlval)) deallocate(gfile%headvarlval)
    if(allocated(gfile%headvarcval)) deallocate(gfile%headvarcval)
    if(allocated(gfile%headaryiname)) deallocate(gfile%headaryiname)
    if(allocated(gfile%headaryrname)) deallocate(gfile%headaryrname)
    if(allocated(gfile%headarycname)) deallocate(gfile%headarycname)
    if(allocated(gfile%headaryival)) deallocate(gfile%headaryival)
    if(allocated(gfile%headaryrval)) deallocate(gfile%headaryrval)
    if(allocated(gfile%headarycval)) deallocate(gfile%headarycval)
    iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine nemsio_axmeta
!------------------------------------------------------------------------------
  subroutine nemsio_setfhead(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: required file header (default)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)     :: gfile
    integer(nemsio_intkind),intent(out)  :: iret
    integer(nemsio_intkind) i,j,k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    iret=-17
    gfile%headvarinum=29
    gfile%headvarrnum=4
    gfile%headvarlnum=1
    gfile%headvarcnum=3
    gfile%headaryinum=2
    gfile%headaryrnum=7
    gfile%headarycnum=2
!
    allocate(gfile%headvariname(gfile%headvarinum),gfile%headvarival(gfile%headvarinum) )
    gfile%headvariname(1)='version'
    gfile%headvarival(1)=gfile%version
    gfile%headvariname(2)='nmeta'
    gfile%headvarival(2)=gfile%nmeta
    gfile%headvariname(3)='lmeta'
    gfile%headvarival(3)=gfile%lmeta
    gfile%headvariname(4)='nrec'
    gfile%headvarival(4)=gfile%nrec
    gfile%headvariname(5)='nfday'
    gfile%headvarival(5)=gfile%nfday
    gfile%headvariname(6)='nfhour'
    gfile%headvarival(6)=gfile%nfhour
    gfile%headvariname(7)='nfminute'
    gfile%headvarival(7)=gfile%nfminute
    gfile%headvariname(8)='nfsecondn'
    gfile%headvarival(8)=gfile%nfsecondn
    gfile%headvariname(9)='nfsecondd'
    gfile%headvarival(9)=gfile%nfsecondd
    gfile%headvariname(10)='dimx'
    gfile%headvarival(10)=gfile%dimx
    gfile%headvariname(11)='dimy'
    gfile%headvarival(11)=gfile%dimy
    gfile%headvariname(12)='dimz'
    gfile%headvarival(12)=gfile%dimz
    gfile%headvariname(13)='nframe'
    gfile%headvarival(13)=gfile%nframe
    gfile%headvariname(14)='nsoil'
    gfile%headvarival(14)=gfile%nsoil
    gfile%headvariname(15)='ntrac'
    gfile%headvarival(15)=gfile%ntrac
    gfile%headvariname(16)='jcap'
    gfile%headvarival(16)=gfile%jcap
    gfile%headvariname(17)='ncldt'
    gfile%headvarival(17)=gfile%ncldt
    gfile%headvariname(18)='idvc'
    gfile%headvarival(18)=gfile%idvc
    gfile%headvariname(19)='idsl'
    gfile%headvarival(19)=gfile%idsl
    gfile%headvariname(20)='idvm'
    gfile%headvarival(20)=gfile%idvm
    gfile%headvariname(21)='idrt'
    gfile%headvarival(21)=gfile%idrt
    gfile%headvariname(22)='nmetavari'
    gfile%headvarival(22)=gfile%nmetavari
    gfile%headvariname(23)='nmetavarr'
    gfile%headvarival(23)=gfile%nmetavarr
    gfile%headvariname(24)='nmetavarl'
    gfile%headvarival(24)=gfile%nmetavarl
    gfile%headvariname(25)='nmetavarc'
    gfile%headvarival(25)=gfile%nmetavarc
    gfile%headvariname(26)='nmetaaryi'
    gfile%headvarival(26)=gfile%nmetaaryi
    gfile%headvariname(27)='nmetaaryr'
    gfile%headvarival(27)=gfile%nmetaaryr
    gfile%headvariname(28)='nmetaaryl'
    gfile%headvarival(28)=gfile%nmetaaryl
    gfile%headvariname(29)='nmetaaryc'
    gfile%headvarival(29)=gfile%nmetaaryc
!
    allocate(gfile%headvarrname(gfile%headvarrnum),gfile%headvarrval(gfile%headvarrnum) )
    gfile%headvarrname(1)='rlon_min'
    gfile%headvarrval(1)=gfile%rlon_min
    gfile%headvarrname(2)='rlon_max'
    gfile%headvarrval(2)=gfile%rlon_max
    gfile%headvarrname(3)='rlat_min'
    gfile%headvarrval(3)=gfile%rlat_min
    gfile%headvarrname(4)='rlat_min'
    gfile%headvarrval(4)=gfile%rlat_min
!
    allocate(gfile%headvarcname(gfile%headvarcnum),gfile%headvarcval(gfile%headvarcnum) )
    gfile%headvarcname(1)='gtype'
    gfile%headvarcval(1)=gfile%gtype
    gfile%headvarcname(2)='modelname'
    gfile%headvarcval(2)=gfile%modelname
    gfile%headvarcname(3)='gdatatype'
    gfile%headvarcval(3)=gfile%gdatatype
!head logic var
!    write(0,*)'before setfhead, headvarl,nrec=',gfile%nrec 
    allocate(gfile%headvarlname(gfile%headvarlnum),gfile%headvarlval(gfile%headvarlnum) )
    gfile%headvarlname(1)='extrameta'
    gfile%headvarlval(1)=gfile%extrameta
!
!--- gfile%head int ary
!    write(0,*)'before setfhead, headaryi,nrec=',gfile%nrec,gfile%headaryinum
    allocate(gfile%headaryiname(gfile%headaryinum) )
    allocate(gfile%headaryival(max(size(gfile%reclev),7),gfile%headaryinum))
    gfile%headaryiname(1)='idate'
    gfile%headaryival(1:7,1)=gfile%idate(1:7)
    gfile%headaryiname(2)='reclev'
    gfile%headaryival(:,2)=gfile%reclev(:)
!
!--- gfile%head real ary
!    write(0,*)'before setfhead, headaryr,',gfile%headaryrnum ,gfile%fieldsize
    allocate(gfile%headaryrname(gfile%headaryrnum) )
    allocate(gfile%headaryrval(max(gfile%fieldsize,(gfile%dimz+1)*6),gfile%headaryrnum))
    gfile%headaryrname(1)='vcoord'
!    print *,'in setfhead, headaryr, before gfile%headaryrval 1',gfile%dimz
    do j=1,2
     do i=1,3
      do k=1,gfile%dimz+1
       gfile%headaryrval(k+((j-1)*3+i-1)*(gfile%dimz+1),1)=gfile%vcoord(k,i,j)
      enddo
     enddo
    enddo
    gfile%headaryrname(2)='lat'
    gfile%headaryrval(:,2)=gfile%lat
    gfile%headaryrname(3)='lon'
    gfile%headaryrval(:,3)=gfile%lon
    gfile%headaryrname(4)='dx'
    gfile%headaryrval(:,4)=gfile%dx
    gfile%headaryrname(5)='dy'
    gfile%headaryrval(:,5)=gfile%dy
    gfile%headaryrname(6)='cpi'
    gfile%headaryrval(1:size(gfile%cpi),6)=gfile%cpi(:)
    gfile%headaryrname(7)='ri'
    gfile%headaryrval(1:size(gfile%ri),7)=gfile%ri(:)
!
!--- gfile%head char var
!    write(0,*)'before setfhead, headaryc,nrec=',gfile%nrec,gfile%headarycnum
    allocate(gfile%headarycname(gfile%headarycnum) )
    allocate(gfile%headarycval(size(gfile%recname),gfile%headarycnum))
    gfile%headarycname(1)='recname'
    gfile%headarycval(:,1)=gfile%recname
    gfile%headarycname(2)='reclevtyp'
    gfile%headarycval(:,2)=gfile%reclevtyp
!
!    write(0,*)'end ef nemsio_setfhead'
    iret=0
  end subroutine nemsio_setfhead
!------------------------------------------------------------------------------
  subroutine nemsio_getrechead(gfile,jrec,name,levtyp,lev,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: given record number, return users record name, lev typ, and levs
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(in)                :: gfile
    integer(nemsio_intkind),intent(in)           :: jrec
    character*(*),intent(out)                   :: name,levtyp
    integer(nemsio_intkind),intent(out)          :: lev
    integer(nemsio_intkind),optional,intent(out) :: iret
    integer :: ios
! - - - - - - - - - - - - - -  - - - - - - - -  - - - - - - - - - - - - - - - -
    if( present(iret)) iret=-6
    if ( jrec.gt.0 .or. jrec.le.gfile%nrec) then
      name=gfile%recname(jrec)
      levtyp=gfile%reclevtyp(jrec)
      lev=gfile%reclev(jrec)
      if(present(iret)) iret=0
!      print *,'in getrechead, nrec=',gfile%nrec,'name=',name,'levtyp=',levtyp,'lev=',lev
      return
    else
      if ( present(iret))  then
       return
      else
        call nemsio_stop
      endif
    endif
  end subroutine nemsio_getrechead
!------------------------------------------------------------------------------
  subroutine nemsio_gfinit(gfile,iret,recname,reclevtyp,reclev)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: set gfile variables to operational model output
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
    implicit none
    type(nemsio_gfile),intent(inout)     :: gfile
    integer(nemsio_intkind),intent(out)  :: iret
    character(nemsio_charkind),optional,intent(in)  :: recname(:)
    character(nemsio_charkind*2),optional,intent(in):: reclevtyp(:)
    integer(nemsio_intkind),optional,intent(in)     :: reclev(:)
    integer  :: i,j,rec,rec3dopt
    real(nemsio_dblekind),allocatable :: slat(:),wlat(:)
    real(nemsio_dblekind),allocatable :: dx(:)
    real(nemsio_dblekind)             :: radi
    logical ::linit
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! set operational format
!
    iret=-8
    gfile%version=200809
    gfile%nfday=0
    gfile%nfhour=0
    gfile%nfminute=0
    gfile%nfsecondn=0
    gfile%nfsecondd=100
    gfile%extrameta=.false.
    gfile%nmetavari=0
    gfile%nmetavarr=0
    gfile%nmetavarl=0
    gfile%nmetavarc=0
    gfile%nmetaaryi=0
    gfile%nmetaaryr=0
    gfile%nmetaaryl=0
    gfile%nmetaaryc=0
!
!    print *,'in gfinit, modelname=',gfile%modelname

!
   iret=0
  end subroutine nemsio_gfinit
!------------------------------------------------------------------------------
  subroutine nemsio_stop(message)
    implicit none
    character(*),optional,intent(in) :: message
    integer ::ierr
!---
     if ( present(message) ) print *,'message'
     call mpi_finalize(ierr)
     stop
!
  end subroutine nemsio_stop
!------------------------------------------------------------------------------
!  temporary subroutines for basio file unit
    subroutine nemsio_getlu(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: set unit number to the first number available between 600-699
!           according to unit number array fileunit
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
     type(nemsio_gfile),intent (inout) :: gfile
     integer,intent(out) :: iret
     integer :: i
!
     iret=-10
     do i=600,699
       if ( fileunit(i) .eq. 0 ) then 
         gfile%flunit=i
         fileunit(i)=i
         iret=0
         exit
       endif
     enddo
    end subroutine nemsio_getlu
!------------------------------------------------------------------------------
!  temporary subroutines for free unit number 
    subroutine nemsio_clslu(gfile,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: free unit number array index corresponding to unit number
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
     type(nemsio_gfile),intent (inout) :: gfile
     integer, intent(out) :: iret
     iret=-10
     if ( fileunit(gfile%flunit) .ne. 0 ) then
       fileunit(gfile%flunit)=0
       gfile%flunit=0
       iret=0
     endif
    end subroutine nemsio_clslu
!------------------------------------------------------------------------------
!
    subroutine nemsio_denseread4(gfile,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: free unit number array index corresponding to unit number
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: ista,iend,jsta,jend
     real(nemsio_realkind),intent(out)   :: data(:)
     integer,optional,intent(out)        :: iret
!--- local vars
     integer                  :: status(MPI_STATUS_SIZE)
     integer                  :: fieldmapsize,nfld,nfldloop,mfldrmd
     integer,allocatable      :: fieldmap(:)
     integer ios,i,j,nfldsize,fieldmapsize1,k,nstt,nend
     integer(8) idispstt
     real(nemsio_dblekind),allocatable :: tmp(:)
!---
     iret=-25
!
!--- set nfld
    if(trim(gfile%gdatatype).eq.'bin4') then
        nfldsize=gfile%fieldsize+2
     elseif (trim(gfile%gdatatype).eq.'bin8') then
        nfldsize=gfile%fieldsize+1
     endif
     nfld=min(gfile%nrec,nemsio_maxint/nfldsize)
     nfldloop=(gfile%nrec-1)/nfld+1
     mfldrmd=mod(gfile%nrec,nfld)
!     write(0,*)'in dense read,nfld=',nfld,'nfldloop=',nfldloop, &
!       'mfldrmd=',mfldrmd
!--- set file map
     fieldmapsize=(iend-ista+1)*(jend-jsta+1)*nfld
     allocate(fieldmap(fieldmapsize) )
     call set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,ios)
     if(ios.ne.0) return
!---
     do k=1,nfldloop
!
       if(k<nfldloop.or.mfldrmd==0) then
         nstt=(k-1)*fieldmapsize+1
         nend=k*fieldmapsize
       elseif(mfldrmd/=0) then
         nstt=(k-1)*fieldmapsize+1
         nend=gfile%nrec*(iend-ista+1)*(jend-jsta+1)
         deallocate(fieldmap)
         fieldmapsize=(iend-ista+1)*(jend-jsta+1)*mfldrmd
         allocate(fieldmap(fieldmapsize) )
!          print *,'bf set_mpa_read,size=',fieldmapsize,'mfldrmd=',mfldrmd
         call set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,ios)
         if(ios.ne.0) return
       endif
!       print *,'bf readmpi,k=',k,'nstt=',nstt,'nend=',nend

       if(trim(gfile%gdatatype)=='bin4') then
         idispstt=int(k-1,8)*int(nfld,8)*int(gfile%fieldsize*4+8,8)
!         print *,'right bf readmpi4'
         call readmpi4(gfile,fieldmapsize,fieldmap,data(nstt:nend),ios,idispstt)
!         print *,'af readmpi4,k=',k,'tmp=',maxval(data(nstt:nend)),&
!           minval(data(nstt:nend)),'ios=',ios
       else if (trim(gfile%gdatatype)=='bin8') then
         allocate(tmp(size(data)))
         idispstt=int(k-1,8)*int(nfld,8)*int(gfile%fieldsize*8+8,8)
         call readmpi8(gfile,fieldmapsize,fieldmap,tmp(nstt:nend),ios,idispstt)
         data=tmp
         deallocate(tmp)
       endif
       if(ios.ne.0) return
!
     enddo
     deallocate(fieldmap)
!
     iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine nemsio_denseread4
!------------------------------------------------------------------------------
    subroutine nemsio_denseread8(gfile,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read all the fields out in real 8 MPI
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: ista,iend,jsta,jend
     real(nemsio_dblekind),intent(out)  :: data(:)
     integer,optional,intent(out)        :: iret
!--- local vars
     integer                  :: fieldmapsize
     integer,allocatable      :: fieldmap(:)
     integer ios,i,j,nfldsize,nfld,nfldloop,mfldrmd,k,nstt,nend
     integer(8) idispstt
     real(nemsio_realkind),allocatable :: tmp(:)
!---
     iret=-25
!
!--- set nfld
    if(trim(gfile%gdatatype).eq.'bin4') then
        nfldsize=gfile%fieldsize+2
     elseif (trim(gfile%gdatatype).eq.'bin8') then
        nfldsize=gfile%fieldsize+1
     endif
     nfld=min(gfile%nrec,nemsio_maxint/nfldsize)
     nfldloop=(gfile%nrec-1)/nfld+1
     mfldrmd=mod(gfile%nrec,nfld)
!     write(0,*)'in dense read,nfld=',nfld,'nfldloop=',nfldloop, &
!       'mfldrmd=',mfldrmd
!--- set file map
     fieldmapsize=(iend-ista+1)*(jend-jsta+1)*nfld
     allocate(fieldmap(fieldmapsize) )
     call set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,ios)
     if(ios.ne.0) return
!---
     do k=1,nfldloop
!
       if(k<nfldloop.or.mfldrmd==0) then
         nstt=(k-1)*fieldmapsize+1
         nend=k*fieldmapsize
       elseif(mfldrmd/=0) then
         nstt=(k-1)*fieldmapsize+1
         nend=gfile%nrec*(iend-ista+1)*(jend-jsta+1)
         deallocate(fieldmap)
         fieldmapsize=(iend-ista+1)*(jend-jsta+1)*mfldrmd
         allocate(fieldmap(fieldmapsize) )
         call set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,ios)
         if(ios.ne.0) return
       endif
!
!---
       if(trim(gfile%gdatatype)=='bin4') then
         allocate(tmp(size(data)))
         idispstt=int(k-1,8)*int(nfld,8)*int(gfile%fieldsize*4+8,8)
         call readmpi4(gfile,fieldmapsize,fieldmap,tmp,ios,idispstt)
         data=tmp
         deallocate(tmp)
       elseif(trim(gfile%gdatatype)=='bin8') then
         idispstt=int(k-1,8)*int(nfld,8)*int(gfile%fieldsize*8+8,8)
         call readmpi8(gfile,fieldmapsize,fieldmap,data,ios,idispstt)
       endif
       if(ios.ne.0) return
     enddo
     deallocate(fieldmap)
!
     iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine nemsio_denseread8
!------------------------------------------------------------------------------
   subroutine readmpi4(gfile,fieldmapsize,fieldmap,data,iret,idispstt)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: read real 4 data out using MPI
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: fieldmapsize
     integer,intent(in)                  :: fieldmap(fieldmapsize)
     real(nemsio_realkind),intent(out)   :: data(fieldmapsize)
     integer,optional,intent(out)        :: iret
     integer(8),optional,intent(in)        :: idispstt
!local vars
     integer(MPI_OFFSET_KIND) :: idisp
     integer     :: filetype,ios
     integer                  :: status(MPI_STATUS_SIZE)
     real(nemsio_dblekind),allocatable   :: tmp(:)
!
!--- set file type
     if(trim(gfile%gdatatype).eq."bin4" ) then
!       print *,'be type create,fieldmapsize=',fieldmapsize,'fieldmap=',maxval(fieldmap),&
!        minval(fieldmap)
       call mpi_type_create_indexed_block(fieldmapsize,1,fieldmap, &
            MPI_REAL,filetype,ios)
!       print *,'af type_create,ios=',ios,'size=',fieldmapsize, &
!          'fieldmap=',maxval(fieldmap),minval(fieldmap)
     else if ( trim(gfile%gdatatype).eq."bin8" ) then
       call mpi_type_create_indexed_block(fieldmapsize,1,fieldmap, &
            MPI_REAL8,filetype,ios)
     endif
     call MPI_TYPE_COMMIT(filetype,iret)
!       print *,'af type_commit,ios=',iret
     if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop('stop at MPI set field map!')
       endif
     endif
!
!--- file set view, and read
     if(trim(gfile%gdatatype).eq."bin4") then

       if(present(idispstt)) then
         idisp=gfile%tlmeta+4+idispstt
       else
         idisp=gfile%tlmeta+4
       endif
       call mpi_file_set_view(gfile%fh,idisp,MPI_REAL4,filetype,'native', &
         MPI_INFO_NULL,ios)
!       print *,'af fiel_setview,ios=',ios
       call MPI_FILE_READ_ALL(gfile%fh,data,fieldmapsize,MPI_REAL4,  &
        status,ios)
!       print *,'af fiel_readall,ios=',ios
       if ( ios.ne.0 ) then
        if ( present(iret))  then
          iret=ios
          return
        else
          call nemsio_stop('stop at MPI read file all for bin4!')
        endif
       endif

      elseif (trim(gfile%gdatatype).eq."bin8") then

       allocate(tmp(fieldmapsize))
       if(present(idispstt)) then
         idisp=gfile%tlmeta+4+idispstt
       else
         idisp=gfile%tlmeta+4
       endif
       call mpi_file_set_view(gfile%fh,idisp,MPI_REAL8,filetype,'native', &
         MPI_INFO_NULL,ios)
       call MPI_FILE_READ_ALL(gfile%fh,tmp,fieldmapsize,MPI_REAL8,  &
         status,ios)
       if ( ios.ne.0 ) then
        if ( present(iret))  then
          iret=ios
          return
        else
          call nemsio_stop('stop at MPI read file all for bin8!')
        endif
       endif
       data=tmp
       deallocate(tmp)
!
      endif
!
!       print *,'end of readmpi4'
      iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine readmpi4
!------------------------------------------------------------------------------
   subroutine readmpi8(gfile,fieldmapsize,fieldmap,data,iret,idispstt)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: free unit number array index corresponding to unit number
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: fieldmapsize
     integer,intent(in)                  :: fieldmap(fieldmapsize)
     real(nemsio_dblekind),intent(out)   :: data(fieldmapsize)
     integer,optional,intent(out)        :: iret
     integer(8),optional,intent(in)        :: idispstt
!local vars
     integer(MPI_OFFSET_KIND) :: idisp
     integer     :: filetype,ios
     integer                  :: status(MPI_STATUS_SIZE)
     real(nemsio_realkind),allocatable   :: tmp(:)
!---
     iret=-25
!
!--- set file type
     if(trim(gfile%gdatatype).eq."bin4" ) then
       call mpi_type_create_indexed_block(fieldmapsize,1,fieldmap, &
            MPI_REAL,filetype,ios)
     else if ( trim(gfile%gdatatype).eq."bin8" ) then
       call mpi_type_create_indexed_block(fieldmapsize,1,fieldmap, &
            MPI_REAL8,filetype,ios)
     endif
     call MPI_TYPE_COMMIT(filetype,iret)
     if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop('stop at MPI set field map!')
       endif
     endif
!
!--- file set view, and read
     if(trim(gfile%gdatatype).eq."bin4") then

       allocate(tmp(fieldmapsize))
       if(present(idispstt)) then
         idisp=gfile%tlmeta+4+idispstt
       else
         idisp=gfile%tlmeta+4
       endif
       call mpi_file_set_view(gfile%fh,idisp,MPI_REAL4,filetype,'native', &
         MPI_INFO_NULL,ios)
       call MPI_FILE_READ_ALL(gfile%fh,tmp,fieldmapsize,MPI_REAL4,  &
        status,ios)
       if ( ios.ne.0 ) then
        if ( present(iret))  then
          iret=ios
          return
        else
          call nemsio_stop('stop at MPI read file all for bin4!')
        endif
       endif
       data(1:fieldmapsize)=tmp(1:fieldmapsize)

      elseif (trim(gfile%gdatatype).eq."bin8") then

       if(present(idispstt)) then
         idisp=gfile%tlmeta+4+idispstt
       else
         idisp=gfile%tlmeta+4
       endif
       call mpi_file_set_view(gfile%fh,idisp,MPI_REAL8,filetype,'native', &
         MPI_INFO_NULL,ios)
       call MPI_FILE_READ_ALL(gfile%fh,data,fieldmapsize,MPI_REAL8,  &
         status,ios)
       if ( ios.ne.0 ) then
        if ( present(iret))  then
          iret=ios
          return
        else
          call nemsio_stop('stop at MPI read file all for bin8!')
        endif
       endif

      endif
!
      iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine readmpi8
!------------------------------------------------------------------------------
    subroutine set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,iret,jrec)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: free unit number array index corresponding to unit number
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(in)     :: gfile
     integer,intent(in)                :: ista,iend,jsta,jend
     integer,intent(out)            :: fieldmap(:)
     integer,intent(out)               :: iret
     integer,optional,intent(in)       :: jrec
!-- local vars
     integer i,j,k,m,jm,km,nfieldsize,nfld,krec,kstart
!---
     iret=-20
!---
     if(trim(gfile%gdatatype).eq.'bin4') then
        nfieldsize=gfile%fieldsize+2
     elseif (trim(gfile%gdatatype).eq.'bin8') then
        nfieldsize=gfile%fieldsize+1
     endif
!---
     if(present(jrec)) then
       krec=jrec
       nfld=1
     else
       krec=1
       nfld=size(fieldmap)/((iend-ista+1)*(jend-jsta+1))
     endif
!--- set file map
     kstart=(krec-1)*nfieldsize
!     write(0,*)'in set_mpimap, kstart=',kstart,' tlmeta=',gfile%tlmeta,  &
!      ' nfieldsize=',nfieldsize,'krec=',krec,'nfld=',nfld,'fldsize=',gfile%fieldsize, &
!      'dimx=',gfile%dimx,'dimy=',gfile%dimy,'nfrmae=',gfile%nframe
!
     if (gfile%nframe.eq.0) then
       m=0
       do k=1,nfld
           km=(k-1)*nfieldsize+kstart-1
           do j=jsta,jend
             jm=(j-1)*gfile%dimx
             do i=ista,iend
               m=m+1
               fieldmap(m)=i+jm+km
             enddo
           enddo
        enddo
     else if(gfile%nframe.gt.0) then
       m=0
       do k=1,nfld
         km=(k-1)*nfieldsize+kstart-1
         do j=jsta,jend
           jm=(j-1)*(gfile%dimx+2*gfile%nframe)
           do i=ista,iend
             m=m+1
             fieldmap(m)=i+jm+km
           enddo
         enddo
       enddo
     endif
!     if (trim(gfile%gdatatype).eq.'bin8') gfile%fieldmap=gfile%fieldmap-1
!
!      print *,'Check field map size,',size(fieldmap), m,'end of set_mpimap'
!
     iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine set_mpimap_read
!------------------------------------------------------------------------------
    subroutine set_mpimap_wrt(gfile,ista,iend,jsta,jend,fieldmap,iret,jrec)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: free unit number array index corresponding to unit number
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(in)     :: gfile
     integer,intent(in)                :: ista,iend,jsta,jend
     integer,intent(out)               :: fieldmap(:)
     integer,intent(out)               :: iret
     integer,optional,intent(in)       :: jrec
!-- local vars
     integer i,j,k,m,jm,km,nfieldsize,nfld,krec,kstart,inum
!---
     iret=-20
!---
     if(present(jrec)) then
       krec=jrec
       nfld=1
     else
       krec=1
       nfld=size(fieldmap)/((iend-ista+1)*(jend-jsta+1))
     endif
!--- set file map
     nfieldsize=gfile%fieldsize+2
     kstart=(krec-1)*nfieldsize
!     print *,'in set_mpimap, kstart=',kstart,' tlmeta=',gfile%tlmeta,  &
!       ' nfieldsize=',nfieldsize,'krec=',krec,'nfld=',nfld,'fldsize=',gfile%fieldsize
!
!     if(gfile%mype.eq.gfile%lead_task ) then
!      if(size(fieldmap) < ((iend-ista+1)*(jend-jsta+1)+2)*nfld) &
!        call nemsio_stop('in set_mpimap, fieldmap size is too small!')
!     endif
!

     if (gfile%nframe.eq.0) then
       inum=gfile%dimx
     elseif(gfile%nframe.gt.0) then
       inum=gfile%dimx+2*gfile%nframe
     endif
!
     m=0
     do k=1,nfld
         km=(k-1)*nfieldsize+kstart
         if(gfile%mype.eq.gfile%lead_task) then
           m=m+1
           fieldmap(m)=km
         endif
         do j=jsta,jend
             jm=(j-1)*inum
             do i=ista,iend
               m=m+1
               fieldmap(m)=i+jm+km
             enddo
         enddo
         if(gfile%mype.eq.gfile%lead_task) then
           m=m+1
           fieldmap(m)=km+nfieldsize-1
         endif
     enddo
!     if (trim(gfile%gdatatype).eq.'bin8') gfile%fieldmap=gfile%fieldmap-1
!     
!      print *,'Check field map size,',size(fieldmap), m,'end of set_mpimap'
!
     iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine set_mpimap_wrt
!------------------------------------------------------------------------------
   subroutine nemsio_densewrite4(gfile,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: free unit number array index corresponding to unit number
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: ista,iend,jsta,jend
     real(nemsio_realkind),intent(out)   :: data(:)
     integer,optional,intent(out)        :: iret
!
     real(nemsio_dblekind),allocatable   :: data8(:)
!
     if(trim(gfile%gdatatype)=='bin4') then
      call  mpi_densewrite4(gfile,ista,iend,jsta,jend,data,iret)
     else if (trim(gfile%gdatatype)=='bin8') then
      allocate(data8(size(data)))
      data8=data
      call  mpi_densewrite8(gfile,ista,iend,jsta,jend,data8,iret)
      deallocate(data8)
     endif
!
   end subroutine nemsio_densewrite4
!------------------------------------------------------------------------------
   subroutine nemsio_densewrite8(gfile,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: free unit number array index corresponding to unit number
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: ista,iend,jsta,jend
     real(nemsio_dblekind),intent(out)   :: data(:)
     integer,optional,intent(out)        :: iret
!
     real(nemsio_realkind),allocatable   :: data4(:)
!
     if(trim(gfile%gdatatype)=='bin4') then
      allocate(data4(size(data)))
      data4=data
      call  mpi_densewrite4(gfile,ista,iend,jsta,jend,data4,iret)
      deallocate(data4)
     else if (trim(gfile%gdatatype)=='bin8') then
      call  mpi_densewrite8(gfile,ista,iend,jsta,jend,data,iret)
     endif
!
   end subroutine nemsio_densewrite8
!
!------------------------------------------------------------------------------
   subroutine mpi_densewrite4(gfile,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: free unit number array index corresponding to unit number
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: ista,iend,jsta,jend
     real(nemsio_realkind),intent(out)   :: data(:)
     integer,optional,intent(out)        :: iret
!--- local vars
     integer                  :: i,ierr,nfldsize,nfld,nfldloop,mfldrmd,k
     integer               :: fieldmapsize,fldmapsize1,fldmapsize,fielddatasize
     integer,allocatable      :: fieldmap(:)
     real(nemsio_realkind),allocatable      :: datatmp(:)
     integer ios,nstt,nend
     integer(8) idispstt
!---
     iret=-25
!
!--- check lead_task and last_task
     if(gfile%mype==gfile%lead_task .and. ista/=1 .and. jsta/=1 ) &
       call nemsio_stop("lead_task\'s subdomain must cover the (1,1) corner")
!
!--- set nfld
     nfldsize=gfile%fieldsize+2
     nfld=min(gfile%nrec,nemsio_maxint/nfldsize)
     nfldloop=(gfile%nrec-1)/nfld+1
     mfldrmd=mod(gfile%nrec,nfld)
!     write(0,*)'in dense read,nfld=',nfld,'nfldloop=',nfldloop, &
!       'mfldrmd=',mfldrmd
!
!--- set file map
     if(gfile%mype==gfile%lead_task) then
      fieldmapsize=((iend-ista+1)*(jend-jsta+1)+2)*nfld
      fldmapsize=(iend-ista+1)*(jend-jsta+1)+2
      fldmapsize1=(iend-ista+1)*(jend-jsta+1)
      fielddatasize=((iend-ista+1)*(jend-jsta+1)+2)*gfile%nrec
     else
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)*nfld
      fielddatasize=(iend-ista+1)*(jend-jsta+1)*gfile%nrec
     endif
!    print *,'in dense write, size(data)=',size(data),'fieldmapsize=',fieldmapsize,gfile%fieldsize_real4
     allocate(datatmp(fielddatasize))
     allocate(fieldmap(fieldmapsize) )
     call set_mpimap_wrt(gfile,ista,iend,jsta,jend,fieldmap,ierr)
     if(ierr.ne.0) return
!
!--- prepare data
     if(gfile%mype.eq.gfile%lead_task) then
      do i=1,gfile%nrec
        datatmp((i-1)*fldmapsize+1)=gfile%fieldsize_real4
        datatmp(i*fldmapsize)=datatmp(1)
        datatmp((i-1)*fldmapsize+2:i*fldmapsize-1)=data((i-1)*fldmapsize1+1:i*fldmapsize1)
      enddo
     else
      datatmp(:)=data(:)
     endif
!
!---
     do k=1,nfldloop
!
       if(k<nfldloop.or.mfldrmd==0) then
         nstt=(k-1)*fieldmapsize+1
         nend=k*fieldmapsize
       elseif(mfldrmd/=0) then
         nstt=(k-1)*fieldmapsize+1
         nend=gfile%nrec*(iend-ista+1)*(jend-jsta+1)
         deallocate(fieldmap)
         fieldmapsize=(iend-ista+1)*(jend-jsta+1)*mfldrmd
         allocate(fieldmap(fieldmapsize) )
         call set_mpimap_read(gfile,ista,iend,jsta,jend,fieldmap,ios)
         if(ios.ne.0) return
       endif
       idispstt=int(k-1,8)*int(nfld,8)*int(nfldsize*4,8)

       call writempi4(gfile,fieldmapsize,fieldmap,datatmp(nstt:nend),iret=iret, &
         idispstt=idispstt)
       if (iret.eq.0) return
!
     enddo
!
     iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine mpi_densewrite4
!------------------------------------------------------------------------------
   subroutine mpi_densewrite8(gfile,ista,iend,jsta,jend,data,iret)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: free unit number array index corresponding to unit number
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: ista,iend,jsta,jend
     real(nemsio_dblekind),intent(out)   :: data(:)
     integer,optional,intent(out)        :: iret
!--- local vars
     integer                  :: i,ierr,nfldsize,nfld,nfldloop,mfldrmd,k
     integer               :: fieldmapsize,fldmapsize,fldmapsize1,fielddatasize
     integer,allocatable      :: fieldmap(:)
     real(nemsio_dblekind),allocatable   :: datatmp(:)
     integer ios,nstt,nend
     integer(8) idispstt
!---
     iret=-25
!
!--- check lead_task and last_task
    if(gfile%mype==gfile%lead_task .and. ista/=1 .and. jsta/=1 ) &
       call nemsio_stop("lead_task\'s subdomain must cover the (1,1) corner")
!
!--- set nfld
    nfldsize=gfile%fieldsize+2
    nfld=min(gfile%nrec,nemsio_maxint/nfldsize)
    nfldloop=(gfile%nrec-1)/nfld+1
    mfldrmd=mod(gfile%nrec,nfld)
!    write(0,*)'in dense read,nfld=',nfld,'nfldloop=',nfldloop, &
!       'mfldrmd=',mfldrmd
!
!--- set file map
    if(gfile%mype==gfile%lead_task) then
      fieldmapsize=((iend-ista+1)*(jend-jsta+1)+2)*nfld
      fldmapsize=(iend-ista+1)*(jend-jsta+1)+2
      fldmapsize1=(iend-ista+1)*(jend-jsta+1)
      fielddatasize=((iend-ista+1)*(jend-jsta+1)+2)*gfile%nrec
    else
      fieldmapsize=(iend-ista+1)*(jend-jsta+1)*nfld
      fielddatasize=(iend-ista+1)*(jend-jsta+1)*gfile%nrec
    endif
    allocate(datatmp(fielddatasize))
    allocate(fieldmap(fieldmapsize) )
    call set_mpimap_wrt(gfile,ista,iend,jsta,jend,fieldmap,ierr)
!    print *,'in writerecbin4, after set_mpimap,ierr=',ierr
    if(ierr.ne.0) return
!
!--- prepare data
    if(gfile%mype.eq.gfile%lead_task) then
      do i=1,gfile%nrec
        datatmp((i-1)*fldmapsize+1)=gfile%fieldsize_real8
        datatmp(i*fldmapsize)=datatmp(1)
        datatmp((i-1)*fldmapsize+2:i*fldmapsize-1)=data((i-1)*fldmapsize1+1:i*fldmapsize1)
      enddo
    else
      datatmp(:)=data(:)
    endif
!
!---
     do k=1,nfldloop
!
       if(k<nfldloop.or.mfldrmd==0) then
         nstt=(k-1)*fieldmapsize+1
         nend=k*fieldmapsize
       elseif(mfldrmd/=0) then
         nstt=(k-1)*fieldmapsize+1
         nend=gfile%nrec*(iend-ista+1)*(jend-jsta+1)
         deallocate(fieldmap)
         fieldmapsize=(iend-ista+1)*(jend-jsta+1)*mfldrmd
         allocate(fieldmap(fieldmapsize) )
         call set_mpimap_wrt(gfile,ista,iend,jsta,jend,fieldmap,ios)
         if(ios.ne.0) return
       endif
       idispstt=int(k-1,8)*int(nfld,8)*int(gfile%fieldsize*8+8,8)
!
       call writempi8(gfile,fieldmapsize,fieldmap,datatmp,iret=iret,idispstt=idispstt)
     enddo
!
     iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine mpi_densewrite8
!------------------------------------------------------------------------------
   subroutine writempi4(gfile,fieldmapsize,fieldmap,data,iret,idispstt) 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: write out real 4 data using MPI
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: fieldmapsize
     integer,intent(in)                  :: fieldmap(:)
     real(nemsio_realkind),intent(in)    :: data(:)
     integer,optional,intent(out)        :: iret
     integer(8),optional,intent(in)      :: idispstt
!--- local vars
     integer(MPI_OFFSET_KIND) :: idisp
     integer                  :: status(MPI_STATUS_SIZE)
     integer ios,filetype
     real(nemsio_dblekind),allocatable :: tmp(:)
!
!--- set file type
     if(trim(gfile%gdatatype).eq."bin4" ) then
       call mpi_type_create_indexed_block(fieldmapsize,1,fieldmap, &
            MPI_REAL,filetype,ios)
     else if ( trim(gfile%gdatatype).eq."bin8" ) then
       call mpi_type_create_indexed_block(fieldmapsize,1,fieldmap, &
            MPI_REAL8,filetype,ios)
     endif
     call MPI_TYPE_COMMIT(filetype,ios)
     if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop('stop: at write set type indexed block')
       endif
     endif
!
!--- file set view, and read
!
       if(present(idispstt))  then
         idisp=gfile%tlmeta+idispstt
       else
         idisp=gfile%tlmeta
       endif
       call mpi_file_set_view(gfile%fh,idisp,MPI_REAL4,filetype,'native', &
         MPI_INFO_NULL,ios)
       call MPI_FILE_WRITE_ALL(gfile%fh,data,fieldmapsize,MPI_REAL4,  &
        status,ios)
       if ( ios.ne.0 ) then
         if ( present(iret))  then
           iret=ios
           return
         else
           call nemsio_stop('stop: at MPI write all for bin4')
         endif
       endif
!
      iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   end subroutine writempi4
!------------------------------------------------------------------------------
  subroutine writempi8(gfile,fieldmapsize,fieldmap,data,iret,idispstt)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
! abstract: write out real 4 data using MPI
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - -
      implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     type(nemsio_gfile),intent(inout)    :: gfile
     integer,intent(in)                  :: fieldmapsize
     integer,intent(in)                  :: fieldmap(:)
     real(nemsio_dblekind),intent(in)    :: data(:)
     integer,optional,intent(out)        :: iret
     integer(8),optional,intent(in)         :: idispstt
!--- local vars
     integer(MPI_OFFSET_KIND) :: idisp
     integer                  :: status(MPI_STATUS_SIZE)
     integer ios,filetype
     real(nemsio_realkind),allocatable :: tmp(:)
!---
     iret=-25
!
!--- set file type
     if(trim(gfile%gdatatype).eq."bin4" ) then
       call mpi_type_create_indexed_block(fieldmapsize,1,fieldmap, &
            MPI_REAL,filetype,ios)
     else if ( trim(gfile%gdatatype).eq."bin8" ) then
       call mpi_type_create_indexed_block(fieldmapsize,1,fieldmap, &
            MPI_REAL8,filetype,ios)
     endif
     call MPI_TYPE_COMMIT(filetype,iret)
     if ( ios.ne.0 ) then
       if ( present(iret))  then
         iret=ios
         return
       else
         call nemsio_stop('stop: at write set type indexed block')
       endif
     endif
!
!--- file set view, and read
!
       if(present(idispstt))  then
         idisp=gfile%tlmeta+idispstt
       else
         idisp=gfile%tlmeta
       endif
       call mpi_file_set_view(gfile%fh,idisp,MPI_REAL8,filetype,'native', &
         MPI_INFO_NULL,ios)
       call MPI_FILE_WRITE_ALL(gfile%fh,data,fieldmapsize,MPI_REAL8,  &
         status,ios)
       if ( ios.ne.0 ) then
         if ( present(iret))  then
           iret=ios
           return
         else
           call nemsio_stop('stop: at MPI write all for bin8')
         endif
       endif
!
      iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    end subroutine writempi8
!
!------------------------------------------------------------------------------
!
     elemental function equal_str_nocase(str1,str2)
!
!-----------------------------------------------------------------------
!
! convert a word to lower case
!
      logical              :: equal_str_nocase
      Character (len=*) , intent(in) :: str1
      Character (len=*) , intent(in) :: str2
      integer :: i,ic1,ic2,nlen
      nlen = len(str2)
!
      if(len(str1)/=nlen)  then
        equal_str_nocase=.false.
        return
      endif
      equal_str_nocase=.false.
      do i=1,nlen
        ic1 = ichar(str1(i:i))
        if (ic1 >= 65 .and. ic1 < 91) ic1 = ic1+32
        ic2 = ichar(str2(i:i))
        if (ic2 >= 65 .and. ic2 < 91) ic2 = ic2+32
        if(ic1/=ic2) then
           equal_str_nocase=.false.
           return
        endif
      end do
      equal_str_nocase=.true.
!
!-----------------------------------------------------------------------
!
      end function equal_str_nocase

!------------------------------------------------------------------------------
  end module module_nemsio_mpi
