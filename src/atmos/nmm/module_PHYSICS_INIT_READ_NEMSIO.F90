!----------------------------------------------------------------------
!
      MODULE MODULE_PHYSICS_INIT_READ_NEMSIO
!
!----------------------------------------------------------------------
!
      USE ESMF_MOD
      USE NEMSIO_MODULE_MPI
      USE MODULE_INCLUDE
      USE MODULE_PHYSICS_INTERNAL_STATE,ONLY: PHYSICS_INTERNAL_STATE
      USE MODULE_DM_PARALLEL,ONLY : IDS,IDE,JDS,JDE                     &
                                   ,IMS,IME,JMS,JME                     &
                                   ,ITS,ITE,JTS,JTE                     &
                                   ,DSTRB,IDSTRB                        &
                                   ,MYPE_SHARE,NUM_TILES
!
      USE MODULE_MICROPHYSICS_NMM
      USE MODULE_CONSTANTS,ONLY : G
      USE MODULE_EXCHANGE
      USE MODULE_CONTROL,ONLY : TIMEF
!
!----------------------------------------------------------------------
!
      IMPLICIT NONE
!
!----------------------------------------------------------------------
!
      PRIVATE
!
      PUBLIC :: PHYSICS_READ_INPUT_NEMSIO                              &
               ,PHYSICS_READ_RESTT_NEMSIO
!
!----------------------------------------------------------------------
!
      CONTAINS
!
!----------------------------------------------------------------------
!######################################################################
!----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_READ_INPUT_NEMSIO(INFILE                       &
                                          ,MYPE,MPI_COMM_COMP           &
                                          ,IDAT,IHRST,PT                &
                                          ,INT_STATE                    &
                                          ,NSOIL,LM                     &
                                          ,IDS,IDE,JDS,JDE              &
                                          ,IMS,IME,JMS,JME              &
                                          ,IRTN )
!
!-----------------------------------------------------------------------
!
        implicit none
!
!------------------------
!***  Argument variables
!------------------------
!
      CHARACTER(LEN=*),INTENT(IN) :: INFILE
      INTEGER,INTENT(IN) :: MYPE,MPI_COMM_COMP
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME,NSOIL
      INTEGER,DIMENSION(3),INTENT(OUT) :: IDAT
!
      INTEGER,INTENT(OUT) :: IHRST
      REAL,INTENT(OUT) :: PT
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE                 
!
      INTEGER,INTENT(OUT) :: IRTN
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: LPT2
      INTEGER :: LDIM1,LDIM2,UDIM1,UDIM2
      INTEGER :: N,I,J,L,K,II,JJ
!
      REAL :: PDTOP
      REAL,DIMENSION(LM+1) :: PSG1,SG1,SG2,SGM
      REAL,DIMENSION(LM) :: DSG1,DSG2,PDSG1,PSGML1,SGML1,SGML2
!
      REAL,DIMENSION(:),ALLOCATABLE :: tmp
!
      LOGICAL :: RUN
!
      TYPE(NEMSIO_GFILE) :: GFILE
!
      integer :: nrec,recn,fldsize,fldst,js,rc
      character(16),allocatable :: recname(:),reclevtyp(:)
      integer,allocatable   :: reclev(:)
!
      real(8) :: stime,etime
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      stime=timef()
!
!-----------------------------------------------------------------------
!***  First we need the value of PT (pressure at top of domain)
!-----------------------------------------------------------------------
!
      CALL NEMSIO_INIT()
!
      CALL NEMSIO_OPEN(gfile,INFILE,'read',mpi_comm_comp,iret=irtn)
      if(irtn/=0)then
        write(0,*)'ERROR: Unable to open file ',trim(infile)            &
                 ,' in PHYSICS_READ_INPUT_NEMSIO'
        write(0,*)' ABORTING!'
        call esmf_finalize(terminationflag=esmf_abort                   &
                          ,rc             =rc)
      endif
!
!-----------------------------------------------------------------------
!***  Vertical layer information is needed in order to send it to
!***  some specific schemes' initialization routines.
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'PT',int_state%PT,iret=irtn)
      PT=int_state%PT
      CALL NEMSIO_GETHEADVAR(gfile,'RUN',run,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'IHRST',ihrst,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'IDAT',idat,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'PDTOP',pdtop,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'LPT2',lpt2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGM',sgm,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SG1',sg1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DSG1',dsg1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGML1',sgml1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SG2',sg2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DSG2',dsg2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGML2',sgml2,iret=irtn)
      if(mype==0)write(0,*)'in phys,pt=',pt,'run=',run,'ihrst=',ihrst,'idat=',idat,   &
         'pdtop=',pdtop,'lpt2=',lpt2
      if(mype==0)write(0,*)'in phys,sgm=',sgm(1:10),'sg1=',sg1(1:10),maxval(sg1),minval(sg1),'sg2=',sg2(1:10), &
         maxval(sg2),minval(sg2),'dsg1=',dsg1(1:10),'dsg2=',dsg2(1:10),'sgml1=',sgml1(1:10),  &
         'sgml2=',sgml2(1:10)
!
!-----------------------------------------------------------------------
!***  Before moving on, transfer values to the internal state.
!-----------------------------------------------------------------------
!
!
      DO L=1,LM+1
        PSG1(L)=SG1(L)*PDTOP+PT
      ENDDO
      DO L=1,LM
        PDSG1(L)=DSG1(L)*PDTOP
        PSGML1(L)=SGML1(L)*PDTOP+PT
      ENDDO
!
      int_state%PDTOP=PDTOP
!
      DO L=1,LM
        int_state%DSG2(L)=DSG2(L)
        int_state%SGML2(L)=SGML2(L)
        int_state%PDSG1(L)=PDSG1(L)
        int_state%PSGML1(L)=PSGML1(L)
      ENDDO
!
      DO L=1,LM+1
        int_state%SG1(L)=SG1(L)
        int_state%PSG1(L)=PSG1(L)
        int_state%SG2(L)=SG2(L)
        int_state%SGM(L)=SGM(L)
      ENDDO
!
      call nemsio_getfilehead(gfile,nrec=nrec,iret=irtn)
      ALLOCATE(recname(nrec),reclevtyp(nrec),reclev(nrec))
!     
      call nemsio_getfilehead(gfile,recname=recname,reclevtyp=reclevtyp,   &
                              reclev=reclev,iret=irtn)
!
!-----------------------------------------------------------------------
!***  Proceed with getting fields from input file.
!***  NOTE: Two records were already read at the top of this routine.
!-----------------------------------------------------------------------
!
      fldsize=(jte-jts+1)*(ite-its+1)
      ALLOCATE(TMP(fldsize*nrec),STAT=I)
!
      call nemsio_denseread(gfile,its,ite,jts,jte,tmp,iret=irtn)
      if(irtn/=0) then
        print *,'WRONG: Could not read all the fields in the file!'
      endif
!
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)
!
!-----------------------------------------------------------------------
!***  FIS (Sfc Geopotential)
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'fis','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%FIS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%FIS,1,3,3)
!      write(0,*)'phycs fis=',maxval(int_state%fis),minval(int_state%fis),  &
!        'fis1=',int_state%fis(its:its+2,jte-2:jte)
!
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sm','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'phycs sm=',maxval(int_state%sm),minval(int_state%sm),  &
!        'sm=',int_state%sm(its:its+2,jte-2:jte)
!
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
      ENDDO
      ENDDO
!
      call getrecn(recname,reclevtyp,reclev,nrec,'dpres','hybrid sig lev',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%PD(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%PD,1,2,2)
!      write(0,*)'phycs pd=',maxval(int_state%pd),minval(int_state%pd),  &
!        'pd=',int_state%pd(its:its+2,jte-2:jte)
!

!-----------------------------------------------------------------------
!***  PINT (test add)
!-----------------------------------------------------------------------

!      DO L=1,LM+1
!      DO J=JMS,JME
!      DO I=IMS,IME
!        int_state%PINT(I,J,L)=0.
!      ENDDO
!      ENDDO
!      ENDDO

!      DO K=1,LM
!
!          do j=jts,jte
!            do i=its,ite
!          if (K .eq. 1) then
!            int_state%PINT(i,j,k)=pt
!          endif
!            int_state%PINT(i,j,k+1)=int_state%PINT(i,j,k)+dsg2(K)*int_state%PD(i,j)+pdsg1(K)
!            enddo
!          enddo
!
!      ENDDO

!      CALL HALO_EXCH(int_state%PINT,LM+1,2,2)

!       if (MYPE .eq. 0) then
!     DO L=1,LM+1
!       write(0,*) 'L, pint(5,5,L): ', L, int_state%PINT(5,5,L)
!     ENDDO
!       endif

!

!-----------------------------------------------------------------------
!***  U, V, T, Q, CW
!-----------------------------------------------------------------------
!
      DO K=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,K)=0.
        ENDDO
        ENDDO
!
        call getrecn(recname,reclevtyp,reclev,nrec,'ugrd','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%U(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      ENDDO
      CALL HALO_EXCH(int_state%U,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,K)=0.
        ENDDO
        ENDDO
!
        call getrecn(recname,reclevtyp,reclev,nrec,'vgrd','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%V(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      ENDDO
      CALL HALO_EXCH(int_state%V,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,K)=0.
        ENDDO
        ENDDO
!
        call getrecn(recname,reclevtyp,reclev,nrec,'tmp','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%T(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      ENDDO
      CALL HALO_EXCH(int_state%T,LM,2,2)
!      write(0,*)'phycs t=',maxval(int_state%t),minval(int_state%t), &
!       int_state%t(its:its+2,jts:jts+2,1)
!-----------------------------------------------------------------------
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)
!
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'albedo','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALBEDO(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'albase','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALBASE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
! **** EPSR
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'epsr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%EPSR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'phycs epsr=',maxval(int_state%epsr),minval(int_state%epsr), &
!       int_state%epsr(its:its+2,jts:jts+2)
!
!-----------------------------------------------------------------------
!*** SNOW ALBEDO
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'mxsnal','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%MXSNAL(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  SST/TSK
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'tskin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TSKIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'tskin=',maxval(int_state%TSKIN),minval(int_state%TSKIN), &
!       int_state%TSKIN(its:its+2,jts:jts+2)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'tsea','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SST(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'phycs sst=',maxval(int_state%sst),minval(int_state%sst), &
!         int_state%sst(its:its+1,jts:jts+1)
!
!-----------------------------------------------------------------------
!***  SNO, SICE, STC, SMC, ISLTYP, IVGTYP, VEGFRC
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sno','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SNO(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
      call getrecn(recname,reclevtyp,reclev,nrec,'si','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SI(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'phycs si=',maxval(int_state%si),minval(int_state%si), &
!         int_state%si(its:its+1,jts:jts+1)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sice','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SICE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'aft sice, sice=',maxval(int_state%SICE),minval(int_state%SICE)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'tg','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TG(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!
      call getrecn(recname,reclevtyp,reclev,nrec,'cmc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CMC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'phycs cmc=',maxval(int_state%cmc),minval(int_state%cmc), &
!         int_state%cmc(its:its+1,jts:jts+1)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'phycs sr=',maxval(int_state%sr),minval(int_state%sr), &
!         int_state%sr(its:its+1,jts:jts+1)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'ustar','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ustar(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'phycs ustar=',maxval(int_state%ustar),minval(int_state%ustar), &
!         int_state%ustar(its:its+1,jts:jts+1)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'zorl','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%Z0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!      write(0,*)'phycs zo=',maxval(int_state%z0),minval(int_state%z0), &
!         int_state%z0(its:its+1,jts:jts+1)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'z0base','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%Z0BASE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%Z0BASE,1,3,3)
!      write(0,*)'phycs zobase=',maxval(int_state%z0base),minval(int_state%z0base), &
!         int_state%z0base(its:its+1,jts:jts+1)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'stdh','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%STDH(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%STDH,1,3,3)
!      write(0,*)'phycs stdh=',maxval(int_state%stdh),minval(int_state%stdh), &
!         int_state%stdh(its:its+1,jts:jts+1)
!
      DO L=1,NSOIL
        call getrecn(recname,reclevtyp,reclev,nrec,'stc','soil layer',l,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%STC(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
        call getrecn(recname,reclevtyp,reclev,nrec,'smc','soil layer',l,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%SMC(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
        call getrecn(recname,reclevtyp,reclev,nrec,'sh2o','soil layer',l,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%SH2O(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!      write(0,*)'phycs stc=',maxval(int_state%stc(:,:,l)),minval(int_state%stc(:,:,l)), &
!         int_state%stc(its:its+1,jts:jts+1,l)
!
      ENDDO
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sltyp','sfc',1,recn)
      int_state%ISLTYP=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ISLTYP(i,j)=INT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!      write(0,*)'phycs isltyp=',maxval(int_state%isltyp),minval(int_state%isltyp), &
!         int_state%isltyp(its:its+1,jts:jts+1)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'vgtyp','sfc',1,recn)
      int_state%IVGTYP=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%IVGTYP(i,j)=INT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!      write(0,*)'phycs ivgtyp=',maxval(int_state%ivgtyp),minval(int_state%ivgtyp), &
!         int_state%ivgtyp(its:its+1,jts:jts+1)
!
      call getrecn(recname,reclevtyp,reclev,nrec,'vegfrc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%VEGFRC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!      write(0,*)'phycs vegfrv=',maxval(int_state%VEGFRC),minval(int_state%VEGFRC), &
!         int_state%VEGFRC(its:its+1,jts:jts+1)
!
!----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'SBD',int_state%SBD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'WBD',int_state%WBD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DPHD',int_state%DPHD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DLMD',int_state%DLMD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TPH0D',int_state%TPH0D,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TLM0D',int_state%TLM0D,iret=irtn)
!     write(0,*)'phycs sbd=',int_state%SBD,int_state%WBD,int_state%DPHD,  &
!      int_state%DLMD,int_state%TPH0D,int_state%TLM0D
!
!----------------------------------------------------------------------
!
      CALL NEMSIO_CLOSE(gfile)
!
      CALL NEMSIO_FINALIZE()
!
!----------------------------------------------------------------------
!
      DEALLOCATE(TMP)
!
       etime=timef()
      if(mype==0)then
        write(0,*)'PHYSICS_READ_INPUT_NEMSIO,time=',etime-stime
      endif
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHYSICS_READ_INPUT_NEMSIO
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_READ_RESTT_NEMSIO(INFILE                       &
                                          ,MYPE,MPI_COMM_COMP           &
                                          ,IYEAR_FCST,IMONTH_FCST       &
                                          ,IDAY_FCST,IHOUR_FCST         &
                                          ,IMINUTE_FCST,SECOND_FCST     &
                                          ,IHRST,IDAT,PT                &
                                          ,INT_STATE                    &
                                          ,NSOIL,LM                     &
                                          ,IDS,IDE,JDS,JDE              &
                                          ,IMS,IME,JMS,JME              &
                                          ,IRTN )
!
!-----------------------------------------------------------------------
!
        implicit none
!
!------------------------
!***  Argument variables
!------------------------
!
      CHARACTER(LEN=*),INTENT(IN) :: INFILE
!
      INTEGER,INTENT(IN) :: MYPE,MPI_COMM_COMP
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE,LM                          &
                           ,IMS,IME,JMS,JME
      INTEGER,INTENT(OUT) :: IYEAR_FCST,IMONTH_FCST,IDAY_FCST           &
                           ,IHOUR_FCST,IMINUTE_FCST,IHRST
      INTEGER,INTENT(OUT) :: NSOIL
!
      REAL,INTENT(OUT) :: SECOND_FCST
!
      INTEGER,DIMENSION(3),INTENT(OUT) :: IDAT
!
      REAL,INTENT(OUT) :: PT
!
      TYPE(PHYSICS_INTERNAL_STATE),POINTER :: INT_STATE
!
      INTEGER,INTENT(OUT) :: IRTN
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: N,I,J,K,L
      INTEGER :: LDIM1,LDIM2,UDIM1,UDIM2
      INTEGER :: LPT2
      INTEGER,DIMENSION(7) :: FCSTDATE
!
      REAL :: PDTOP
      REAL,DIMENSION(LM) :: DSG1,DSG2,SGML1,SGML2
      REAL,DIMENSION(LM+1) :: SG1,SG2,SGM
      REAL,ALLOCATABLE,DIMENSION(:) :: SLDPTH
      REAL,DIMENSION(LM) :: PDSG1,PSGML1
      REAL,DIMENSION(LM+1) :: PSG1
      REAL,DIMENSION(:),ALLOCATABLE :: tmp
!
      LOGICAL :: RUN
!
      CHARACTER(10) :: VARNAME
!
      TYPE(NEMSIO_GFILE) :: GFILE
!
      integer :: fldsize,fldst,js,nrec,rc,recn
      character(16),allocatable :: recname(:),reclevtyp(:)
      integer,allocatable   :: reclev(:)
      
      real(8) :: stime,etime,stime1,stime2
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
        stime=timef()

       CALL NEMSIO_INIT()
!
       CALL NEMSIO_OPEN(gfile,INFILE,'read',mpi_comm_comp,iret=irtn)
        if(irtn/=0)then
          write(0,*)'ERROR: Unable to open file ',trim(infile)          &
                   ,' in PHYSICS_READ_RESTT_NEMSIO'
          write(0,*)' ABORTING!'
          call esmf_finalize(terminationflag=esmf_abort                 &
                            ,rc             =rc)
        endif
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'FCSTDATE',FCSTDATE,iret=irtn)
      IYEAR_FCST=FCSTDATE(1)
      IMONTH_FCST=FCSTDATE(2)
      IDAY_FCST=FCSTDATE(3)
      IHOUR_FCST=FCSTDATE(4)
      IMINUTE_FCST=FCSTDATE(5)
      SECOND_FCST=0.
      if(FCSTDATE(7)/=0) SECOND_FCST=FCSTDATE(6)/(FCSTDATE(7)*1.)
!
      CALL NEMSIO_GETHEADVAR(gfile,'IHRST',IHRST,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'LPT2',LPT2,iret=irtn)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 1D arrays
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'IDAT',IDAT,iret=irtn)
!
      IF(MYPE==0)THEN
        write(0,*)'**** read in physics ****************'
        write(0,*)' Restart year =',iyear_fcst
        write(0,*)' Restart month=',imonth_fcst
        write(0,*)' Restart day  =',iday_fcst
        write(0,*)' Restart hour =',ihour_fcst
        write(0,*)' Original start year =',idat(3)
        write(0,*)' Original start month=',idat(2)
        write(0,*)' Original start day  =',idat(1)
        write(0,*)' Original start hour =',ihrst
        write(0,*)' Timestep   =',int_state%dt
        write(0,*)' Steps/hour =',3600./int_state%dt
        write(0,*)'*************************************'
      ENDIF
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'NSOIL',NSOIL,iret=irtn)
!      write(0,*)'in rst, fcstdate=',fcstdate,'ihrst=',ihrst,'lpt2=',lpt2, &
!        'idat=',idat,'nsoil=',nsoil
      ALLOCATE(SLDPTH(1:NSOIL))
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real scalars
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'PDTOP',PDTOP,iret=irtn)
!
!-----------------------------------------------------------------------
      CALL NEMSIO_GETHEADVAR(gfile,'PT',PT,iret=irtn)
      int_state%PT=PT
!
      CALL NEMSIO_GETHEADVAR(gfile,'TLM0D',int_state%TLM0D,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TPH0D',int_state%TPH0D,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DPHD',int_state%DPHD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DLMD',int_state%DLMD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SBD',int_state%SBD,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'WBD',int_state%WBD,iret=irtn)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 1D arrays
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'SG1',SG1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SG2',SG2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DSG1',DSG1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'DSG2',DSG2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGML1',SGML1,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGML2',SGML2,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SGM',SGM,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'SLDPTH',SLDPTH,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'MP_RESTART',int_state%MP_RESTART_STATE,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TBPVS_STAT',int_state%TBPVS_STATE,iret=irtn)
      CALL NEMSIO_GETHEADVAR(gfile,'TBPVS0_STA',int_state%TBPVS0_STATE,iret=irtn)
!
      DO N=1,NSOIL
        int_state%SLDPTH(N)=SLDPTH(N)
      ENDDO
!
!-----------------------------------------------------------------------
!*** Read from restart file: Logical
!-----------------------------------------------------------------------
!
      CALL NEMSIO_GETHEADVAR(gfile,'RUN',RUN,iret=irtn)
!
!-----------------------------------------------------------------------
!
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!
!-----------------------------------------------------------------------
!
      DO L=1,LM
        PDSG1(L)=DSG1(L)*PDTOP
        PSGML1(L)=SGML1(L)*PDTOP+PT
      ENDDO
!
      DO L=1,LM+1
        PSG1(L)=SG1(L)*PDTOP+PT
      ENDDO

!
!-----------------------------------------------------------------------
!***  Transfer values to the internal state.
!-----------------------------------------------------------------------
!
      int_state%PDTOP=PDTOP
!
      DO L=1,LM
        int_state%DSG2(L)=DSG2(L)
        int_state%PDSG1(L)=PDSG1(L)
        int_state%PSGML1(L)=PSGML1(L)
        int_state%SGML2(L)=SGML2(L)
      ENDDO
!
      DO L=1,LM+1
        int_state%SG1(L)=SG1(L)
        int_state%PSG1(L)=PSG1(L)
        int_state%SG2(L)=SG2(L)
        int_state%SGM(L)=SGM(L)
      ENDDO
!
!-----------------------------------------------------------------------
!***  Read from restart file: Record info including name,levtyp, and lev
!-----------------------------------------------------------------------
!
      call nemsio_getfilehead(gfile,nrec=nrec,iret=irtn)
      ALLOCATE(recname(nrec),reclevtyp(nrec),reclev(nrec))
!
      call nemsio_getfilehead(gfile,recname=recname,reclevtyp=reclevtyp,   &
                              reclev=reclev,iret=irtn)
!
!-----------------------------------------------------------------------
!***  Proceed with getting all the fields from input file.
!-----------------------------------------------------------------------
!
      fldsize=(jte-jts+1)*(ite-its+1)
      ALLOCATE(TMP(fldsize*nrec),STAT=I)
!
      stime1=timef()
      call nemsio_denseread(gfile,its,ite,jts,jte,tmp,iret=irtn)
!     write(0,*)'after dense read, time=',timef()-stime1,'tmp=',maxval(tmp),minval(tmp)
!
!-----------------------------------------------------------------------
!***  close nemsio file
!-----------------------------------------------------------------------
!
      CALL NEMSIO_CLOSE(GFILE)
!-----------------------------------------------------------------------
!
      CALL NEMSIO_FINALIZE()
!
!-----------------------------------------
!***  I and J limits for tracer variables
!-----------------------------------------
!
      LDIM1=LBOUND(int_state%Q,1)
      UDIM1=UBOUND(int_state%Q,1)
      LDIM2=LBOUND(int_state%Q,2)
      UDIM2=UBOUND(int_state%Q,2)
!
!-----------------------------------------------------------------------
!***  assign data: Integer 2D arrays
!-----------------------------------------------------------------------
!
      call getrecn(recname,reclevtyp,reclev,nrec,'sltyp','sfc',1,recn)
      int_state%ISLTYP=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ISLTYP(i,j)=NINT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
      call getrecn(recname,reclevtyp,reclev,nrec,'vgtyp','sfc',1,recn)
      int_state%IVGTYP=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%IVGTYP(i,j)=NINT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
      call getrecn(recname,reclevtyp,reclev,nrec,'cfrcv','sfc',1,recn)
      int_state%NCFRCV=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%NCFRCV(i,j)=NINT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
      call getrecn(recname,reclevtyp,reclev,nrec,'cfrst','sfc',1,recn)
      int_state%NCFRST=0
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%NCFRST(i,j)=NINT(tmp(i-its+1+js+fldst))
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Assign data: Real 2D arrays
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  FIS (Sfc Geopotential)
!-----------------------------------------------------------------------
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%FIS(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'fis','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%FIS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%FIS,1,3,3)
!
!-----------------------------------------------------------------------
!***  PD
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%PD(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'dpres','hybrid sig lev',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%PD(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%PD,1,2,2)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays (contd.)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ACFRCV
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRCV(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acfrcv','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACFRCV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACFRST
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRST(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acfrst','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACFRST(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACPREC
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACPREC(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acprec','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACPREC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACSNOM
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOM(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acsnom','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACSNOM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACSNOW
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOW(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'acsnow','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACSNOW(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AKHS_OUT
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKHS_OUT(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'akhs_out','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AKHS_OUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------

!***  AKMS_OUT
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKMS_OUT(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'akms_out','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AKMS_OUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ALBASE(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'albase','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALBASE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'albedo','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALBEDO(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALWIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'alwin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALWIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALWOUT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'alwout','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALWOUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ALWTOA
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'alwtoa','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ALWTOA(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ASWIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'aswin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ASWIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ASWOUT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'aswout','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ASWOUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ASWTOA
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'aswtoa','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ASWTOA(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  BGROFF
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'bgroff','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%BGROFF(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CFRACH
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cfrach','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CFRACH(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CFRACL
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cfracl','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CFRACL(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CFRACM
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cfracm','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CFRACM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CLDEFI
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cldefi','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CLDEFI(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cmc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CMC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CNVBOT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cnvbot','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CNVBOT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CNVTOP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cnvtop','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CNVTOP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CPRATE
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cprate','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CPRATE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CUPPT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cuppt','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CUPPT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CUPREC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'cuprec','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CUPREC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CZEN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'czen','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CZEN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  CZMEAN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'czmean','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%CZMEAN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  EPSR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'epsr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%EPSR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  GRNFLX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'grnflx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%GRNFLX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HBOTD
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'hbotd','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HBOTD(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HBOTS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'hbots','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HBOTS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HTOPD
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'htopd','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HTOPD(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HTOPS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'htops','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HTOPS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SNOW ALBEDO
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'mxsnal','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%MXSNAL(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  PBLH
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'pblh','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%PBLH(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  POTEVP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'potevp','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%POTEVP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  PREC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'prec','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%PREC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  PSHLTR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'pshltr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%PSHLTR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  Q10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'q10','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%Q10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  QSH
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'qsh','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%QSH(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  QSHLTR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'qshltr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%QSHLTR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  QWBS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'qwbs','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%QWBS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  QZ0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'qz0','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%QZ0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RADOT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'radot','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RADOT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RLWIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rlwin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RLWIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RLWTOA
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rlwtoa','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RLWTOA(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RSWIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rswin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RSWIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RSWINC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rswinc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%rswinc(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RSWOUT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rswout','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RSWOUT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SFCEVP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sfcevp','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SFCEVP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SFCEXC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sfcexc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SFCEXC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SFCLHX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sfclhx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SFCLHX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SFCSHX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sfcshx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SFCSHX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'si','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SI(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sice','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SICE(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SIGT4
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sigt4','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SIGT4(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      call getrecn(recname,reclevtyp,reclev,nrec,'sm','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SMSTAV
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'smstav','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SMSTAV(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SMSTOT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'smstot','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SMSTOT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SNO
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sno','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SNO(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SNOPCX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'snopcx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SNOPCX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SOILTB
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'soiltb','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SOILTB(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'sr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SSROFF
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'ssroff','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SSROFF(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SST
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tsea','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SST(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  SUBSHX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'subshx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%SUBSHX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TG
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tg','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TG(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TH10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'th10','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TH10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  THS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'ths','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%THS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  THZ0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'thz0','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%THZ0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TSHLTR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tshltr','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TSHLTR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TWBS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'twbs','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TWBS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  U10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'u10','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%U10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'uustar','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%USTAR(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  UZ0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'uz0','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%UZ0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%UZ0,1,3,3)
!-----------------------------------------------------------------------
!***  V10
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'v10','10 m above gnd',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%V10(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'vegfrc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%VEGFRC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  VZ0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'vz0','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%VZ0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%VZ0,1,3,3)
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'zorl','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%Z0(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!-----------------------------------------------------------------------
!***  TSKIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tskin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TSKIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AKHS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'akhs','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AKHS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AKMS
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'akms','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AKMS(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HBOT
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'hbot','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HBOT(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  HTOP
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'htop','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%HTOP(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RSWTOA
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rswtoa','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RSWTOA(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  POTFLX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'potflx','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%POTFLX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  RMOL
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'rmol','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%RMOL(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  T2
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'t2','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%T2(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'z0base','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%z0base(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TLMIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tlmin','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TLMIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  TLMAX
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'tlmax','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%TLMAX(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ACUTIM
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'acutim','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ACUTIM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  APHTIM
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'aphtim','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%APHTIM(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ARDLW
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'ardlw','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ARDLW(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ARDSW
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'ardsw','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ARDSW(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  ASRFC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'asrfc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%ASRFC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AVRAIN
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'avrain','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AVRAIN(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!***  AVCNVC
!-----------------------------------------------------------------------
      call getrecn(recname,reclevtyp,reclev,nrec,'avcnvc','sfc',1,recn)
      if(recn>0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            int_state%AVCNVC(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!              ASSIGN DATA: REAL 3D ARRAYS
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  U, V, T, Q, Q2, CW, F_ICE, F_RIMEF, F_RAIN
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%CLDFRA(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'cldfra','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%CLDFRA(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%Q2(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'q2','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%Q2(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif

      ENDDO
      CALL HALO_EXCH(int_state%Q2,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RLWTT(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'rlwtt','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%RLWTT(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%RSWTT(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'rswtt','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%RSWTT(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%T(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'tmp','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%T(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
      CALL HALO_EXCH(int_state%T,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TCUCN(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'tcucn','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%TCUCN(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRAIN(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'train','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%TRAIN(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%U(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'ugrd','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%U(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
      CALL HALO_EXCH(int_state%U,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%V(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'vgrd','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%V(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
      CALL HALO_EXCH(int_state%V,LM,2,2)
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%XLEN_MIX(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'xlen_mix','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%XLEN_MIX(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_ICE(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'f_ice','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%F_ICE(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RIMEF(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'f_rimef','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%F_RIMEF(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RAIN(I,J,K)=0.
        ENDDO
        ENDDO
        call getrecn(recname,reclevtyp,reclev,nrec,'f_rain','mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%F_RAIN(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
!-----------------------------------------------------------------------
!***  SH2O, SMC, STC
!-----------------------------------------------------------------------
!
      DO K=1,NSOIL
!
        call getrecn(recname,reclevtyp,reclev,nrec,'sh2o','soil layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%SH2O(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif

        call getrecn(recname,reclevtyp,reclev,nrec,'smc','soil layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%SMC(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
        call getrecn(recname,reclevtyp,reclev,nrec,'stc','soil layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%STC(i,j,k)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif

      ENDDO
!
!-----------------------------------------------------------------------
!***  TRACERS
!-----------------------------------------------------------------------
      DO N=int_state%INDX_O3+1,int_state%NUM_TRACERS_TOTAL                !<-- The first 'INDX_O3' arrays are unallocated pointers
      DO K=1,LM
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%TRACERS(I,J,K,N)=0.
        ENDDO
        ENDDO
!
        write(varname,'(a8,I2.2)')'tracers_',N
        call getrecn(recname,reclevtyp,reclev,nrec,varname,'mid layer',k,recn)
        if(recn>0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              int_state%TRACERS(i,j,k,n)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
!
      ENDDO
      ENDDO
      CALL HALO_EXCH(int_state%TRACERS,LM,int_state%NUM_TRACERS_TOTAL,1,2,2)
!-----------------------------------------------------------------------
!
      DEALLOCATE(TMP)
      DEALLOCATE(SLDPTH)
!
!-----------------------------------------------------------------------
      CALL MPI_BARRIER(MPI_COMM_COMP,IRTN)
!-----------------------------------------------------------------------
      etime=timef()
!     write(0,*)'PHYSICS_READ_RESTT_NEMSIO,time=',etime-stime
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      END SUBROUTINE PHYSICS_READ_RESTT_NEMSIO
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      SUBROUTINE getrecn(recname,reclevtyp,reclev,nrec,fldname,          &
                         fldlevtyp,fldlev,recn)
!-----------------------------------------------------------------------
!-- this subroutine searches the field list to find out a specific field,
!-- and return the field number for that field
!-----------------------------------------------------------------------
!
        implicit none
!
        integer,intent(in)      :: nrec
        character(*),intent(in) :: recname(nrec)
        character(*),intent(in) :: reclevtyp(nrec)
        integer,intent(in)      :: reclev(nrec)
        character(*),intent(in) :: fldname
        character(*),intent(in) :: fldlevtyp
        integer,intent(in)      :: fldlev
        integer,intent(out)     :: recn
!
        integer i
!
        recn=0
        do i=1,nrec
          if(trim(recname(i))==trim(fldname).and.                         &
            trim(reclevtyp(i))==trim(fldlevtyp) .and.                    &
            reclev(i)==fldlev) then
            recn=i
            return
          endif
        enddo
!
        if(recn==0) print *,'WARNING: field ',trim(fldname),' ',         &
          trim(fldlevtyp),' ',fldlev,' is not in the nemsio file!'

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------

      END SUBROUTINE getrecn
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
      END MODULE MODULE_PHYSICS_INIT_READ_NEMSIO
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
