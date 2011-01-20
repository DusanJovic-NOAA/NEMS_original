       SUBROUTINE PARA_FIXIO_W(IOPROC,sfc_fld, cfile,xhour,idate,
     &         lats_nodes_r,global_lats_r,lonsperlar,
     &         phy_f3d,phy_f2d,ngptc,nblck,ens_nam)
!!
! !revision history:
!  Sep      2010   Jun Wang  change to nemsio file
!  Dec      2010   Jun Wang  change to nemsio library
!
      use resol_def,     ONLY: latr, lonr, levs, lsoil, ivssfc_restart,
     &                         num_p2d, num_p3d
      use layout1,       ONLY: lats_node_r, me,ipt_lats_node_r,nodes
      use nemsio_module
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      USE machine,   ONLY: kind_io4, kind_ior, kind_io8,kind_phys
      implicit none
!!
      TYPE(Sfc_Var_Data),intent(in)    :: sfc_fld
      integer,intent(in)               :: idate(4),ioproc
      real(kind=kind_io8),intent(in)   :: xhour
      character*(*),intent(in)         :: cfile,ens_nam
      INTEGER,intent(in)               :: lats_nodes_r(nodes)
      INTEGER,intent(in)               :: GLOBAL_LATS_R(latr)
      INTEGER,intent(in)               :: lonsperlar(latr)
!
      integer,intent(in) :: ngptc, nblck
      REAL (KIND=KIND_phys),intent(in) ::
     &            phy_f3d(ngptc,levs,nblck,LATS_NODE_R,num_p3d)
     &,           phy_f2d(LONR,LATS_NODE_R,num_p2d)
!
!!
      real(kind=kind_io8),allocatable:: bfo(:,:)
      integer k,lan,i,nphyfld,fieldsize
!!
      integer,save:: version
!
      CHARACTER*2 nump3d,nump2d
      type(nemsio_gfile)  :: gfile
!
      integer iret,nrec,nmetavari,nmetavarr8,nmetaaryi,nmetaaryr8,
     &    idate7(7),nmeta,jrec,l,nsfcrec,nrecs,nrecs1
      integer nfhour,nfminute,nfsecondn,nfsecondd,nsoil
      character(16),allocatable :: variname(:),varr8name(:)
      character(16),allocatable :: aryiname(:),aryr8name(:)
      integer,allocatable :: varival(:),aryilen(:),aryr8len(:),
     &     aryival(:,:)
      real(kind=kind_io8),allocatable :: varr8val(:),aryr8val(:,:)
      character(16),allocatable :: recname(:),reclevtyp(:)
      integer,allocatable :: reclev(:)
      logical first
      sAve first,nsoil,nmetavari,variname,varival,
     &     nmetavarr8,varr8name,varr8val,nmetaaryi,aryiname,aryilen,
     &     aryival,nmetaaryr8,aryr8name,aryr8len,aryr8val,
     &     nmeta,nrec,nrecs,nrecs1,recname,reclevtyp,reclev
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!      print *,'in para_fix_w'
      nsfcrec=31
      fieldsize=sum(lonsperlar)
      nphyfld=nsfcrec+3*lsoil+num_p2d+num_p3d*levs
      allocate(bfo(fieldsize,nphyfld))
!
      call fld_collect(sfc_fld,phy_f2d,phy_f3d,ngptc,nblck,
     &  fieldsize,nphyfld,bfo,lats_nodes_r,global_lats_r,lonsperlar,
     &  ioproc,nsfcrec)
!       write(0,*)'after fld_collect'
!
      if (me.eq.ioproc) then
        if (first) then
          nsoil   = lsoil
          nmeta=5

          nmetavari=2
          allocate(variname(nmetavari),varival(nmetavari))
          variname=(/'ivs   ','irealf'/)
          varival=(/ivssfc_restart, 2/)
!          print *,'after vari'

          nmetaaryi=1
          allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
          aryiname(1)='lpl'
          aryilen(1)=latr/2
          allocate(aryival(maxval(aryilen),nmetaaryi))
          aryival(1:latr/2,1)=lonsperlar(1:latr/2)
!          print *,'after aryi'
!
          nmetaaryr8=1
          allocate(aryr8name(nmetaaryr8),aryr8len(nmetaaryr8))
          aryr8name(1)='zsoil'
          aryr8len(1)=lsoil
          allocate(aryr8val(maxval(aryr8len),nmetaaryr8))
          if (lsoil .eq. 4) then
            aryr8val(1:lsoil,1) = (/-0.1,-0.4,-1.0,-2.0/)
          elseif (lsoil .eq. 2) then
            aryr8val(1:lsoil,1) = (/-0.1,-2.0/)
          endif
!          print *,'after aryr8'
!
          nmetavarr8=1
          allocate(varr8name(nmetavarr8),varr8val(nmetavarr8))
          varr8name=(/'fhour'/)
          print *,'after varr'
!
          nrecs=nsfcrec
          nrecs1=nrecs+1
          nrec=nrecs+3*nsoil+num_p3d*levs+num_p2d
!          write(0,*)'after nrec=',nrec,'nsoil=',nsoil,'num_p3d=',
!     &      num_p3d,'num_p2d=',num_p2d
          allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
!record name
          RECNAME(1)='tmp'
          RECNAME(2)='weasd'
          RECNAME(3)='tg3'
          RECNAME(4)='sfcr'
!          RECNAME(5)='tcdc'
!          RECNAME(6)='pres'
!          RECNAME(7)='pres'
          RECNAME(5)='alvsf'
          RECNAME(6)='alvwf'
          RECNAME(7)='alnsf'
          RECNAME(8)='alnwf'
          RECNAME(9)='land'
          RECNAME(10)='veg'
          RECNAME(11)='cnwat'
          RECNAME(12)='f10m'
!          RECNAME(13)='tmp'
!          RECNAME(14)='spfh'
          RECNAME(13)='vtype'
          RECNAME(14)='sotyp'
          RECNAME(15)='facsf'
          RECNAME(16)='facwf'
          RECNAME(17)='fricv'
          RECNAME(18)='ffhh'
          RECNAME(19)='ffmm'
          RECNAME(20)='icetk'
          RECNAME(21)='icec'
          RECNAME(22)='tisfc'
          RECNAME(23)='tprcp'
          RECNAME(24)='crain'
          RECNAME(25)='snod'
          RECNAME(26)='shdmin'
          RECNAME(27)='shdmax'
          RECNAME(28)='sltyp'
          RECNAME(29)='salbd'
          RECNAME(30)='orog'
          RECNAME(31)='sncovr'
!
          RECNAME(nrecs1:nrecs+NSOIL)='smc'
          RECNAME(NSOIL+nrecs1:nrecs+2*NSOIL)='stc'
          RECNAME(2*NSOIL+nrecs1:nrecs+3*NSOIL)='slc'
          DO k=1,num_p2d
            write(nump2d,'(I2.2)')k
            RECNAME(3*NSOIL+nrecs+k)='phyf2d_'//nump2d
          enddo
          DO k=1,num_p3d
            write(nump3d,'(I2.2)')k
            RECNAME(3*NSOIL+num_p2d+nrecs1+(k-1)*levs:
     &        nrecs+3*NSOIL+num_p2d+k*levs)='phyf3d_'//nump3d
          enddo
!          print *,'after recname=',recname(36:40),recname(46:50)
!
          RECLEVTYP(1:nrecs)='sfc'
          RECLEVTYP(12)='10 m above gnd'
          RECLEVTYP(nrecs1:nrecs+3*nsoil)='soil layer'
          RECLEVTYP(nrecs1+3*nsoil:nrecs+3*nsoil+num_p2d)='sfc'
          RECLEVTYP(nrecs1+3*nsoil+num_p2d:nrecs+
     &       3*nsoil+num_p3d*levs+num_p2d)='mid layer'
!          print *,'after reclevtyp=',reclevtyp(36:40),reclevtyp(46:50)
!
          RECLEV(1:nrecs)=1
          DO K=1,NSOIL
            RECLEV(K+nrecs)=K
            RECLEV(NSOIL+K+nrecs)=K
            RECLEV(2*NSOIL+K+nrecs)=K
          ENDDO
          RECLEV(nrecs1+3*nsoil:nrecs+3*nsoil+num_p2d)=1
          do l=1,num_p3d
          DO K=1,levs
            RECLEV(nrecs+3*nsoil+num_p2d+(l-1)*levs+k)=k
          enddo
          enddo
!          print *,'after reclev=',reclev(36:40),reclev(46:50)
!
!endif first 
        endif
!
        varr8val(1)=xhour
!
        idate7(1:6)=0;idate7(7)=1
        idate7(1)=idate(4)
        idate7(2:3)=idate(2:3)
        idate7(4)=idate(1)
        nfhour=int(xhour)
        nfminute=int((xhour-nfhour)*60.)
        nfsecondn=int((xhour-nfhour)*3600.-nfminute*60.)
        nfsecondd=1
!
        PRINT 99,xhour,IDATE
99      FORMAT(1H ,'in fixio HOUR=',f8.2,3x,'IDATE=',
     &  4(1X,I4))
!!
! open nemsio sfc restart file
        call nemsio_init()
!
        write(0,*)'before nemsio_open for restart file'
        call nemsio_open(gfile,trim(cfile),'write',iret,   
     & modelname='GFS',gdatatype='bin8',idate=idate7,nfhour=nfhour, 
     & nfminute=nfminute,nfsecondn=nfsecondn,nfsecondd=nfsecondd,   
     & dimx=fieldsize,dimy=1,dimz=levs,nsoil=nsoil,nrec=nrec,         
     & nmeta=5,recname=recname,reclevtyp=reclevtyp,reclev=reclev,   
     & extrameta=.true.,nmetavari=nmetavari,nmetaaryi=nmetaaryi,    
     & nmetaaryr8=nmetaaryr8,nmetavarr8=nmetavarr8,variname=variname,   
     & varival=varival,varr8name=varr8name,varr8val=varr8val,         
     & aryiname=aryiname,aryilen=aryilen,aryival=aryival, 
     & aryr8name=aryr8name,aryr8len=aryr8len,aryr8val=aryr8val) 

       print *,'after restart nemsio_open, iret=',iret

        do jrec=1,nrec
          call nemsio_writerec(gfile,jrec,bfo(:,jrec),iret=iret)
        enddo
        deallocate(bfo)
!
        call nemsio_close(gfile)
        call nemsio_finalize()
!
! endof ioproc
      ENDIF
!
      if(first) first=.false.
         
      return
      end
