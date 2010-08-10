
      subroutine wrtout_physics(phour,fhour,zhour,idate,
     &                  sl,si,
     &                  sfc_fld, flx_fld, g2d_fld, 
     &                  fluxr,
     &                  lats_nodes_r,global_lats_r,lonsperlar,nblck,
     &                  colat1,cfhour1,pl_coeff)
!!
!
! May 2009 Jun Wang, modified to use write grid component
! Jan 2010 Sarah Lu, AOD added to flx files
! Feb 2010 Jun Wang, write out restart file
! Jul 2010 Sarah Lu, write out aerosol diag files (for g2d_fld)
! Aug 2010 Sarah Lu, scale the 2d_aer_diag by 1.e6
!

      use resol_def,               ONLY: latr, levs, levp1, lonr, nfxr
      use layout1,                 ONLY: me, nodes, lats_node_r, 
     &                                   nodes_comp
      use namelist_physics_def,    ONLY: gen_coord_hybrid, ldiag3d, 
     &                                   hybrid, fhlwr, fhswr, ens_nam
      use mpi_def,                 ONLY: liope, info, mpi_comm_all, 
     &                                   mc_comp, mpi_comm_null,quilting
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data, Flx_Var_Data
      use gfs_physics_g2d_mod,     ONLY: G2D_Var_Data
      USE machine,                 ONLY: kind_evod, kind_io8
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      TYPE(G2D_Var_Data)        :: g2d_fld
      CHARACTER(16)             :: CFHOUR1    ! for ESMF Export State Creation
      integer ixgr, pl_coeff
      real(kind=kind_evod) phour,fhour,zhour
!     real(kind=kind_evod) phour,fhour,zhour, xgf
!!
      integer              idate(4),nblck,km,iostat,no3d,ks
      logical lfnhr
      real colat1
      real(kind=8) t1,t2,t3,t4,t5,ta,tb,tc,td,te,tf,rtc,tx,ty
      real timesum
!!
      real(kind=kind_evod) sl(levs), si(levp1)
!!
      integer              lats_nodes_r(nodes)
!!
      integer              ierr,j,k,l,lenrec,locl,n,node
      integer nosfc,noflx,noaer,nfill
      character*16 cosfc
      data timesum/0./
!!
!!
      character CFHOUR*40,CFORM*40
      integer jdate(4),ndigyr,ndig,kh,IOPROC
!!
      REAL (KIND=KIND_IO8) GESHEM(LONR,LATS_NODE_R)
      INTEGER              GLOBAL_LATS_R(LATR),   lonsperlar(LATR)
!
      REAL (KIND=kind_io8) fluxr(nfxr,LONR,LATS_NODE_R)
!!
      real(kind=kind_evod) secphy,secswr,seclwr
      real(kind=8) tba,tbb,tbc,tbd
      integer iret
!
!     print *,' in wrtout_phyiscs me=',me
      t3=rtc()
      call mpi_barrier(mpi_comm_all,ierr)
      t4=rtc()
      tba=t4-t3
      if(nodes_comp .lt. 1 .or. nodes_comp .gt. nodes) then
        print *, '  NODES_COMP UNDEFINED, CANNOT DO I.O '
        call mpi_finalize()
         stop 333
      endif
!
      ioproc=nodes_comp-1
       
      t1=rtc()
!!
!!
!     CREATE CFHOUR
      JDATE=IDATE
      ndigyr=4
      IF(NDIGYR.EQ.2) THEN
        JDATE(4)=MOD(IDATE(4)-1,100)+1
      ENDIF

!sela set lfnhr to false for writing one step output etc.
      lfnhr=.true.    ! no output
      lfnhr=3600*abs(fhour-nint(fhour)).le.1.or.phour.eq.0
      IF(LFNHR) THEN
        KH=NINT(FHOUR)
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
      ELSE
        KS=NINT(FHOUR*3600)
        KH=KS/3600
        KM=(KS-KH*3600)/60
        KS=KS-KH*3600-KM*60
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,
     &      '("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH,':',KM,':',KS
      ENDIF
      IF(nfill(ens_nam) == 0) THEN
      CFHOUR = CFHOUR(1:nfill(CFHOUR))
      ELSE
      CFHOUR = CFHOUR(1:nfill(CFHOUR)) // ens_nam(1:nfill(ens_nam))
      END IF
!jfe
      nosfc=62
      noflx=63
      noaer=65               ! for g2d_fld
!!
      t3=rtc()
      call MPI_BARRIER(mpi_comm_all,ierr)
      t4=rtc()
      tbd=t4-t3
      t3=rtc()
      SECPHY=(FHOUR-ZHOUR)*3600.
      SECSWR=MAX(SECPHY,FHSWR*3600.)
      SECLWR=MAX(SECPHY,FHLWR*3600.)
!
!*** BUILD STATE ON EACH NODE ********
! build state on each node.   COMP tasks only
! assemble spectral state first then sfc state,
! then (only if liope)  flux state.
! finally (only if gocart is turned on) aer_diag state 
! 
!
      print *,'---- start sfc collection section -----'
      t3=rtc()
      if(mc_comp .ne. MPI_COMM_NULL) then
        print *,' wrtout_physics call sfc_collect '
        CALL sfc_collect(sfc_fld,global_lats_r,lonsperlar)
!
! collect flux grids as was done with sfc grids above.
! but only if liope is true.  If liope is false,
! the fluxes are handled by the original wrtsfc
! predating the I/O task updates.
!
            call   wrtflx_a
     &             (IOPROC,noflx,ZHOUR,FHOUR,IDATE,colat1,SECSWR,SECLWR,
     &              sfc_fld, flx_fld, fluxr, global_lats_r,lonsperlar)

            call   wrtaer
     &             (IOPROC,noaer,ZHOUR,FHOUR,IDATE,
     &              sfc_fld, g2d_fld, global_lats_r, lonsperlar)

      endif                 ! comp node
      t4=rtc()
      td=t4-t3
!
!  done with state build
!  NOW STATE IS ASSEMBLED ON EACH NODE.  GET EVERYTHING OFF THE COMPUTE
!  NODES (currently done with a send to the I/O task_
!  send state to I/O task.  All tasks
!
!jw      if(.not.quilting) then
!jw          print *,'---- start sfc.f section -----'
!jw          call sfc_only_move(ioproc)
!jw          cosfc='SFC.F'//CFHOUR
!jw          call sfc_wrt(ioproc,cosfc,fhour,jdate
!jw     &,                global_lats_r,lonsperlar)
!
!jw          print *,' wrtout_physics call wrtsfc to write out flx'
!jw          call FLX_ONLY_MOVE(ioproc)
!jw          cosfc='FLX.F'//CFHOUR
!jw          call  flx_wrt
!jw     &          (IOPROC,cosfc,ZHOUR,FHOUR,IDATE,
!jw     &           global_lats_r,lonsperlar)
!jw      endif          !  quilting
!
      t4=rtc()
      te=t4-t3
!
!jw      print *,'---- start diag3d.f section -----'
!jw        IF (LDIAG3D) THEN
!jw          print *,' wrtout_physics ldiag3d on so wrt3d '
!jw          no3d=64
!jw          if(icolor.eq.2.and.me.eq.IOPROC)
!jw     &    call BAOPENWT(NO3D,'D3D.F'//CFHOUR,iostat)
!jw          if (hybrid .or. gen_coord_hybrid) then
!     print *,' pl_coeff bef call wrt3d_hyb=',pl_coeff
!jw            call WRT3D_hyb(IOPROC,no3d,nblck,ZHOUR,FHOUR,IDATE,colat1,
!jw     .                     global_lats_r,lonsperlar,pl_coeff,
!jw     &                     SECSWR,SECLWR,sfc_fld%slmsk,flx_fld%psurf)
!jw          else
!jw            call WRT3D(IOPROC,no3d,nblck,ZHOUR,FHOUR,IDATE,colat1,
!jw     .                 global_lats_r,lonsperlar,pl_coeff,
!jw     &                 SECSWR,SECLWR,sl,si,sfc_fld%slmsk,flx_fld%psurf)
!jw          endif
!jw        ENDIF
!
!      if(me .eq. ioproc)  call wrtlog_physics(phour,fhour,idate)
      tb=rtc()
      tf=tb-ta
      t2=rtc()
      print 1011,tf
 1011 format(' WRTOUT_PHYSICS TIME ',f10.4)
      timesum=timesum+(t2-t1)
!     print 1012,timesum,t2-t1,td,te,tf,t4-t3,tba,tbb,tbc,tbd
 1012 format(
     1 ' WRTOUT_PHYSICS TIME ALL TASKS  ',f10.4,f10.4,
     1 ' state, send, io  iobarr, (beginbarr),
     1 spectbarr,open, openbarr )  ' ,
     1  8f9.4)
!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE wrtout_restart_physics(
     &        sfc_fld, fhour,idate,
     &        lats_nodes_r,global_lats_r,lonsperlar,
     &        phy_f3d, phy_f2d, ngptc, nblck, ens_nam)
!!
      use resol_def,               ONLY: latr, levp1, levs, lonr,
     &                                   num_p2d, num_p3d
      use layout1,                 ONLY: me, nodes, lats_node_r
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data, Flx_Var_Data
      USE machine,                 ONLY: kind_evod, kind_phys
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      real(kind=kind_evod) fhour
      character (len=*)  :: ens_nam
!!
      integer              idate(4), ixgr
!
      integer              ngptc, nblck
      REAL (KIND=KIND_phys)
     &            phy_f3d(ngptc,levs,nblck,LATS_NODE_R,num_p3d)
     &,           phy_f2d(LONR,LATS_NODE_R,num_p2d)
!!
      real(kind=kind_evod) sl(levs)
      real(kind=kind_evod) si(levp1)
!!
      integer igen
!!
      INTEGER              lats_nodes_r(nodes)
      INTEGER              GLOBAL_LATS_R(LATR)
      INTEGER              lonsperlar(LATR)
      integer IOPROC, IPRINT
      integer needoro, iret, nfill
!
!!
      integer n3,n4,nflop
      character*20 cfile
      integer nn
!!
      IPRINT = 0
!
      cfile='SFCR'
      print *,' cfile=',cfile,'ens_nam=',ens_nam(1:nfill(ens_nam))
!
      print *,' in rest write fhour=',fhour,
     &  'idate=',idate,' before para_fixio_w'
!
      IOPROC=nodes-1
      CALL para_fixio_w(ioproc,sfc_fld,cfile,fhour,idate,
     &  lats_nodes_r,global_lats_r,lonsperlar,
     &  phy_f3d, phy_f2d, ngptc, nblck, ens_nam)
      print *,'leave restart'
!
      return
      end
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE wrtlog_physics(phour,fhour,idate)
      use namelist_physics_def, ONLY: ens_nam
      implicit none

      integer idate(4),ndigyr,nolog
      integer ks,kh,km,ndig,nfill
      character CFHOUR*40,CFORM*40
      logical lfnhr
      real phour,fhour
!
!     CREATE CFHOUR

!sela set lfnhr to false for writing one step output etc.
      lfnhr=.true.    ! no output
!!mr  lfnhr=.false.   !    output
      lfnhr=3600*abs(fhour-nint(fhour)).le.1.or.phour.eq.0
      IF(LFNHR) THEN
        KH=NINT(FHOUR)
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
        WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH
      ELSE
        KS=NINT(FHOUR*3600)
        KH=KS/3600
        KM=(KS-KH*3600)/60
        KS=KS-KH*3600-KM*60
        NDIG=MAX(LOG10(KH+0.5)+1.,2.)
        WRITE(CFORM,
     &      '("(I",I1,".",I1,",A1,I2.2,A1,I2.2)")') NDIG,NDIG
        WRITE(CFHOUR,CFORM) KH,':',KM,':',KS
      ENDIF
      IF(nfill(ens_nam) == 0) THEN
      CFHOUR = CFHOUR(1:nfill(CFHOUR))
      ELSE
      CFHOUR = CFHOUR(1:nfill(CFHOUR)) // ens_nam(1:nfill(ens_nam))
      END IF

      nolog=99
      OPEN(NOlog,FILE='LOG.F'//CFHOUR,FORM='FORMATTED')
      write(nolog,100)fhour,idate
100   format(' completed mrf fhour=',f10.3,2x,4(i4,2x))
      CLOSE(NOlog)

      RETURN
      END



      SUBROUTINE sfc_collect (sfc_fld,global_lats_r,lonsperlar)
!!
      use resol_def,               ONLY: latr, lonr, ngrids_sfcc, 
     &                                   ngrids_sfcc2d,ngrids_sfcc3d,
     &                                   ngrids_flx, lsoil
      use mod_state,               ONLY:
     &                                   buff_mult_piecea2d,ngrid2d,
     &                                   buff_mult_piecea3d,ngrid3d
      use layout1,                 ONLY: lats_node_r,lats_node_r_max
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      USE machine,                 ONLY: kind_io8, kind_io4
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
!
      INTEGER              GLOBAL_LATS_R(latr)
      INTEGER              lonsperlar(latr)
!!
!!!   real(kind=kind_io4) buff4(lonr,latr,4),bufsave(lonr,lats_node_r)
      real(kind=kind_io8) buffo(lonr,lats_node_r)
      real(kind=kind_io8) buffi(lonr,lats_node_r)
      integer kmsk(lonr,lats_node_r),kmskcv(lonr,lats_node_r)
      integer k,il
       integer ubound
       integer icount
        integer  ierr
!!
      CHARACTER*8 labfix(4)
      real(kind=kind_io4) yhour
      integer,save:: version
      data version/200004/
      data  icount/0/
      integer maxlats_comp
!
      ngrid2d=1
      ngrid3d=1
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!!
!jw      if(allocated(buff_mult_piece)) then
!jw         continue
!jw      else
!jw         allocate(buff_mult_piece(lonr,ngrids_sfcc,lats_node_r))
!jw         allocate(buff_mult_piecef(lonr,0:ngrids_flx,lats_node_r))
!jw         allocate
!jw     1 (buff_mult_piecea(lonr,1:ngrids_flx+ngrids_sfcc+1,lats_node_r))
!jw      endif
      if(allocated(buff_mult_piecea2d)) then
         continue
      else
         allocate
     1 (buff_mult_piecea2d(lonr,lats_node_r_max,1:ngrids_sfcc2d+1),
     1  buff_mult_piecea3d(lonr,lats_node_r_max,1:ngrids_sfcc3d+1))
      endif
!
      kmsk= nint(sfc_fld%slmsk)
!
      ngrid2d=1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%tsea,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
! ngrid=2 here
                                                                                                        
!
      ngrid3d=0
      DO k=1,LSOIL
        buffi(:,:) = sfc_fld%SMC(k,:,:)
        ngrid3d=ngrid3d+1
        CALL uninterprez(1,kmsk,buffo,buffi,global_lats_r,lonsperlar,
     &        buff_mult_piecea3d(1,1,ngrid3d))
      ENDDO
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SHELEG,
     &   global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      DO k=1,LSOIL
        buffi(:,:) = sfc_fld%STC(k,:,:)
!
        ngrid3d=ngrid3d+1
        CALL uninterprez(1,kmsk,buffo,buffi,global_lats_r,lonsperlar,
     &         buff_mult_piecea3d(1,1,ngrid3d))
      ENDDO
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%TG3,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ZORL,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!!
!     where(CV.gt.0.)
!         kmskcv=1
!     elsewhere
!         kmskcv=0
!     endwhere
!
!*********************************************************************
!   Not in version 200501
!     CALL uninterprez(1,kmskcv,buffo,CV,global_lats_r,lonsperlar)
!     CALL uninterprez(1,kmskcv,buffo,CVB,global_lats_r,lonsperlar)
!     CALL uninterprez(1,kmskcv,buffo,CVT,global_lats_r,lonsperlar)
!*********************************************************************
!jws
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALVSF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALVWF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALNSF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALNWF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SLMSK,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%VFRAC,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%CANOPY,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%F10M,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
! T2M
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%T2M,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
! Q2M
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%Q2M,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%VTYPE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%STYPE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))

!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FACSF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FACWF,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%UUSTAR,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FFMM,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FFHH,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
!c-- XW: FOR SEA-ICE Nov04
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%HICE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FICE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%TISFC,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!      write(0,*)'in wrt phy,tisfc=',
!     &    maxval(buff_mult_piecea2d(:,:,ngrid2d)),
!     &    minval(buff_mult_piecea2d(:,:,ngrid2d))

!c-- XW: END SEA-ICE Nov04
!
!lu: the addition of 8 Noah-related records starts here ........................
!tprcp
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%TPRCP,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!srflag
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SRFLAG,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!snwdph
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SNWDPH,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!slc
!      write(0,*)'in wrt phy, before stc,ngrid2d=',ngrid2d,'ngrid3d=',
!     &   ngrid3d,'lsoil=',lsoil
      DO k=1,LSOIL
        buffi(:,:) = sfc_fld%SLC(k,:,:)
        ngrid3d=ngrid3d+1
        CALL uninterprez(1,kmsk,buffo,buffi,global_lats_r,lonsperlar,
     &         buff_mult_piecea3d(1,1,ngrid3d))
      ENDDO
!shdmin
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SHDMIN,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!shdmax
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SHDMAX,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!slope
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SLOPE,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!snoalb
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SNOALB,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!lu: the addition of 8 Noah records ends here .........................
!
! Oro
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ORO,
     &       global_lats_r,lonsperlar,buff_mult_piecea2d(1,1,ngrid2d))
!
!      write(0,*)' finished sfc_collect for  ngrid2d=',ngrid2d,ngrid3d
  999 continue
      return
      end
       subroutine sfc_only_move(ioproc)
!
!***********************************************************************
!
      use resol_def, ONLY: ngrids_flx, ngrids_sfcc, lonr,latr
     &                    ,ngrids_sfcc2d,ngrids_sfcc3d
      use mod_state, ONLY: buff_mult_pieces,buff_mult_piece,
     &                     buff_mult_piecea2d,
     &                     buff_mult_piecea3d, 
     &                     ivar_global_a, ivar_global
      use layout1,   ONLY: nodes, ipt_lats_node_r, lats_node_r, 
     &                     lats_node_r_max, me, nodes_comp
      use mpi_def,   ONLY: mpi_comm_null, mpi_r_io, mc_comp, 
     &                     mpi_integer, mpi_comm_all, liope, 
     &                     info, stat
      implicit none
!
      integer ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
!     integer lats_nodes_r(nodes),ipt,maxfld,ioproc,nproct
      integer ioproc
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer illen,ubound,nd1
      integer icount
      data icount/0/
      save icount
      integer maxlats_comp
!  allocate the data structures
!
      if(icount .eq. 0) then
         allocate(ivar_global(10))
         allocate(ivar_global_a(10,nodes))
         ivar_global(1)=ipt_lats_node_r
         ivar_global(2)= lats_node_r
         ivar_global(3)=lats_node_r_max
         call mpi_gather(ivar_global,10,MPI_INTEGER,
     1       ivar_global_a,10,MPI_INTEGER,ioproc,mc_comp,ierr)
         if(me==ioproc) write(0,*)'in sfc_only_move, ivar_global_a=',
     &     ivar_global_a(1:3,1:nodes)
         icount=icount+1
      endif
!!
      if(allocated(buff_mult_pieces)) then
          deallocate(buff_mult_pieces)
      else
          maxlats_comp=lats_node_r_max
          if(me .eq. ioproc) then
            maxlats_comp=ivar_global_a(3,1)
           endif
      endif
      if(me .eq. ioproc) then
!gwv watch this!!
          allocate
     1  (buff_mult_pieces(lonr*latr*ngrids_sfcc))
         buff_mult_pieces=0.
      endif

      if(allocated(buff_mult_piece)) then
         continue
      else
         allocate(buff_mult_piece(lonr*lats_node_r*ngrids_sfcc))
      endif                                                   
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!   SENDLOOP of grids from comp processors to I/O task.  The
!   I/O task may or may not be a comp task also.  The
!   send logic on that task is different for these two cases
!
!  big send
!     if(me .gt. -1) return
!
       buff_mult_piece=0.
       buff_mult_piece(1:lonr*lats_node_r*ngrids_sfcc2d)=
     1 reshape(buff_mult_piecea2d(1:lonr,1:lats_node_r,1:ngrids_sfcc2d),
     1   (/lonr*lats_node_r*ngrids_sfcc2d/)) 
       buff_mult_piece(lonr*lats_node_r*ngrids_sfcc2d+1:
     1    lonr*lats_node_r*ngrids_sfcc)=
     1 reshape(buff_mult_piecea3d(1:lonr,1:lats_node_r,1:ngrids_sfcc3d),
     1   (/lonr*lats_node_r*ngrids_sfcc3d/) )
!
      IF (ME .ne. ioproc) THEN    !   Sending the data
         msgtag=me
         illen=lats_node_r
         CALL mpi_send            !  send the local grid domain
     &(buff_mult_piece,illen*lonr*ngrids_sfcc,MPI_R_IO,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
      ELSE
        if( MC_COMP .ne. MPI_COMM_NULL) then
!
c iotask is also a compute task.  send is replaced with direct
c  array copy
!
         if(nodes_comp==1) then
           buff_mult_pieces(1:lonr*lats_node_r*ngrids_sfcc)=
     1     buff_mult_piece(1:lonr*lats_node_r*ngrids_sfcc)
!                              END COMPUTE TASKS PORTION OF LOGIC
         else
!
!  END COMPUTE TASKS PORTION OF LOGIC
!  receiving part of I/O task, ioproc is the last fcst pe
!
!!
!!      for pes ioproc
!jw        nd1=lonr*lats_node_r*ngrids_sfcc
        nd1=0
        DO proc=1,nodes_comp
          illen=ivar_global_a(2,proc)
          if (proc.ne.ioproc+1) then
            msgtag=proc-1
            CALL mpi_recv(buff_mult_pieces(nd1+1),
     1        illen*lonr*ngrids_sfcc
     1        ,MPI_R_IO,proc-1,
     &                msgtag,MPI_COMM_ALL,stat,info)
          else
           buff_mult_pieces(nd1+1:nd1+lonr*illen*ngrids_sfcc)=
     1       buff_mult_piece(1:lonr*illen*ngrids_sfcc)
          endif
          nd1=nd1+illen*lonr*ngrids_sfcc
        enddo
        endif

       Endif
!end ioproc
      ENDIF
!!
      return
      end
      SUBROUTINE sfc_wrt(IOPROC,cfile,xhour,idate
     &,                  global_lats_r,lonsperlar)
!!
      use module_nemsio
      use resol_def,    ONLY: lonr, latr, levs,ngrids_sfcc,
     &   ncld,ntrac,ntcw,ntoz,lsoil, ivssfc,thermodyn_id,sfcpress_id
      use layout1,      ONLY: me,idrt
      USE machine,      ONLY: kind_io8, kind_io4
!jw
      use gfs_physics_output, only : PHY_INT_STATE_ISCALAR,
     &    PHY_INT_STATE_RSCALAR,
     &    PHY_INT_STATE_1D_I,PHY_INT_STATE_1D_R,
     &    PHY_INT_STATE_2D_R_SFC,PHY_INT_STATE_3D_R
      implicit none
!!
      integer IOPROC
      character*16 cfile
      real(kind=kind_io8) xhour
      integer idate(4),k,il, ngridss
!
      integer i,j,ndim3,N2DR,idate7(7),nrec,kount
      integer nfhour,nfminute,nfsecondd,nfsecondn
      logical  :: outtest
      integer ::nmetavari,nmetavarr,nmetavarl,nmetaaryi,nmetaaryr
      character(16),allocatable :: recname(:),reclevtyp(:)
      integer,allocatable :: reclev(:)
      character(16),allocatable :: variname(:),varrname(:),
     &    aryiname(:),aryrname(:)
      integer,allocatable :: varival(:),aryilen(:),aryival(:,:)
      real(kind_io4),allocatable    :: varrval(:),aryrval(:,:)
      real(kind_io4),allocatable :: buff_mult(:,:,:),tmp(:)
      type(nemsio_gfile) gfileout
!
!!
      CHARACTER*8 labfix(4)
      real(kind=kind_io4) yhour
      integer,save:: version
      data version/200501/
      INTEGER              GLOBAL_LATS_R(latr), lonsperlar(latr)
!
      integer iret
      logical first
      save first
      save  recname, reclevtyp, reclev
      save nrec,nmetavari,nmetavarr,nmetaaryi,nmetaaryr,
     &     variname,varrname,aryiname,
     &     varival,varrval,aryilen,aryival
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build surface fields in to buff_mult
!
      if (me.eq.ioproc) then
!
!
         allocate(buff_mult(lonr,latr,ngrids_sfcc))
         do ngridss=1,ngrids_sfcc
           call unsplit2z(ioproc,ngridss,ngrids_sfcc,
     &       buff_mult(1,1,ngridss),global_lats_r)
         enddo
!
!    Building surface field is done
!
         if (first) then
!write out nemsio sfc file:
          nrec=ngrids_sfcc
          kount=size(PHY_INT_STATE_ISCALAR,2)
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_SFC')
     &        nmetavari=nmetavari+1
          enddo
          allocate(variname(nmetavari),varival(nmetavari))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY' .or.
     &      trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_SFC' )then
            j=j+1
            variname(j)=trim(PHY_INT_STATE_ISCALAR(1,i))
            if(trim(variname(j))=='latr') varival(j)=latr
            if(trim(variname(j))=='lonr') varival(j)=lonr
            if(trim(variname(j))=='levs') varival(j)=levs
            if(trim(variname(j))=='ntoz') varival(j)=ntoz
            if(trim(variname(j))=='ntcw') varival(j)=ntcw
            if(trim(variname(j))=='ncld') varival(j)=ncld
            if(trim(variname(j))=='ntrac') varival(j)=ntrac
            if(trim(variname(j))=='thermodyn_id')varival(j)=thermodyn_id
            if(trim(variname(j))=='sfcpress_id') varival(j)=sfcpress_id
            if(trim(variname(j))=='lsoil') varival(j)=lsoil
            if(trim(variname(j))=='idrt') varival(j)=idrt
            if(trim(variname(j))=='ivssfc') varival(j)=ivssfc
           endif
          enddo
!!for real var::
          nmetavarr=0
          do i=1,kount
           if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_SFC')
     &     nmetavarr=nmetavarr+1
          enddo
          allocate(varrname(nmetavarr),varrval(nmetavarr))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_SFC')then
             j=j+1
             varrname(j)=trim(PHY_INT_STATE_RSCALAR(1,i))
             if(trim(varrname(j))=='fhour') varrval(j)=xhour
           endif
          enddo
!!for 1D ary::
          nmetaaryi=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_SFC')
     &     nmetaaryi=nmetaaryi+1
          enddo
          allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_SFC')then
             j=j+1
             aryiname(j)=trim(PHY_INT_STATE_1D_I(1,i))
             if(aryiname(j)=='IDATE') aryilen(j)=size(idate)
           endif
          enddo
          allocate(aryival(maxval(aryilen),nmetaaryi) )
          aryival(1:aryilen(1),1)=idate(1:aryilen(1))
!
!!for record name, levtyp and lev
          allocate (recname(nrec),reclevtyp(nrec),reclev(nrec))
          N2DR=0
          do i=1,kount
           if(trim(PHY_INT_STATE_2D_R_SFC(2,i)).eq.'OGFS_SFC')then
            N2DR=N2DR+1
            recname(N2DR)=trim(PHY_INT_STATE_2D_R_SFC(1,i))
            reclevtyp(N2DR)=trim(trim(PHY_INT_STATE_2D_R_SFC(3,i)))
            reclev(N2DR)=1
           endif
          enddo
!
          do i=1,kount
           if(trim(PHY_INT_STATE_3D_R(2,i)).eq.'OGFS_SFC')then
            ndim3=0
            if(trim(PHY_INT_STATE_3D_R(4,i)).eq.'lsoil') then
             ndim3=lsoil
            endif
            if(ndim3>0) then
             do j=1,ndim3
              N2DR=N2DR+1
              recname(N2DR)=trim(PHY_INT_STATE_3D_R(1,i))
              reclevtyp(N2DR)=trim(trim(PHY_INT_STATE_3D_R(3,i)) )
              if(trim(PHY_INT_STATE_3D_R(4,i)).eq.'lsoil') then
                reclev(N2DR)=j
              endif
             enddo
            endif
!
           endif
          enddo
!end first
          first=.false.
         endif
     
        idate7=0
        idate7(1)=idate(4)
        idate7(2)=idate(2)
        idate7(3)=idate(3)
        idate7(4)=idate(1)
        idate7(7)=100           !denominator for second
!
        nfhour=int(xhour)
        nfminute=int((xhour-nfhour)*60)
        nfsecondn=int(((xhour-nfhour)*3600-nfminute*60)*100)
        nfsecondd=100
!
        call nemsio_init()
!
        call nemsio_open(gfileout,trim(cfile),'write',
     &    iret = iret,
     &    modelname='GFS',gdatatype='bin4',
     &    idate=idate7,nrec=nrec,
     &    dimx=lonr,dimy=latr,dimz=levs,ncldt=ncld,nmeta=5,
     &    nfhour=nfhour,nfminute=nfminute,nfsecondn=nfsecondn,
     &    nfsecondd=nfsecondd,
     &    extrameta=.true.,nmetavari=nmetavari,
     &    nmetavarr=nmetavarr,
     &    nmetaaryi=nmetaaryi,nmetaaryr=nmetaaryr,
     &    variname=variname,varival=varival,varrname=varrname,
     &    varrval=varrval,
     &    aryiname=aryiname,aryilen=aryilen,aryival=aryival,
     &    ntrac=ntrac,nsoil=lsoil,idrt=idrt,
     &    recname=recname,reclevtyp=reclevtyp,reclev=reclev)
!
        allocate(tmp(lonr*latr))
        do i=1,nrec
         tmp(:)=reshape(buff_mult(:,:,i),(/lonr*latr/) )
         call nemsio_writerec(gfileout,i,tmp,iret=iret)
        enddo
        deallocate(tmp)
        deallocate(buff_mult)
!
        call nemsio_close(gfileout)
!end write pe
      endif
!
      return
      end
      SUBROUTINE wrtflx_a(IOPROC,noflx,ZHOUR,FHOUR,IDATE,colat1,
     &                  SECSWR,SECLWR, sfc_fld, flx_fld, fluxr,
     &                  global_lats_r,lonsperlar)
!!
      use resol_def,               ONLY: lonr, latr, levp1, lsoil, nfxr,
     *                                   ngrids_sfcc
      use mod_state,               ONLY: buff_mult_piecef
      use layout1,                 ONLY: me, lats_node_r
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data, Flx_Var_Data
      USE machine,             ONLY: kind_io8, kind_io4,grib_undef
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      INTEGER              GLOBAL_LATS_R(LATR)
      INTEGER              lonsperlar(LATR)
      integer   IOPROC
!!
      integer LEN,NFLD
      integer j,i,k,itop,ibot,k4,l,noflx,nundef,ngrid2d
!*    PARAMETER(NFLD=18)
      PARAMETER(NFLD=18+6)      ! 550nm AOD added
       integer ilpds,iyr,imo,ida,ihr,ifhr,ithr,lg,ierr
       real (kind=kind_io8) RTIMER(NFLD),rtime,rtimsw,rtimlw
       real (kind=kind_io8) colat1
       real (kind=kind_io8) cl1,secswr,zhour,fhour,seclwr
C

      real(kind=kind_io4) wrkga(lonr*latr),wrkgb(lonr*latr)
      real(kind=kind_io8) slmskful(lonr,lats_node_r)
      real(kind=kind_io8) slmskloc(LONR,LATS_NODE_R)
!
      INTEGER     IDATE(4), IDS(255),IENS(5)
!     INTEGER     IDATE(4)
      real (kind=kind_io8) SI(LEVP1)
!
!sela..................................................................
!* change rflux 3rd dimension from 27 to nfxr (Sarah Lu)
!*    real (kind=kind_io8)   rflux(lonr,LATS_NODE_R,27)     
      real (kind=kind_io8)   rflux(lonr,LATS_NODE_R,nfxr)
      real (kind=kind_io8)   glolal(lonr,LATS_NODE_R)
      real (kind=kind_io8)   buffo(lonr,LATS_NODE_R)
      real (kind=kind_io4)   buff1(lonr,latr)
      real (kind=kind_io4)   buff1l(lonr*latr)
!sela..................................................................
      real (kind=kind_io8)  FLUXR(nfxr,LONR,LATS_NODE_R)
!sela..................................................................
      integer kmsk(lonr,lats_node_r),kmsk0(lonr,lats_node_r)
      integer kmskcv(lonr,LATS_NODE_R)
!jws
      integer kmskgrib(lonr,lats_node_r)
      real(kind=kind_io4) buff_max
!jwe
!
!!
      kmsk=nint(sfc_fld%slmsk)
      kmsk0=0
!
      kmskgrib=0
      ngrid2d=1
!
      CALL uninterprez(1,kmsk,glolal,sfc_fld%slmsk,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
      slmskloc=glolal
      slmskful=buff_mult_piecef(1:lonr,1:lats_node_r,ngrid2d)
c
      do k=1,nfxr
       do j=1,LATS_NODE_R
        do i=1,lonr
         rflux(i,j,k)=fluxr(k,i,j)
        enddo
       enddo
      enddo
!!
!
      IF(FHOUR.GT.ZHOUR) THEN
        RTIME=1./(3600.*(FHOUR-ZHOUR))
      ELSE
        RTIME=0.
      ENDIF
      IF(SECSWR.GT.0.) THEN
        RTIMSW=1./SECSWR
      ELSE
        RTIMSW=1.
      ENDIF
      IF(SECLWR.GT.0.) THEN
        RTIMLW=1./SECLWR
      ELSE
        RTIMLW=1.
      ENDIF
      RTIMER=RTIMSW
      RTIMER(1)=RTIMLW
      CL1=colat1
!!
!..........................................................
      glolal=flx_fld%DUSFC*RTIME
!jw
      ngrid2d=1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &    buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '01)Zonal compt of momentum flux (N/m**2) land and sea surface '

!..........................................................
      glolal=flx_fld%DVSFC*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '02)Merid compt of momentum flux (N/m**2) land and sea surface '
!..........................................................
      glolal=flx_fld%DTSFC*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '03)Sensible heat flux (W/m**2) land and sea surface           '
!..........................................................
      glolal=flx_fld%DQSFC*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '04)Latent heat flux (W/m**2) land and sea surface             '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%tsea,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribsn ierr=',ierr,'  ',
!    x '05)Temperature (K) land and sea surface                       '
!..........................................................
      glolal(:,:) = sfc_fld%SMC(1,:,:)
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribsn ierr=',ierr,'  ',
!    x '06)Volumetric soil moist content (frac) layer 10cm and 0cm    '
!..........................................................
      glolal(:,:) = sfc_fld%SMC(2,:,:)
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!lu  x '07)Volumetric soil moist content (frac) layer 200cm and 10cm  '
!    + '07)Volumetric soil moist content (frac) layer 40cm and 10cm  '
!..........................................................
      glolal(:,:) = sfc_fld%STC(1,:,:)
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
      nundef=0
      buff_max=0.
      do j=1,lats_node_r
      do i=1,lonr
       if(buff_mult_piecef(i,j,ngrid2d)/=grib_undef) then
         if(buff_mult_piecef(i,j,ngrid2d) >buff_max)
     &   buff_max=buff_mult_piecef(i,j,ngrid2d)
         nundef=nundef+1
       endif
      enddo
      enddo
!      write(0,*)'in wrtsfc_a, max stc=',buff_max,' grib_undef=',
!     &   grib_undef,'nundef=',nundef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '08)Temp (K) layer betw two depth below land sfc 10cm and 0cm  '
!..........................................................
      glolal(:,:) = sfc_fld%STC(2,:,:)
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(slmskful/=1._kind_io8)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!lu  x '09)Temp (K) layer betw two depth below land sfc 200cm and 10cm'
!    + '09)Temp (K) layer betw two depth below land sfc 40cm and 10cm'
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,sfc_fld%sheleg,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '10)Water equiv of accum snow depth (kg/m**2) land sea surface '
c..........................................................
!      write(0,*)'before DLWSFC'
      glolal = flx_fld%DLWSFC*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '11)Downward long wave radiation flux (W/m**2) land sea surface'
!..........................................................
      glolal = flx_fld%ULWSFC*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '12)Upward long wave radiation flux (W/m**2) land sea surface  '
!..........................................................
!.......  FIX FLUXES FOR APPROX DIURNAL CYCLE
      DO 113 K=1,4
       do j=1,LATS_NODE_R
        do i=1,lonr
         glolal(i,j) = rflux(i,j,k)*RTIMER(k)
        enddo
       enddo
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0.and.k.eq.1)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '13)Upward long wave radiation flux (W/m**2) top of atmosphere '
!     if(ierr.ne.0.and.k.eq.2)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '14)Upward solar radiation flux (W/m**2) top of atmosphere     '
!     if(ierr.ne.0.and.k.eq.3)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '15)Upward solar radiation flux (W/m**2) land and sea surface  '
!     if(ierr.ne.0.and.k.eq.4)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '16)Downward solar radiation flux (W/m**2) land and sea surface'
  113 CONTINUE
!..........................................................
!
!     For UV-B fluxes
!
      do j=1,LATS_NODE_R
        do i=1,lonr
          glolal(i,j) = rflux(i,j,21)*rtimsw
        enddo
      enddo
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '17)UV-B Downward solar flux (W/m**2) land sea surface'
      do j=1,LATS_NODE_R
        do i=1,lonr
          glolal(i,j) = rflux(i,j,22)*rtimsw
        enddo
      enddo
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '18)clear sky UV-B Downward solar flux (W/m**2) land sea surface'
!
!     End UV-B fluxes
!
!..........................................................
!..........................................................
      DO 813 K=5,7
!
       do j=1,LATS_NODE_R
        do i=1,lonr
         glolal(i,j) = rflux(i,j,k)*100.*rtimsw
        enddo
       enddo
      where(glolal.ge.0.5)
        kmskcv = 1
      elsewhere
        kmskcv = 0
      endwhere
!!
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!      where(buff_mult_piecef(:,:,ngrid2d)<=0.5_kind_io4)
!     &    buff_mult_piecef(:,:,ngrid2d)=grib_undef
      kmskgrib=0
      where(buff_mult_piecef(:,:,ngrid2d)<=0.5_kind_io4)
     &    kmskgrib=1
!     if(ierr.ne.0.and.k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '19)Total cloud cover (percent) high cloud layer               '
!     if(ierr.ne.0.and.k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '23)Total cloud cover (percent) middle cloud layer             '
!     if(ierr.ne.0.and.k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '27)Total cloud cover (percent) low cloud layer                '
!
        K4=4+(K-5)*4
        L=K4+1
!
       do j=1,LATS_NODE_R
        do i=1,lonr
         if(rflux(i,j,k).gt.0.)then
          glolal(i,j) = rflux(i,j,k+3)*1000./rflux(i,j,k)
         else
          glolal(i,j) = 0.
         endif
        enddo
       enddo
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(kmskgrib==1) buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0.and.k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '20)Pressure (Pa) high cloud top level                         '
!     if(ierr.ne.0.and.k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '24)Pressure (Pa) middle cloud top level                       '
!     if(ierr.ne.0.and.k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '28)Pressure (Pa) low cloud top level                          '
        L=K4+2
!
       do j=1,LATS_NODE_R
        do i=1,lonr
         if(rflux(i,j,k).gt.0.)then
          glolal(i,j) = rflux(i,j,k+6)*1000./rflux(i,j,k)
         else
          glolal(i,j) = 0.
         endif
        enddo
       enddo
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(kmskgrib==1) buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0.and.k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '21)Pressure (Pa) high cloud bottom level                      '
!     if(ierr.ne.0.and.k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '25)Pressure (Pa) middle cloud bottom level                    '
!     if(ierr.ne.0.and.k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '29)Pressure (Pa) low cloud bottom level                       '
        L=K4+3
!
       do j=1,LATS_NODE_R
        do i=1,lonr
         if(rflux(i,j,k).gt.0.)then
          glolal(i,j) = rflux(i,j,k+9)/rflux(i,j,k)
         else
          glolal(i,j) = 0.
         endif
        enddo
       enddo
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(kmskgrib==1) buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0.and.k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '22)Temperature (K) high cloud top level                       '
!     if(ierr.ne.0.and.k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '26)Temperature (K) middle cloud top level                     '
!     if(ierr.ne.0.and.k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '30)Temperature (K) low cloud top level                        '
        L=K4+4
!
  813 CONTINUE
!!
!...................................................................
      glolal = flx_fld%GESHEM*1.E3*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '31)Precipitation rate (kg/m**2/s) land and sea surface        '
c...................................................................
      glolal = flx_fld%BENGSH*1.E3*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '32)Convective precipitation rate (kg/m**2/s) land sea surface '
!...................................................................
      glolal = flx_fld%GFLUX*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(slmskful==0._kind_io8)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '33)Ground heat flux (W/m**2) land and sea surface             '
!...................................................................
!     buffo=MOD(slmskloc,2._kind_io8)
!gwv   add something here
!     do j=1,lats_node_r
!       do i=1,lonr
!         buff_mult_piecea(i,ngrid,j)=buffo(i,j)
!       end do
!     end do
!     ngrid=ngrid+1
!...................................................................
!     Add land/sea mask here
      ngrid2d=ngrid2d+1
      buffo=MOD(slmskloc,2._kind_io8)
      do j=1,lats_node_r
        do i=1,lonr
!jw          buff_mult_piecea(i,ngrid,j)=buffo(i,j)
          buff_mult_piecef(i,j,ngrid2d)=buffo(i,j)
        end do
      end do
!jw        ngrid=ngrid+1
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '34)Land-sea mask (1=land; 0=sea) (integer) land sea surface   '
!gwv   add something here
!
!c-- XW: FOR SEA-ICE Nov04
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%fice,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '35)Ice concentration (ice>0; no ice=0) (1/0) land sea surface '
!c-- XW: END SEA-ICE
!...................................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%u10m,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '36)u wind (m/s) height above ground                           '
!...................................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%v10m,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '37)v wind (m/s) height above ground                           '
!...................................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%t2m,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '38)Temperature (K) height above ground                        '
!...................................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%q2m,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '39)Specific humidity (kg/kg) height above ground              '
!...................................................................
      glolal = flx_fld%PSURF*1.E3
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '40)Pressure (Pa) land and sea surface                         '
!...................................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%tmpmax,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '41)Maximum temperature (K) height above ground                '
!...................................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%tmpmin,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '42)Minimum temperature (K) height above ground                '
!...................................................................
      glolal = flx_fld%RUNOFF * 1.E3
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(slmskful==0._kind_io8)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '43)Runoff (kg/m**2) land and sea surface                      '
!...................................................................
      glolal = flx_fld%EP * RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(slmskful==0._kind_io8)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '44)Potential evaporation rate (w/m**/) land and sea surface   '
!...................................................................
      glolal = flx_fld%CLDWRK * RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '45)Cloud work function (J/Kg) total atmospheric column        '
!...................................................................
      glolal = flx_fld%DUGWD*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '46)Zonal gravity wave stress (N/m**2) land and sea surface    '
!...................................................................
      glolal = flx_fld%DVGWD*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '47)Meridional gravity wave stress (N/m**2) land sea surface   '
!...................................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%hpbl,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '48)Boundary layer height '
!...................................................................
!hmhj CALL uninterprez(2,kmsk0,buffo,flx_fld%pwat,
!hmhj&                 global_lats_r,lonsperlar)
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%pwat,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '49)Precipitable water (kg/m**2) total atmospheric column      '
!...................................................................
!
       do j=1,LATS_NODE_R
        do i=1,lonr
         if (rflux(i,j,4).GT.0.) then
          glolal(i,j) = rflux(i,j,3)/rflux(i,j,4) * 100.
          if (glolal(i,j).GT.100.) glolal(i,j) = 100.
         else
          glolal(i,j) = 0.
         endif
        enddo
       enddo
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '50)Albedo (percent) land and sea surface                      '
!
       do j=1,LATS_NODE_R
        do i=1,lonr
         glolal(i,j) = rflux(i,j,26)*100.*rtimsw
        enddo
       enddo
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '51)Total cloud cover (percent) total atmospheric column       '
!
! CONVECTIVE CLOUDS
! LABELED INSTANTANEOUS BUT ACTUALLY AVERAGED OVER FHSWR HOURS
!
      glolal = sfc_fld%CV*1.E2
      where(glolal.ge.0.5)
        kmskcv = 1
      elsewhere
        kmskcv = 0
      endwhere
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      kmskgrib=0
      where(buff_mult_piecef(:,:,ngrid2d)<0.5_kind_io8)
     &     kmskgrib=1
!      where(buff_mult_piecef(:,:,ngrid2d)<0.5_kind_io8)
!     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!      write(0,*)'52after cloud cover,buff=',maxval(buff_mult_piecef(
!     &   1:lonr,1:lats_node_r,ngrid2d)), minval(buff_mult_piecef(
!     &   1:lonr,1:lats_node_r,ngrid2d)),'ngrid2d=',ngrid2d

!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '52)Total cloud cover (percent) convective cloud layer         '
!.................................................
       do j=1,LATS_NODE_R
        do i=1,lonr
        glolal(i,j) = 0.
        IF(sfc_fld%CV(i,j).GT.0.) THEN
!        ITOP=NINT(CVT(i,j))
!        IF(ITOP.GE.1.AND.ITOP.LE.LEVS)
!    &   glolal(i,j)=SI(ITOP+1)*PSURF(i,j)*1.E3
!...      cvt already a pressure (cb)...convert to Pa
         glolal(i,j) = sfc_fld%CVT(i,j)*1.E3
        END IF
       ENDDO
      ENDDO
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(kmskgrib==1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '53)Pressure (Pa) convective cloud top level                   '
!.................................................
       do j=1,LATS_NODE_R
        do i=1,lonr
        glolal(i,j) = 0.
        IF(sfc_fld%CV(i,j).GT.0.) THEN
!        Ibot=NINT(CVB(i,j))
!        IF(Ibot.GE.1.AND.Ibot.LE.LEVS)
!    &   glolal(i,j)=SI(IBOT)*PSURF(i,j)*1.E3
!...      cvb already a pressure (cb)...convert to Pa
         glolal(i,j) = sfc_fld%CVB(i,j)*1.E3
        END IF
       ENDDO
      ENDDO
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(kmskgrib==1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '54)Pressure (Pa) convective cloud bottom level                '
!.................................................
!...   SAVE B.L. CLOUD AMOUNT
!
       do j=1,LATS_NODE_R
        do i=1,lonr
         glolal(i,j) = rflux(i,j,27)*100.*rtimsw
        enddo
       enddo
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '55)Total cloud cover (percent) boundary layer cloud layer     '
!c-- XW: FOR SEA-ICE Nov04
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%hice,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '56)Sea ice thickness (m) category 1'
!c-- XW: END SEA-ICE
!.................................................
!lu: add smc(3:4), stc(3:4), slc(1:4), snwdph, canopy
!lu: addition of 10 records starts here -------------------------------
      if(lsoil.gt.2)then
        glolal(:,:) = sfc_fld%SMC(3,:,:)
      ngrid2d=ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '57)Volumetric soil moist content (frac) layer 100cm and 40cm '
!..........................................................
        glolal(:,:) = sfc_fld%SMC(4,:,:)
      ngrid2d=ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '58)Volumetric soil moist content (frac) layer 200cm and 100cm '
!..........................................................
        glolal(:,:) = sfc_fld%STC(3,:,:)
      ngrid2d=ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '59)Temp (K) layer betw two depth below land sfc 100cm and 40cm'
!..........................................................
        glolal(:,:) = sfc_fld%STC(4,:,:)
      ngrid2d=ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '60)Temp (K) layer betw two depth below land sfc 200cm and 100cm'
      endif
!..........................................................
      glolal(:,:) = sfc_fld%SLC(1,:,:)
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '61)Liquid soil moist content (frac) layer 10cm and 0cm  '
!..........................................................
      glolal(:,:) = sfc_fld%SLC(2,:,:)
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '62)Liquid soil moist content (frac) layer 40cm and 10cm '
!..........................................................
      if(lsoil.gt.2)then
        glolal(:,:) = sfc_fld%SLC(3,:,:)
      ngrid2d=ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '63)Liquid soil moist content (frac) layer 100cm and 40cm'
!..........................................................
        glolal(:,:) = sfc_fld%SLC(4,:,:)
      ngrid2d=ngrid2d+1
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '64)Liquid soil moist content (frac) layer 200cm and 100cm'
      endif
!..........................................................
      glolal = sfc_fld%SNWDPH / 1.E3       !! convert from mm to m
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '65)Snow depth (m) land surface 
c..........................................................
!     LBM=slmskful.EQ.1._kind_io8
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,sfc_fld%canopy,
     &       global_lats_r,lonsperlar,
     &       buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '66)Canopy water content (kg/m^2) land surface      '
!lu: addition of 10 records ends here -------------------------------
!
!wei: addition of 30 records starts here -------------------------------
      glolal = sfc_fld%ZORL / 1.E2       !! convert from cm to m
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '67)Surface roughness (m)       '
!..........................................................
      glolal = sfc_fld%vfrac*100.
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '68)Vegetation fraction (fractional) land surface      '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,glolal,sfc_fld%vtype,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
      buffo=MOD(glolal,2._kind_io8)
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '69)Vegetation type land surface      '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,glolal,sfc_fld%stype,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
      buffo=MOD(glolal,2._kind_io8)
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '70)Soil type land surface      '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,glolal,sfc_fld%slope,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
      buffo=MOD(glolal,2._kind_io8)
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '71)Slope type land surface      '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%uustar,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '72)Frictional velocity (m/s)     '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%oro,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '73)Surface height (m)       '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(1,kmsk,buffo,sfc_fld%srflag,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '74)Freezing precip flag land surface      '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%chh,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '75)Exchange coefficient CH(m/s)       '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%cmm,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '76)Exchange coefficient CM(m/s)       '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,flx_fld%EPI,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '77)Potential evaporation rate (w/m**2) land and sea surface   '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%DLWSFCI,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '78)Downward long wave radiation flux (W/m**2) '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%ULWSFCI,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '79)Upward long wave radiation flux (W/m**2)  '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%USWSFCI,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '80)Upward short wave radiation flux (W/m**2)  '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%DSWSFCI,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '81)Downward short wave radiation flux (W/m**2)   '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%DTSFCI,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '82)Sensible heat flux (W/m**2) land and sea surface       '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%DQSFCI,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '83)Latent heat flux (W/m**2) land and sea surface         '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,flx_fld%GFLUXI,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '84)Ground heat flux (W/m**2) land and sea surface         '
!..........................................................
      glolal = flx_fld%SRUNOFF * 1.E3
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &    buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '85)Surface runoff (kg/m^2) land surface      '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%t1,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '86)Lowest model level Temp (K)      '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%q1,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '87)Lowest model specific humidity (kg/kg)    '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%u1,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '88)Lowest model u wind (m/s)      '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk0,buffo,flx_fld%v1,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '89)Lowest model v wind (m/s)       '
!..........................................................
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,flx_fld%zlvl,
     &       global_lats_r,lonsperlar,buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '90)Lowest model level height (m) land surface      '
!..........................................................
      glolal = flx_fld%EVBSA*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '91)Direct evaporation from bare soil(W/m^2) land surface      '
!..........................................................
      glolal = flx_fld%EVCWA*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '92)Canopy water evaporation(W/m^2) land surface      '
!..........................................................
      glolal = flx_fld%TRANSA*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '93)Transpiration (W/m^2) land surface      '
!..........................................................
      glolal = flx_fld%SBSNOA*RTIME
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '94)Snow Sublimation (W/m^2) land surface      '
!..........................................................
      glolal = flx_fld%SNOWCA*RTIME*100.
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '95)Snow Cover (fraction) land surface      '
!..........................................................
      glolal = flx_fld%soilm*1.E3       !! convert from m to (mm)kg/m^2
      ngrid2d=ngrid2d+1
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
      where(nint(slmskful)/=1)
     &     buff_mult_piecef(:,:,ngrid2d)=grib_undef
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '96)Total column soil moisture (Kg/m^2) land surface      '

Cwei: addition of 30 records ends here -------------------------------

Clu: addition of 6 aod fields starts here ---------------------------
       do k = nfxr-5, nfxr
         do j=1,LATS_NODE_R
           do i=1,lonr
             glolal(i,j) = rflux(i,j,k)*RTIMER(k-9)
           enddo
         enddo
       ngrid2d=ngrid2d+1
       CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar,
     &     buff_mult_piecef(1,1,ngrid2d))
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '97)DU Aerosol optical depth at 550nm land sea surface'
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '98)BC Aerosol optical depth at 550nm land sea surface'
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '99)OC Aerosol optical depth at 550nm land sea surface'
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '100)SU Aerosol optical depth at 550nm land sea surface'
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '101)SS Aerosol optical depth at 550nm land sea surface'
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '102)Total Aerosol optical depth at 550nm land sea surface'
       enddo
Clu: addition of 6 aod fields ends here -----------------------------


      if(me.eq.ioproc)
     &   PRINT *,'(wrtflx_a) GRIB FLUX FILE WRITTEN ',FHOUR,IDATE,noflx
!!
      RETURN
      END

!!*********************************************************************
!! This routine is added to output 2d aerosol diag fields (Sarah Lu)

      SUBROUTINE wrtaer(IOPROC,noaer,ZHOUR,FHOUR,IDATE,
     &             sfc_fld, g2d_fld,global_lats_r, lonsperlar)
!!
      use resol_def,               ONLY: lonr, latr, ngrids_aer
      use mod_state,               ONLY: buff_mult_pieceg
      use layout1,                 ONLY: me, lats_node_r
      use gfs_physics_sfc_flx_mod, ONLY: Sfc_Var_Data
      use gfs_physics_g2d_mod,     ONLY: G2D_Var_Data
      USE machine,                 ONLY: kind_io8, kind_io4
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(G2D_Var_Data)        :: g2d_fld
      INTEGER                   GLOBAL_LATS_R(LATR)
      INTEGER                   lonsperlar(LATR)
      integer                   IOPROC
!!
      integer                   i,j,k,l,noaer,ngrid2d,ierr
      real (kind=kind_io8)      rtime
      real (kind=kind_io8)      zhour,fhour
C

!     real(kind=kind_io8) slmskful(lonr,lats_node_r)
!     real(kind=kind_io8) slmskloc(LONR,LATS_NODE_R)
!
      INTEGER     IDATE(4), IDS(255),IENS(5)
!
      real (kind=kind_io8)   glolal(lonr,LATS_NODE_R)
      real (kind=kind_io8)   buffo(lonr,LATS_NODE_R)
      integer kmsk  (lonr,lats_node_r),kmsk0(lonr,lats_node_r)
!
      kmsk=nint(sfc_fld%slmsk)
      kmsk0=0
!
      IF(FHOUR.GT.ZHOUR) THEN
        RTIME=1./(3600.*(FHOUR-ZHOUR))
      ELSE
        RTIME=0.
      ENDIF
!
!..........................................................
!
!      ngrid2d=1
!      CALL uninterprez(1,kmsk,glolal,sfc_fld%slmsk,
!     &       global_lats_r,lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
!      slmskloc=glolal
!      slmskful=buff_mult_pieceg(1:lonr,1:lats_node_r,ngrid2d)
!
!..........................................................
!
      ngrid2d=0
      if ( g2d_fld%du%nfld > 0 ) then
        do  k = 1, g2d_fld%du%nfld
          glolal=1.e6*g2d_fld%du%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo

      endif
!
!..........................................................
!
      if ( g2d_fld%su%nfld > 0 ) then
        do  k = 1, g2d_fld%su%nfld
          glolal=1.e6*g2d_fld%su%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif
!
!..........................................................
!
      if ( g2d_fld%ss%nfld > 0 ) then
        do  k = 1, g2d_fld%ss%nfld
          glolal=1.e6*g2d_fld%ss%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif
!
!..........................................................
!
      if ( g2d_fld%oc%nfld > 0 ) then
        do  k = 1, g2d_fld%oc%nfld
          glolal=1.e6*g2d_fld%oc%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif
!
!..........................................................
!
      if ( g2d_fld%bc%nfld > 0 ) then
        do  k = 1, g2d_fld%bc%nfld
          glolal=1.e6*g2d_fld%bc%diag(k)%flds
          ngrid2d=ngrid2d+1
          CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,
     &                   lonsperlar,buff_mult_pieceg(1,1,ngrid2d))
        enddo
      endif
!!

      if(me.eq.ioproc)
     &   PRINT *,'(wrtaer) GRIB AER FILE WRITTEN ',FHOUR,IDATE,noaer
!!
      RETURN
      END
!!****

       subroutine flx_only_move(ioproc)
!
!***********************************************************************
!
      use resol_def, ONLY: ngrids_flx, ngrids_sfcc, lonr,latr
      use mod_state, ONLY: buff_mult_pieces, buff_mult_piecef,
     &                     ivar_global_a, ivar_global
      use layout1,   ONLY: me, nodes, ipt_lats_node_r, lats_node_r,
     &                     lats_node_r_max, nodes_comp
      use mpi_def,   ONLY: mpi_r_io, stat, mpi_comm_null, info, 
     &                     mc_comp, mpi_integer, mpi_comm_all, liope
      implicit none
!
      integer ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
!      integer lats_nodes_r(nodes),ipt,maxfld,ioproc,nproct
      integer ioproc
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer illen,ubound,nd1
       integer icount
       data icount/0/
       integer maxlats_comp
       save maxlats_comp,icount
       integer kllen
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         if(icount .eq. 0) then
          if(.not.allocated(ivar_global)) allocate(ivar_global(10))
          if(.not.allocated(ivar_global_a)) 
     &       allocate(ivar_global_a(10,nodes))
         ivar_global(1)=ipt_lats_node_r
         ivar_global(2)= lats_node_r
         ivar_global(3)=lats_node_r_max
         call mpi_gather(ivar_global,10,MPI_INTEGER,
     1 ivar_global_a,10,MPI_INTEGER,ioproc,MPI_COMM_ALL,ierr)
         icount=icount+1
         endif
!!
       if(allocated(buff_mult_pieces)) then
          deallocate(buff_mult_pieces)
       else
          maxlats_comp=lats_node_r_max
          if(me .eq. ioproc) then
            maxlats_comp=ivar_global_a(3,1)
           endif
       endif
       if(me .eq. ioproc) then
!gwv watch this!!
          allocate
     1  (buff_mult_pieces(lonr*latr*ngrids_flx))
         buff_mult_pieces=0.
       endif
!
!  big send
       IF (me.ne.ioproc) THEN
!
!         Sending the data
         msgtag=me
         illen=lats_node_r
         kllen=illen*lonr*ngrids_flx
! send the local grid domain
         CALL mpi_send
     &(buff_mult_piecef,kllen,MPI_R_IO,ioproc,
     &                  msgtag,mc_comp,info)
!jw     &                  msgtag,MPI_COMM_ALL,info)
      ELSE
         if( MC_COMP .ne. MPI_COMM_NULL) then
! iotask is also a compute task.  send is replaced with direct
!  array copy
 
         if(nodes_comp==1) then
           buff_mult_pieces(1:lonr*lats_node_r*ngrids_flx)=
     1   reshape(buff_mult_piecef(1:lonr,1:lats_node_r,1:ngrids_flx),
     1     (/lonr*lats_node_r*ngrids_flx/) )
         else

!  END COMPUTE TASKS PORTION OF LOGIC
!  receiving part of I/O task
 
!!
!!     for pes ioproc
        nd1=0
        DO proc=1,nodes_comp
         illen=ivar_global_a(2,proc)
         if (proc.ne.ioproc+1) then
           msgtag=proc-1
           kllen=illen*lonr*ngrids_flx
           CALL mpi_recv
     1       (buff_mult_pieces(nd1+1),kllen,MPI_R_IO,proc-1,
     &                msgtag,mc_comp,stat,info)
!     &                msgtag,MPI_COMM_ALL,stat,info)
         else
           buff_mult_pieces(nd1+1:nd1+lonr*illen*ngrids_flx)=
     1      reshape(buff_mult_piecef(1:lonr,1:illen,1:ngrids_flx),
     1       (/lonr*illen*ngrids_flx/) )
         endif
         nd1=nd1+illen*lonr*ngrids_flx
        enddo
       endif

      endif
!end ioproc
      ENDIF
!
      return
      end
 
      SUBROUTINE flx_wrt(IOPROC,cfile,ZHOUR,FHOUR,idate
     &,                  global_lats_r,lonsperlar)
!!
      use module_nemsio, only: nemsio_open,nemsio_writerec,nemsio_close
     &  ,nemsio_gfile, nemsio_init,nemsio_finalize
      use resol_def,    ONLY: lonr, latr, levs,ngrids_flx,
     & ncld,ntrac,ntcw,ntoz,lsoil, ivssfc,thermodyn_id,sfcpress_id
      use layout1,      ONLY: me,idrt
      USE machine,      ONLY: kind_io8, kind_io4
!
      use gfs_physics_output, only : PHY_INT_STATE_ISCALAR,
     &    PHY_INT_STATE_RSCALAR,
     &    PHY_INT_STATE_1D_I,PHY_INT_STATE_1D_R,
     &    PHY_INT_STATE_2D_R_FLX
      implicit none
!!
      integer nw,IOPROC
      character*16 cfile,NAME2D
      real(kind=kind_io8) zhour,fhour
      integer idate(4),k,il, ngridss
!
      integer i,j,ndim3,N2DR,INDX,idate7(7),kount,nrec
      integer nfhour,nfminute,nfsecondn,nfsecondd
      logical  :: outtest
      integer ::nmetavari,nmetavarr,nmetavarl,nmetaaryi,nmetaaryr
      character(16),allocatable :: recname(:),reclevtyp(:)
      integer,allocatable :: reclev(:),itr(:)
      character(16),allocatable :: variname(:),varrname(:),
     &    aryiname(:),aryrname(:)
      integer,allocatable :: varival(:),aryilen(:),
     &    aryival(:,:)
      real(kind=kind_io4),allocatable    :: varrval(:)
      real(kind=kind_io4),allocatable    :: buff_mult(:,:,:),tmp(:)
      type(nemsio_gfile) gfileout
!

!!
      CHARACTER*8 labfix(4)
      real(kind=kind_io4) yhour
      integer,save:: version
      data version/200501/
      INTEGER              GLOBAL_LATS_R(latr), lonsperlar(latr)
!
      integer iret
      logical first
      save first
      save  recname, reclevtyp, reclev
      save nrec,nmetavari,nmetavarr,nmetaaryi,nmetaaryr,
     &     variname,varrname,aryiname,
     &     varival,varrval,aryilen,aryival
!jw     &     variname,varrname,aryiname,aryrname,
!jw     &     varival,aryilen,aryrlen,aryival,aryrval,varrval
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build surface fields in to buff_mult
!
      if (me.eq.ioproc) then
!
      print *,' begin of flx_wrt '

        allocate(buff_mult(lonr,latr,ngrids_flx))
        buff_mult=0.
        do ngridss=1,ngrids_flx
          print *,' inside flx_wrt calling unsp ngridss=',ngridss
          call unsplit2z(ioproc,ngridss,ngrids_flx,
     &      buff_mult(1,1,ngridss),global_lats_r)
        enddo
!    Building surface field is done
!
        if (first) then
!write out nemsio sfc file:
          nrec=ngrids_flx
          kount=size(PHY_INT_STATE_ISCALAR,2)
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_FLX')
     &        nmetavari=nmetavari+1
          enddo
          allocate(variname(nmetavari),varival(nmetavari))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_PHY' .or.
     &      trim(PHY_INT_STATE_ISCALAR(2,i)).eq.'OGFS_FLX' )then
            j=j+1
            variname(j)=trim(PHY_INT_STATE_ISCALAR(1,i))
            if(trim(variname(j))=='latr') varival(j)=latr
            if(trim(variname(j))=='lonr') varival(j)=lonr
            if(trim(variname(j))=='levs') varival(j)=levs
            if(trim(variname(j))=='ntoz') varival(j)=ntoz
            if(trim(variname(j))=='ntcw') varival(j)=ntcw
            if(trim(variname(j))=='ncld') varival(j)=ncld
            if(trim(variname(j))=='ntrac') varival(j)=ntrac
            if(trim(variname(j))=='thermodyn_id')varival(j)=thermodyn_id
            if(trim(variname(j))=='sfcpress_id') varival(j)=sfcpress_id
            if(trim(variname(j))=='lsoil') varival(j)=lsoil
            if(trim(variname(j))=='idrt') varival(j)=idrt
           endif
          enddo
!!for real var::
          nmetavarr=0
          do i=1,kount
           if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_FLX')
     &     nmetavarr=nmetavarr+1
          enddo
          if(nmetavarr>0) then
            allocate(varrname(nmetavarr),varrval(nmetavarr))
            j=0
            do i=1,kount
             if(trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_PHY'
     &       .or.trim(PHY_INT_STATE_RSCALAR(2,i)).eq.'OGFS_FLX')then
               j=j+1
               varrname(j)=trim(PHY_INT_STATE_RSCALAR(1,i))
               if(trim(varrname(j))=='fhour') varrval(j)=fhour
               if(trim(varrname(j))=='zhour') varrval(j)=zhour
             endif
            enddo
          endif
!!for 1D ary::
          nmetaaryi=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_FLX')
     &     nmetaaryi=nmetaaryi+1
          enddo
          allocate(aryiname(nmetaaryi),aryilen(nmetaaryi))
          j=0
          do i=1,kount
           if(trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_PHY'
     &     .or.trim(PHY_INT_STATE_1D_I(2,i)).eq.'OGFS_FLX')then
             j=j+1
             aryiname(j)=trim(PHY_INT_STATE_1D_I(1,i))
             if(trim(aryiname(j))=='IDATE') aryilen(j)=size(idate)
           endif
          enddo
          allocate(aryival(maxval(aryilen),nmetaaryi) )
          aryival(1:aryilen(1),1)=idate(:)
!!!for 1D real ary::
!          nmetaaryr=0
!          do i=1,kount
!           if(trim(PHY_INT_STATE_1D_R(2,i)).eq.'OGFS_PHY'
!     &     .or.trim(PHY_INT_STATE_1D_R(2,i)).eq.'OGFS_FLX')
!     &     nmetaaryr=nmetaaryr+1
!          enddo
!          allocate(aryrname(nmetaaryr),aryrlen(nmetaaryr))
!          do i=1,kount
!           if(trim(PHY_INT_STATE_1D_R(2,i)).eq.'OGFS_PHY')
!     &     .or.trim(PHY_INT_STATE_1D_R(2,i)).eq.'OGFS_FLX')then
!             aryrname(i)=trim(PHY_INT_STATE_1D_R(1,i))
!             if(i==1) aryrlen(i)=size(ak5)
!             if(i==2) aryrlen(i)=size(bk5)
!             if(i==3) aryrlen(i)=size(ck5)
!           endif
!          enddo
!          allocate(aryrval(maxval(aryrlen),nmetaaryr)
!          aryrval(1:aryrlen(1),1)=ak5(:)
!          aryrval(1:aryrlen(2),2)=bk5(:)
!          aryrval(1:aryrlen(3),2)=ck5(:)
!
!!for record name, levtyp and lev
          allocate (recname(nrec),reclevtyp(nrec),reclev(nrec))
          allocate (itr(nrec))
          N2DR=0
          itr=-99
          do i=1,kount
           if(trim(PHY_INT_STATE_2D_R_FLX(2,i)).eq.'OGFS_FLX')then
            N2DR=N2DR+1
            NAME2D=trim(PHY_INT_STATE_2D_R_FLX(1,i))
            INDX=INDEX(NAME2D,"_")
            if(indx>0) then
              recname(N2DR)=NAME2D(1:INDX-1)
            else
              recname(N2DR)=NAME2D
            endif
!
            reclevtyp(N2DR)=trim(trim(PHY_INT_STATE_2D_R_FLX(3,i)))
            reclev(N2DR)=1
!
!check time average
           if(INDEX(NAME2D,"_ave") >0) then
               itr(N2DR)=3
            elseif(INDEX(NAME2D,"_acc") >0) then
               itr(N2DR)=4
            elseif(INDEX(NAME2D,"_win") >0) then
               itr(N2DR)=2
            endif

           endif
          enddo
!
!end first
          first=.false.
         endif
!
        idate7=0
        idate7(1)=idate(4)
        idate7(2)=idate(2)
        idate7(3)=idate(3)
        idate7(4)=idate(1)
        idate7(7)=100           !denominator for second
!
        nfhour=int(fhour)     
        nfminute=int((fhour-nfhour)*60)
        nfsecondn=int(((fhour-nfhour)*3600-nfminute*60)*100)
        nfsecondd=100
!
        call nemsio_init()
!
        call nemsio_open(gfileout,trim(cfile),'write',
     &    iret = iret,
     &    modelname='GFS',gdatatype='grib',
     &    idate=idate7,nrec=nrec,
     &    dimx=lonr,dimy=latr,dimz=levs,ncldt=ncld,nmeta=5,
     &    nfhour=nfhour,nfminute=nfminute,nfsecondn=nfsecondn,
     &    nfsecondd=nfsecondd,
     &    extrameta=.true.,nmetavari=nmetavari,
     &    nmetavarr=nmetavarr,
     &    nmetaaryi=nmetaaryi,
     &    variname=variname,varival=varival,varrname=varrname,
     &    varrval=varrval,
     &    aryiname=aryiname,aryilen=aryilen,aryival=aryival,
     &    ntrac=ntrac,nsoil=lsoil,idrt=idrt,
     &    recname=recname,reclevtyp=reclevtyp,reclev=reclev)
!
        allocate(tmp(lonr*latr))
        yhour=zhour
        do i=1,nrec
          tmp(:)=reshape(buff_mult(:,:,i),(/lonr*latr/) )
          if(itr(i)==-99) then
            call nemsio_writerec(gfileout,i,tmp,iret=iret)
          else
            call nemsio_writerec(gfileout,i,tmp,iret=iret,itr=itr(i),
     &        zhour=yhour)
          endif
        enddo
        deallocate(tmp)
        deallocate(buff_mult)
!
        call nemsio_close(gfileout)
!end write pe
        call nemsio_finalize()
      endif
!
      print *,' end of flx_wrt '
      return
      end
!
      INTEGER FUNCTION nfill(C)
      implicit none
      integer j
      CHARACTER*(*) C
      NFILL=LEN(C)
      DO J=1,NFILL
        IF(C(J:J).EQ.' ') THEN
          NFILL=J-1
          RETURN
        ENDIF
      ENDDO
      RETURN
      END
 
 
