      module mod_state
!
!    New module to supply domain information to the GFS output routines
!    called by wrtout.
!
      use machine
      use resol_def
      implicit none
!
!
      real(kind=kind_io4), allocatable :: buff_mult_pieceg(:,:,:),
     1                                    buff_mult_piecesg(:,:,:,:)
!
      real(kind=kind_io4), allocatable :: buff_mult_piece(:,:,:),
     1                                    buff_mult_pieces(:,:,:,:)
      real(kind=kind_io4), allocatable :: buff_mult_piecef(:,:,:),
     1                                    buff_mult_piecesf(:,:,:,:)
      real(kind=kind_io4), allocatable :: buff_mult_piecea(:,:,:),
     1                                    buff_mult_piecesa(:,:,:,:)
      integer , allocatable :: ivar_global(:),ivar_global_a(:,:)
     &,                        ivarg_global(:),ivarg_global_a(:,:)
      integer , allocatable :: maskss(:,:,:)
!
      integer ngrid ,ngrida,ngridg
      save ngrid,ngrida,buff_mult_piece,buff_mult_pieces,ivar_global
     &,    ngridg,buff_mult_pieceg,buff_mult_piecesg,ivarg_global
      save maskss
      end module mod_state
      subroutine wrtout_physics(phour,fhour,zhour,idate,
     &                  sl,si,
     &                  sfc_fld, flx_fld,
     &                  fluxr,
     &                  lats_nodes_r,global_lats_r,lonsperlar,nblck,
     &                  colat1,cfhour1,pl_coeff)
!!
      use resol_def
      use layout1
      use namelist_physics_def
      use mpi_def
      use gfs_physics_sfc_flx_mod
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
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
      integer nosfc,noflx,nfill
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
!     real(kind=kind_rad) zsg(lonr,lats_node_r)
!     real(kind=kind_rad) psg(lonr,lats_node_r)
!     real(kind=kind_rad) uug(lonr,lats_node_r,levs)
!     real(kind=kind_rad) vvg(lonr,lats_node_r,levs)
!     real(kind=kind_rad) teg(lonr,lats_node_r,levs)
!     real(kind=kind_rad) rqg(lonr,lats_node_r,levh)
!     real(kind=kind_rad) dpg(lonr,lats_node_r,levs)
!!
      real secphy,secswr,seclwr
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
      if(liope) ioproc=nodes_comp
       
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
      CFHOUR = CFHOUR(1:nfill(CFHOUR)) // ens_nam(1:nfill(ens_nam))
      print *,' in wrtout_physics cfhour=',cfhour,' ens_nam=',ens_nam
!jfe
      nosfc=62
      noflx=63
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
!
      print *,'---- start sfc collection section -----'
      t3=rtc()
      if(mc_comp .ne. MPI_COMM_NULL) then
        print *,' wrtout_physics call sfc_collect '
        CALL sfc_collect(sfc_fld,global_lats_r,lonsperlar)
        IF(LIOPE) then
!
! collect flux grids as was done with sfc grids above.
! but only if liope is true.  If liope is false,
! the fluxes are handled by the original wrtsfc
! predating the I/O task updates.
!
            print *,' wrtout_physics call wrtflx_a '
            call   wrtflx_a
     &             (IOPROC,noflx,ZHOUR,FHOUR,IDATE,colat1,SECSWR,SECLWR,
     &              sfc_fld, flx_fld, fluxr, global_lats_r,lonsperlar)
        endif             ! liope
      endif                 ! comp node
      t4=rtc()
      td=t4-t3
!
!  done with state build
!  NOW STATE IS ASSEMBLED ON EACH NODE.  GET EVERYTHING OFF THE COMPUTE
!  NODES (currently done with a send to the I/O task_
!  send state to I/O task.  All tasks
!
      print *,'---- start flx.f section -----'
        IF(LIOPE) then                    ! move all grids (fluxes and sfc)
          print *,' wrtout_physics call grids_move '
          call GRIDS_MOVE(ioproc )
        ELSE             ! move sfc grids only,  handle fluxes in original wrtsfc
          print *,' wrtout_physics call sfc_only_move '
          call SFC_ONLY_MOVE(ioproc)
          if(me .eq. ioproc) then
            call BAOPENWT(NOFLX,'FLX.F'//CFHOUR,iostat)
          endif
          print *,' wrtout_physics call wrtsfc '
          call  WRTSFC
     &          (IOPROC,noflx,ZHOUR,FHOUR,IDATE,colat1,SECSWR,SECLWR,
     &           sfc_fld, flx_fld, fluxr, global_lats_r,lonsperlar)
        ENDIF          !  LIOPE
!
        t4=rtc()
        te=t4-t3
!
      print *,'---- start diag3d.f section -----'
        IF (LDIAG3D) THEN
          print *,' wrtout_physics ldiag3d on so wrt3d '
          no3d=64
          if(icolor.eq.2.and.me.eq.IOPROC)
     &    call BAOPENWT(NO3D,'D3D.F'//CFHOUR,iostat)
          if (hybrid .or. gen_coord_hybrid) then
!     print *,' pl_coeff bef call wrt3d_hyb=',pl_coeff
            call WRT3D_hyb(IOPROC,no3d,nblck,ZHOUR,FHOUR,IDATE,colat1,
     .                     global_lats_r,lonsperlar,pl_coeff,
     &                     SECSWR,SECLWR,sfc_fld%slmsk,flx_fld%psurf)
          else
            call WRT3D(IOPROC,no3d,nblck,ZHOUR,FHOUR,IDATE,colat1,
     .                 global_lats_r,lonsperlar,pl_coeff,
     &                 SECSWR,SECLWR,sl,si,sfc_fld%slmsk,flx_fld%psurf)
          endif
        ENDIF
!
      print *,'---- start sfc.f section -----'
! ioproc only
      CFHOUR1 = CFHOUR          !for the ESMF Export State Creation
      ta=rtc()
      if (me .eq. ioproc) then
          if(liope) call BAOPENWT(NOFLX,'SFC.F'//CFHOUR,iostat)
!
!  Now write the surface file
!
          print *,' wrtout_physics call sfc_wrt '
          cosfc='SFC.F'//CFHOUR
          call sfc_wrt(ioproc,nosfc,cosfc,fhour,jdate
     &,                global_lats_r,lonsperlar)
          CLOSE(NOSFC)
      endif
!
      tc=rtc()
      if(me .eq. 0) t2=rtc()
!gwv  t2=rtc()
      t3=rtc()
      if(MC_COMP   .ne. MPI_COMM_NULL) then
        call mpi_barrier(mc_comp,info)
      endif
      t4=rtc()
      if(liope) then                     !  WRITE THE FLUXES
        if(me .eq. ioproc) then
          print *,' wrtout_physics call wrtflx_w '
          call  WRTFLX_w
     &      (IOPROC,noflx,ZHOUR,FHOUR,IDATE,colat1,SECSWR,SECLWR,
     &       sfc_fld%slmsk, global_lats_r,lonsperlar)
        endif
      endif
!                                        !  FLUX WRITE DONE
!
      if(me .eq. ioproc) then
          call baclose(noflx,iostat)
          print *,' iostat after baclose of noflx ',iostat,noflx
      endif
!
      if(me .eq. ioproc)  call wrtlog_physics(phour,fhour,idate)
      tb=rtc()
      tf=tb-ta
      t2=rtc()
      print 1011,tf
 1011 format(' WRTOUT_PHYSICS TIME ',f10.4)
      timesum=timesum+(t2-t1)
      print 1012,timesum,t2-t1,td,te,tf,t4-t3,tba,tbb,tbd
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
      SUBROUTINE wrt_restart_physics(
     &        sfc_fld,
     &        SI,SL,fhour,idate,
     &        igen,
     &        global_lats_r,lonsperlar,
     &        phy_f3d, phy_f2d, ngptc, nblck, ens_nam)
!!
      use resol_def
      use layout1
      use mpi_def
      use gfs_physics_sfc_flx_mod
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
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
      INTEGER              GLOBAL_LATS_R(LATR)
      INTEGER              lonsperlar(LATR)
      integer IOPROC, IPRINT
      integer needoro, iret, nfill
!
!!
      real runid,usrid
      integer n3,n4,nflop
      character*20 cflop
      integer nn
!!
      IPRINT = 0
!
      print *,' cflop=',cflop,'ens_nam=',ens_nam(1:nfill(ens_nam))
!
!     print *,' in rest write fhour=',fhour
      IF (icolor.eq.2) then
         IOPROC=nodes-1
      else
         IOPROC=nodes
      endif
!
      ixgr = 0
        if (num_p3d .eq. 4) then          ! Zhao microphysics
!         ixgr = 2
          ixgr = 4
        elseif (num_p3d .eq. 3) then      ! Ferrier microphysics
!         ixgr = 3
          ixgr = 5
        endif
!     xgf = ixgf
!
      IF (icolor.eq.2.and.me.eq.ioproc) print *,' closed ',n3
!
      IF (icolor.eq.2) then
         IOPROC=nodes-1
      else
         IOPROC=nodes
      endif

      ixgr = 0

      nflop=53
!     cflop='fort.53'
      IF (icolor.eq.2) then
         IOPROC=nodes-1
      else
         IOPROC=nodes
      endif
        CALL para_fixio_w(ioproc,sfc_fld, nflop,cflop,fhour,idate,
     &                    global_lats_r,lonsperlar)
!
      return
      end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE wrtlog_physics(phour,fhour,idate)
      use resol_def
      use namelist_physics_def
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
      CFHOUR = CFHOUR(1:nfill(CFHOUR)) // ens_nam(1:nfill(ens_nam))

      nolog=99
      OPEN(NOlog,FILE='LOG.F'//CFHOUR,FORM='FORMATTED')
      write(nolog,100)fhour,idate
100   format(' completed mrf fhour=',f10.3,2x,4(i4,2x))
      CLOSE(NOlog)

      RETURN
      END



      SUBROUTINE sfc_collect (sfc_fld,global_lats_r,lonsperlar)
!!
      use resol_def
      use mod_state
      use layout1
      use mpi_def
      use gfs_physics_sfc_flx_mod
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
      ngrid=1
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!!
      if(allocated(buff_mult_piece)) then
         continue
      else
         allocate(buff_mult_piece(lonr,ngrids_sfcc,lats_node_r))
         allocate(buff_mult_piecef(lonr,0:ngrids_flx,lats_node_r))
         allocate
     1 (buff_mult_piecea(lonr,1:ngrids_flx+ngrids_sfcc+1,lats_node_r))
      endif
!
      kmsk= nint(sfc_fld%slmsk)
      CALL uninterprez(1,kmsk,buffo,sfc_fld%tsea,
     &                global_lats_r,lonsperlar)
!
! ngrid=2 here
                                                                                                        
!
      DO k=1,LSOIL
        buffi(:,:) = sfc_fld%SMC(k,:,:)
        CALL uninterprez(1,kmsk,buffo,buffi,global_lats_r,lonsperlar)
      ENDDO
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SHELEG,
     &                 global_lats_r,lonsperlar)
!
      DO k=1,LSOIL
        buffi(:,:) = sfc_fld%STC(k,:,:)
        CALL uninterprez(1,kmsk,buffo,buffi,global_lats_r,lonsperlar)
      ENDDO
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%TG3,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ZORL,
     &                 global_lats_r,lonsperlar)
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
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALVSF,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALVWF,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALNSF,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ALNWF,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SLMSK,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%VFRAC,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%CANOPY,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%F10M,
     &                 global_lats_r,lonsperlar)
! T2M
      CALL uninterprez(1,kmsk,buffo,sfc_fld%T2M,
     &                 global_lats_r,lonsperlar)
! Q2M
      CALL uninterprez(1,kmsk,buffo,sfc_fld%Q2M,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%VTYPE,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%STYPE,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FACSF,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FACWF,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%UUSTAR,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FFMM,
     &                 global_lats_r,lonsperlar)
!
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FFHH,
     &                 global_lats_r,lonsperlar)
!
!c-- XW: FOR SEA-ICE Nov04
      CALL uninterprez(1,kmsk,buffo,sfc_fld%HICE,
     &                 global_lats_r,lonsperlar)
      CALL uninterprez(1,kmsk,buffo,sfc_fld%FICE,
     &                 global_lats_r,lonsperlar)
      CALL uninterprez(1,kmsk,buffo,sfc_fld%TISFC,
     &                 global_lats_r,lonsperlar)
!c-- XW: END SEA-ICE Nov04
!
!lu: the addition of 8 Noah-related records starts here ........................
!tprcp
      CALL uninterprez(1,kmsk,buffo,sfc_fld%TPRCP,
     &                 global_lats_r,lonsperlar)
!srflag
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SRFLAG,
     &                 global_lats_r,lonsperlar)
!snwdph
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SNWDPH,
     &                 global_lats_r,lonsperlar)
!slc
      DO k=1,LSOIL
        buffi(:,:) = sfc_fld%SLC(k,:,:)
        CALL uninterprez(1,kmsk,buffo,buffi,global_lats_r,lonsperlar)
!       buffo(:,:)=buff_mult_piece(:,k+3+lsoil,:)
      ENDDO
!shdmin
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SHDMIN,
     &                 global_lats_r,lonsperlar)
!shdmax
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SHDMAX,
     &                 global_lats_r,lonsperlar)
!slope
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SLOPE,
     &                 global_lats_r,lonsperlar)
!snoalb
      CALL uninterprez(1,kmsk,buffo,sfc_fld%SNOALB,
     &                 global_lats_r,lonsperlar)
!lu: the addition of 8 Noah records ends here .........................
!
! Oro
      CALL uninterprez(1,kmsk,buffo,sfc_fld%ORO,
     &                 global_lats_r,lonsperlar)
!
!     print *,' finished sfc_collect for  ngrid=',ngrid
  999 continue
      ngrid=1
      return
      end
       subroutine sfc_only_move(ioproc)
!
!***********************************************************************
!
      use resol_def
      use mod_state
      use layout1
      use mpi_def
      implicit none
!
      integer ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
!     integer lats_nodes_r(nodes),ipt,maxfld,ioproc,nproct
      integer ioproc
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer illen,ubound
      integer icount
      data icount/0/
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
     1       ivar_global_a,10,MPI_INTEGER,ioproc,MPI_COMM_ALL,ierr)
         icount=icount+1
      endif
!!
      if(allocated(buff_mult_pieces)) then
          continue
      else
          maxlats_comp=lats_node_r_max
          if(.not. liope .or. me .ne. ioproc) then
            continue
          else
!           maxlats_comp=ivar_global_a(3,ioproc)
            maxlats_comp=ivar_global_a(3,1)
          endif
          print *,' INDEX FOR MAXLAT SET ',ioproc
!gwv watch this!!
!         print *,' allocating ', lonr,ngrids_sfcc,maxlats_comp,nodes
          allocate
     1    (buff_mult_pieces(lonr,ngrids_sfcc,maxlats_comp,nodes))
!         print *,' allocated', lonr,ngrids_sfcc,maxlats_comp,nodes
          allocate
     1    (buff_mult_piecesf(lonr,0:ngrids_flx,maxlats_comp,nodes))
          allocate
     1    (buff_mult_piecesa(lonr,1:ngrids_flx+1+ngrids_sfcc,
     1     maxlats_comp,nodes))
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
       buff_mult_piece(:,1:ngrids_sfcc,:)=
     1 buff_mult_piecea(:,1:ngrids_sfcc,:)
!
      IF (ME .ne. ioproc) THEN    !   Sending the data
         msgtag=me
         illen=lats_node _r
         CALL mpi_send            !  send the local grid domain
     &(buff_mult_piece,illen*lonr*ngrids_sfcc,MPI_R_IO,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
      ELSE
        if( MC_COMP .ne. MPI_COMM_NULL) then
!
c iotask is also a compute task.  send is replaced with direct
c  array copy
!
                                                                                                        
           buff_mult_pieces(:,:,1:lats_node_r,ioproc+1)=
     1     buff_mult_piece(:,:,1:lats_node_r)
!                              END COMPUTE TASKS PORTION OF LOGIC
        endif
!
!  END COMPUTE TASKS PORTION OF LOGIC
!  receiving part of I/O task
!
!!
!!      for pes ioproc
        DO proc=1,nodes_comp
          if (proc.ne.ioproc+1) then
            msgtag=proc-1
            illen=ivar_global_a(2,proc)
!           print *,' pux target ',ubound(buff_mult_pieces)
            CALL mpi_recv(buff_mult_pieces(1,1,1,proc),
     1        illen*lonr*ngrids_sfcc
     1        ,MPI_R_IO,proc-1,
     &                msgtag,MPI_COMM_ALL,stat,info)
          endif
        enddo
        buff_mult_piecesa(:,1:ngrids_sfcc,:,:)=
     1 buff_mult_pieces(:,1:ngrids_sfcc,:,:)
      ENDIF
!!
      return
      end
      SUBROUTINE sfc_wrt(IOPROC,nw,cfile,xhour,idate
     &,                  global_lats_r,lonsperlar)
!!
      use sfcio_module
      use resol_def
      use mod_state
      use layout1
      use mpi_def
!     use mod_state , only : ngrids_sfcc
      implicit none
!!
      integer nw,IOPROC
      character*16 cfile
      real(kind=kind_io8) xhour
!!!   real(kind=kind_io4) buff4(lonr,latr,4)
      integer idate(4),k,il, ngridss
!     integer idate(4),k,il, ngrid, ngridss
!!
      CHARACTER*8 labfix(4)
      real(kind=kind_io4) yhour
      integer,save:: version
      data version/200501/
      INTEGER              GLOBAL_LATS_R(latr), lonsperlar(latr)
!
      type(sfcio_head) head
      type(sfcio_data) data
      integer iret
      logical first
      save head, first
      data first /.true./
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!    Build surface fields in to buff_mult
!
      print *,' begin of sfc_wrt '

      ngrid=1
      do ngridss=1,ngrids_sfcc
        print *,' inside sfc_wrt calling unsp ngridss=',ngridss
        call unsplit2z(ioproc,buff_mult(1,1,ngridss),global_lats_r)
      enddo
!    Building surface field is done
!
      if (me.eq.ioproc) then
!
        if (first) then
          head%clabsfc = CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)//
     &                   CHAR(0)//CHAR(0)//CHAR(0)//CHAR(0)
          head%latb    = latr
          head%lonb    = lonr
          head%ivs     = ivssfc
!         head%irealf  = 1
          head%lsoil   = lsoil
          call sfcio_alhead(head,iret)
          head%lpl     = lonsperlar(1:latr/2)
          if (lsoil .eq. 4) then
            head%zsoil   = (/-0.1,-0.4,-1.0,-2.0/)
          elseif (lsoil .eq. 2) then
            head%zsoil   = (/-0.1,-2.0/)
          endif
          first = .false.
        endif
        head%fhour   = xhour
        head%idate   = idate
!
        PRINT 99,nw,xhour,IDATE
99      FORMAT(1H ,'in fixio nw=',i7,2x,'HOUR=',f8.2,3x,'IDATE=',
     &  4(1X,I4))
!
        ngrid = 1

        data%tsea=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%smc=>buff_mult(:,:,ngrid:ngrid+lsoil-1)
        ngrid=ngrid+lsoil
        data%sheleg=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%stc=>buff_mult(:,:,ngrid:ngrid+lsoil-1)
        ngrid=ngrid+lsoil
        data%tg3=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%zorl=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%alvsf=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%alvwf=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%alnsf=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%alnwf=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%slmsk=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%vfrac=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%canopy=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%f10m=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%t2m=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%q2m=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%vtype=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%stype=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%facsf=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%facwf=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%uustar=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%ffmm=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%ffhh=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
!c-- XW: FOR SEA-ICE Nov04
        data%hice=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%fice=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%tisfc=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
!c-- XW: END SEA-ICE Nov04
!
!lu: the addition of 8 Noah-related records starts here ...............
        data%tprcp=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%srflag=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%snwdph=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%slc=>buff_mult(:,:,ngrid:ngrid+lsoil-1)
        ngrid=ngrid+lsoil
        data%shdmin=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%shdmax=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%slope=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
        data%snoalb=>buff_mult(:,:,ngrid)
        ngrid=ngrid+1
!lu: the addition of 8 Noah records ends here .........................
!        
   
        data%orog=>buff_mult(:,:,ngrid)      ! Orography
!
!       ngrid=ngrid+1
!
! Not needed for version 200501
!       data%cv=>buff_mult(:,:,ngrid)
!       data%cvb=>buff_mult(:,:,ngrid)
!       data%cvt=>buff_mult(:,:,ngrid)
!
        call sfcio_swohdc(nw,cfile,head,data,iret)
!
      endif
      print *,' end of sfc_wrt '
      return
      end
      SUBROUTINE wrtflx_a(IOPROC,noflx,ZHOUR,FHOUR,IDATE,colat1,
     &                  SECSWR,SECLWR, sfc_fld, flx_fld, fluxr,
     &                  global_lats_r,lonsperlar)
!!
      use resol_def
      use mod_state
      use layout1
      use namelist_physics_def
      use gfs_physics_sfc_flx_mod
      implicit none
!!
      TYPE(Sfc_Var_Data)        :: sfc_fld
      TYPE(Flx_Var_Data)        :: flx_fld
      INTEGER              GLOBAL_LATS_R(LATR)
      INTEGER              lonsperlar(LATR)
      integer   IOPROC
!!
      integer LEN,NFLD
      integer j,i,k,itop,ibot,k4,l,noflx
      PARAMETER(NFLD=18)
       integer ilpds,iyr,imo,ida,ihr,ifhr,ithr,lg,ierr
       real (kind=kind_io8) RTIMER(NFLD),rtime,rtimsw,rtimlw
       real (kind=kind_io8) colat1
       real (kind=kind_io8) cl1,secswr,zhour,fhour,seclwr
C

      real(kind=kind_io4) wrkga(lonr*latr),wrkgb(lonr*latr)
      real(kind=kind_io8) slmskful(lonr*latr)
      real(kind=kind_io8) slmskloc(LONR,LATS_NODE_R)
!
      INTEGER     IDATE(4), IDS(255),IENS(5)
!     INTEGER     IDATE(4)
      real (kind=kind_io8) SI(LEVP1)
!
!sela..................................................................
      real (kind=kind_io8)   rflux(lonr,LATS_NODE_R,27)
      real (kind=kind_io8)   glolal(lonr,LATS_NODE_R)
      real (kind=kind_io8)   buffo(lonr,LATS_NODE_R)
      real (kind=kind_io4)   buff1(lonr,latr)
      real (kind=kind_io4)   buff1l(lonr*latr)
!sela..................................................................
      real (kind=kind_io8)  FLUXR(nfxr,LONR,LATS_NODE_R)
!sela..................................................................
      integer kmsk(lonr,lats_node_r),kmsk0(lonr,lats_node_r)
      integer kmskcv(lonr,LATS_NODE_R)
!
      ngrid=ngrids_sfcc+1
!
!!
      kmsk=nint(sfc_fld%slmsk)
      kmsk0=0
      CALL uninterprez(1,kmsk,glolal,sfc_fld%slmsk,
     &                global_lats_r,lonsperlar)
      slmskloc=glolal
      slmskful=buff1l
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
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '01)Zonal compt of momentum flux (N/m**2) land and sea surface '
 
!..........................................................
      glolal=flx_fld%DVSFC*RTIME
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '02)Merid compt of momentum flux (N/m**2) land and sea surface '
!..........................................................
      glolal=flx_fld%DTSFC*RTIME
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '03)Sensible heat flux (W/m**2) land and sea surface           '
!..........................................................
      glolal=flx_fld%DQSFC*RTIME
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '04)Latent heat flux (W/m**2) land and sea surface             '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%tsea,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribsn ierr=',ierr,'  ',
!    x '05)Temperature (K) land and sea surface                       '
!..........................................................
      glolal(:,:) = sfc_fld%SMC(1,:,:)
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribsn ierr=',ierr,'  ',
!    x '06)Volumetric soil moist content (frac) layer 10cm and 0cm    '
!..........................................................
      glolal(:,:) = sfc_fld%SMC(2,:,:)
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!lu  x '07)Volumetric soil moist content (frac) layer 200cm and 10cm  '
!    + '07)Volumetric soil moist content (frac) layer 40cm and 10cm  '
!..........................................................
      glolal(:,:) = sfc_fld%STC(1,:,:)
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '08)Temp (K) layer betw two depth below land sfc 10cm and 0cm  '
!..........................................................
      glolal(:,:) = sfc_fld%STC(2,:,:)
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!lu  x '09)Temp (K) layer betw two depth below land sfc 200cm and 10cm'
!    + '09)Temp (K) layer betw two depth below land sfc 40cm and 10cm'
!..........................................................
      CALL uninterprez(2,kmsk,buffo,sfc_fld%sheleg,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '10)Water equiv of accum snow depth (kg/m**2) land sea surface '
c..........................................................
      glolal = flx_fld%DLWSFC*RTIME
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '11)Downward long wave radiation flux (W/m**2) land sea surface'
!..........................................................
      glolal = flx_fld%ULWSFC*RTIME
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '17)UV-B Downward solar flux (W/m**2) land sea surface'
      do j=1,LATS_NODE_R
        do i=1,lonr
          glolal(i,j) = rflux(i,j,22)*rtimsw
        enddo
      enddo
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '31)Precipitation rate (kg/m**2/s) land and sea surface        '
c...................................................................
      glolal = flx_fld%BENGSH*1.E3*RTIME
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '32)Convective precipitation rate (kg/m**2/s) land sea surface '
!...................................................................
      glolal = flx_fld%GFLUX*RTIME
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
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
      buffo=MOD(slmskloc,2._kind_io8)
      do j=1,lats_node_r
        do i=1,lonr
          buff_mult_piecea(i,ngrid,j)=buffo(i,j)
        end do
      end do
        ngrid=ngrid+1
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '34)Land-sea mask (1=land; 0=sea) (integer) land sea surface   '
!gwv   add something here
!
!c-- XW: FOR SEA-ICE Nov04
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%fice,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '35)Ice concentration (ice>0; no ice=0) (1/0) land sea surface '
!c-- XW: END SEA-ICE
!...................................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%u10m,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '36)u wind (m/s) height above ground                           '
!...................................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%v10m,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '37)v wind (m/s) height above ground                           '
!...................................................................
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%t2m,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '38)Temperature (K) height above ground                        '
!...................................................................
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%q2m,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '39)Specific humidity (kg/kg) height above ground              '
!...................................................................
      glolal = flx_fld%PSURF*1.E3
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '40)Pressure (Pa) land and sea surface                         '
!...................................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%tmpmax,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '41)Maximum temperature (K) height above ground                '
!...................................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%tmpmin,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '42)Minimum temperature (K) height above ground                '
!...................................................................
      glolal = flx_fld%RUNOFF * 1.E3
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '43)Runoff (kg/m**2) land and sea surface                      '
!...................................................................
      glolal = flx_fld%EP * RTIME
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '44)Potential evaporation rate (w/m**/) land and sea surface   '
!...................................................................
      glolal = flx_fld%CLDWRK * RTIME
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '45)Cloud work function (J/Kg) total atmospheric column        '
!...................................................................
      glolal = flx_fld%DUGWD*RTIME
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '46)Zonal gravity wave stress (N/m**2) land and sea surface    '
!...................................................................
      glolal = flx_fld%DVGWD*RTIME
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '47)Meridional gravity wave stress (N/m**2) land sea surface   '
!...................................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%hpbl,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '48)Boundary layer height '
!...................................................................
!hmhj CALL uninterprez(2,kmsk0,buffo,flx_fld%pwat,
!hmhj&                 global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '50)Albedo (percent) land and sea surface                      '
!
       do j=1,LATS_NODE_R
        do i=1,lonr
         glolal(i,j) = rflux(i,j,26)*100.*rtimsw
        enddo
       enddo
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmskcv,buffo,glolal,global_lats_r,lonsperlar)
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
      CALL uninterprez(2,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '55)Total cloud cover (percent) boundary layer cloud layer     '
!c-- XW: FOR SEA-ICE Nov04
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%hice,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '56)Sea ice thickness (m) category 1'
!c-- XW: END SEA-ICE
!.................................................
!lu: add smc(3:4), stc(3:4), slc(1:4), snwdph, canopy
!lu: addition of 10 records starts here -------------------------------
      if(lsoil.gt.2)then
        glolal(:,:) = sfc_fld%SMC(3,:,:)
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '57)Volumetric soil moist content (frac) layer 100cm and 40cm '
!..........................................................
        glolal(:,:) = sfc_fld%SMC(4,:,:)
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '58)Volumetric soil moist content (frac) layer 200cm and 100cm '
!..........................................................
        glolal(:,:) = sfc_fld%STC(3,:,:)
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '59)Temp (K) layer betw two depth below land sfc 100cm and 40cm'
!..........................................................
        glolal(:,:) = sfc_fld%STC(4,:,:)
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '60)Temp (K) layer betw two depth below land sfc 200cm and 100cm'
      endif
!..........................................................
      glolal(:,:) = sfc_fld%SLC(1,:,:)
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '61)Liquid soil moist content (frac) layer 10cm and 0cm  '
!..........................................................
      glolal(:,:) = sfc_fld%SLC(2,:,:)
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '62)Liquid soil moist content (frac) layer 40cm and 10cm '
!..........................................................
      if(lsoil.gt.2)then
        glolal(:,:) = sfc_fld%SLC(3,:,:)
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '63)Liquid soil moist content (frac) layer 100cm and 40cm'
!..........................................................
        glolal(:,:) = sfc_fld%SLC(4,:,:)
        CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    &   '64)Liquid soil moist content (frac) layer 200cm and 100cm'
      endif
!..........................................................
      glolal = sfc_fld%SNWDPH / 1.E3       !! convert from mm to m
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '65)Snow depth (m) land surface                  '
c..........................................................
!     LBM=slmskful.EQ.1._kind_io8
      CALL uninterprez(2,kmsk,buffo,sfc_fld%canopy,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '66)Canopy water content (kg/m^2) land surface      '
!lu: addition of 10 records ends here -------------------------------
!
!wei: addition of 30 records starts here -------------------------------
      glolal = sfc_fld%ZORL / 1.E2       !! convert from cm to m
      CALL uninterprez(1,kmsk0,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '67)Surface roughness (m)       '
!..........................................................
      glolal = sfc_fld%vfrac*100.
      CALL uninterprez(1,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '68)Vegetation fraction (fractional) land surface      '
!..........................................................
      CALL uninterprez(1,kmsk,glolal,sfc_fld%vtype,
     &                 global_lats_r,lonsperlar)
      buffo=MOD(glolal,2._kind_io8)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '69)Vegetation type land surface      '
!..........................................................
      CALL uninterprez(1,kmsk,glolal,sfc_fld%stype,
     &                 global_lats_r,lonsperlar)
      buffo=MOD(glolal,2._kind_io8)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '70)Soil type land surface      '
!..........................................................
      CALL uninterprez(1,kmsk,glolal,sfc_fld%slope,
     &                 global_lats_r,lonsperlar)
      buffo=MOD(glolal,2._kind_io8)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '71)Slope type land surface      '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,sfc_fld%uustar,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '72)Frictional velocity (m/s)     '
!..........................................................
      CALL uninterprez(1,kmsk,buffo,sfc_fld%oro,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '73)Surface height (m)       '
!..........................................................
      CALL uninterprez(1,kmsk,buffo,sfc_fld%srflag,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '74)Freezing precip flag land surface      '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%chh,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '75)Exchange coefficient CH(m/s)       '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%cmm,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '76)Exchange coefficient CM(m/s)       '
!..........................................................
      CALL uninterprez(2,kmsk,buffo,flx_fld%EPI,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '77)Potential evaporation rate (w/m**2) land and sea surface   '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%DLWSFCI,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '78)Downward long wave radiation flux (W/m**2) '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%ULWSFCI,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '79)Upward long wave radiation flux (W/m**2)  '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%USWSFCI,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '80)Upward short wave radiation flux (W/m**2)  '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%DSWSFCI,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '81)Downward short wave radiation flux (W/m**2)   '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%DTSFCI,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '82)Sensible heat flux (W/m**2) land and sea surface       '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%DQSFCI,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '83)Latent heat flux (W/m**2) land and sea surface         '
!..........................................................
      CALL uninterprez(2,kmsk,buffo,flx_fld%GFLUXI,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '84)Ground heat flux (W/m**2) land and sea surface         '
!..........................................................
      glolal = flx_fld%SRUNOFF * 1.E3
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '85)Surface runoff (kg/m^2) land surface      '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%t1,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '86)Lowest model level Temp (K)      '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%q1,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '87)Lowest model specific humidity (kg/kg)    '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%u1,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '88)Lowest model u wind (m/s)      '
!..........................................................
      CALL uninterprez(2,kmsk0,buffo,flx_fld%v1,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '89)Lowest model v wind (m/s)       '
!..........................................................
      CALL uninterprez(2,kmsk,buffo,flx_fld%zlvl,
     &                 global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    x '90)Lowest model level height (m) land surface      '
!..........................................................
      glolal = flx_fld%EVBSA*RTIME
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '91)Direct evaporation from bare soil(W/m^2) land surface      '
!..........................................................
      glolal = flx_fld%EVCWA*RTIME
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '92)Canopy water evaporation(W/m^2) land surface      '
!..........................................................
      glolal = flx_fld%TRANSA*RTIME
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '93)Transpiration (W/m^2) land surface      '
!..........................................................
      glolal = flx_fld%SBSNOA*RTIME
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '94)Snow Sublimation (W/m^2) land surface      '
!..........................................................
      glolal = flx_fld%SNOWCA*RTIME*100.
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '95)Snow Cover (fraction) land surface      '
!..........................................................
      glolal = flx_fld%soilm*1.E3       !! convert from m to (mm)kg/m^2
      CALL uninterprez(2,kmsk,buffo,glolal,global_lats_r,lonsperlar)
!     if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '96)Total column soil moisture (Kg/m^2) land surface      '



Cwei: addition of 30 records ends here -------------------------------
      if(me.eq.ioproc)
     &   PRINT *,'(wrtflx_a) GRIB FLUX FILE WRITTEN ',FHOUR,IDATE,noflx
!!
      RETURN
      END
      SUBROUTINE wrtflx_w(IOPROC,noflx,ZHOUR,FHOUR,IDATE,colat1,SECSWR,
     &                  SECLWR,slmsk, global_lats_r,lonsperlar)
!
      use resol_def
      use mod_state
      use layout1
      use namelist_physics_def
      implicit none
!!
      INTEGER              GLOBAL_LATS_R(LATR)
      INTEGER              lonsperlar(LATR)
      integer   IOPROC
!!
      integer   IPRS,ITEMP,IZNLW,IMERW,ISPHUM,IPWAT,
     $          IPCPR,ISNOWD,ICLDF,ICCLDF,
     $          ISLMSK,IZORL,IALBDO,ISOILM,ICEMSK,
     $          ILHFLX,ISHFLX,IZWS,IMWS,IGHFLX,
     $          IUSWFC,IDSWFC,IULWFC,IDLWFC,
     $          INSWFC,INLWFC,
     $          IDSWVB,IDSWVD,IDSWNB,IDSWND,
     $          ITMX,ITMN,IRNOF,IEP,
     &          ICLDWK,IZGW,IMGW,IHPBL,
     $          IDSWF,IDLWF,IUSWF,IULWF,ICPCPR,
     $          ISFC,ITOA,IELEV,
     $          ISGLEV,IDBLS,I2DBLS,ICOLMN,
     $          IBLBL,IBLTL,IBLLYR,
     $          ILCBL,ILCTL,ILCLYR,
     $          IMCBL,IMCTL,IMCLYR,
     $          IHCBL,IHCTL,IHCLYR,
     $          ICVBL,ICVTL,ICVLYR,
     $          INST,IWIN,IAVG,IACC,
     $          IFHOUR,IFDAY,
!    $          LEN,NFLD,
     $          NFLD,
     $          IUVBF,IUVBFC,
     $   j,i,k,itop,ibot,k4,l,noflx
     &,  isik, islc, isnod, icnp
     &,  iveg, ivtp, istp, islo,iust,ihgt,irst,ichh
     &,  icmm,isrf,ievbs,ievcw,itran,isbs,isnc,istc
!     PARAMETER(NFLD=16)
      PARAMETER(NFLD=18)
       integer ilpds,iyr,imo,ida,ihr,ifhr,ithr,lg,ierr
       real (kind=kind_io8) RTIMER(NFLD),rtime,rtimsw,rtimlw
       real (kind=kind_io8) colat1
       real (kind=kind_io8) cl1,secswr,zhour,fhour,seclwr
!
      PARAMETER(IPRS=1,ITEMP=11,IZNLW=33,IMERW=34,ISPHUM=51,IPWAT=54,
     $          IPCPR=59,ISNOWD=65,ICLDF=71,ICCLDF=72,
     $          ISLMSK=81,IZORL=83,IALBDO=84,ISOILM=144,ICEMSK=91,
     $          ISIK=92,                                ! FOR SEA-ICE - XW Nov04
     $          ILHFLX=121,ISHFLX=122,IZWS=124,IMWS=125,IGHFLX=155,
     $          IUSWFC=160,IDSWFC=161,IULWFC=162,IDLWFC=163,
     $          INSWFC=164,INLWFC=165,
     $          IDSWVB=166,IDSWVD=167,IDSWNB=168,IDSWND=169,
     $          ITMX=15,ITMN=16,IRNOF=90,IEP=145,
     &          ICLDWK=146,IZGW=147,IMGW=148,IHPBL=221,
     $          IDSWF=204,IDLWF=205,IUSWF=211,IULWF=212,ICPCPR=214,
     &          IUVBF=200,IUVBFC=201)
      PARAMETER(ISFC=1,ITOA=8,IELEV=105,
     $          ISGLEV=109,IDBLS=111,I2DBLS=112,ICOLMN=200,
!Cwei    $          ISGLEV=107,IDBLS=111,I2DBLS=112,ICOLMN=200,
     $          IBLBL=209,IBLTL=210,IBLLYR=211,
     $          ILCBL=212,ILCTL=213,ILCLYR=214,
     $          IMCBL=222,IMCTL=223,IMCLYR=224,
     $          IHCBL=232,IHCTL=233,IHCLYR=234,
     $          ICVBL=242,ICVTL=243,ICVLYR=244)

!Clu [+1L]: define parameter index, using Table 130
      PARAMETER(ISLC=160,ISNOD=66)
!Cwei
      PARAMETER(ISLO=222,ISBS=198,ISNC=238,ICMM=179)
!Clu [+1L]: define parameter index, using Table 2
      PARAMETER(ICNP=223)
!Cwei
      PARAMETER(IVEG=87,IVTP=225,ISTP=224,IUST=253,IHGT=7,
     $          IRST=140,ICHH=208,ISRF=235,IEVBS=199,
     $          IEVCW=200,ITRAN=210,ISTC=86)

      PARAMETER(INST=10,IWIN=2,IAVG=3,IACC=4)
      PARAMETER(IFHOUR=1,IFDAY=2)
!     PARAMETER(LEN=lonr*latr)
      real(kind=kind_io4) wrkga(lonr*latr),wrkgb(lonr*latr)
      real(kind=kind_io8) slmskful(lonr*latr)
      real(kind=kind_io8) slmskloc(LONR,LATS_NODE_R)
c
      LOGICAL(1) LBM(lonr*latr)
      CHARACTER G(200+lonr*latr*(16+1)/8)
      INTEGER   IPUR(NFLD),ITLR(NFLD)
      DATA      IPUR/IULWF , IUSWF , IUSWF , IDSWF ,  ICLDF,   IPRS,
     $                 IPRS, ITEMP ,  ICLDF,   IPRS,   IPRS, ITEMP ,
     $                ICLDF,   IPRS,   IPRS, ITEMP ,  IUVBF, IUVBFC /
!    $                ICLDF,   IPRS,   IPRS, ITEMP /
      DATA      ITLR/ITOA  , ITOA  , ISFC  , ISFC  , IHCLYR, IHCTL ,
     $               IHCBL , IHCTL , IMCLYR, IMCTL , IMCBL , IMCTL ,
     $               ILCLYR, ILCTL , ILCBL , ILCTL , ISFC  , ISFC /
!    $               ILCLYR, ILCTL , ILCBL , ILCTL /
!     INTEGER     IDATE(4), IDS(255)
      INTEGER     IDATE(4), IDS(255),IENS(5),iensi,ienst,icen,icen2
      real (kind=kind_io8) SI(LEVP1)
C
csela..................................................................
!     real (kind=kind_io8)   rflux(lonr,LATS_NODE_R,27)
!     real (kind=kind_io8)   glolal(lonr,LATS_NODE_R)
!     real (kind=kind_io8)   buffo(lonr,LATS_NODE_R)
!     real (kind=kind_io4)   buff1(lonr,latr)
      real (kind=kind_io4)   buff1l(lonr*latr)
csela..................................................................
!     real (kind=kind_io8)  FLUXR(nfxr,LONR,LATS_NODE_R)
      REAL (KIND=KIND_IO8) SLMSK (LONR,LATS_NODE_R)
csela..................................................................
      integer kmsk(lonr,lats_node_r),kmsk0(lonr,lats_node_r)
      integer kmskcv(lonr,LATS_NODE_R),il
        ngrid=0
        ngrid=0+ngrids_sfcc+1
cjfe
      IDS=0
      G=' '
cjfe
!!
      kmsk=nint(slmsk)
      kmsk0=0
      call unsplit2z(ioproc,buff1l,global_lats_r)
      slmskful=buff1l
!
!     do k=1,nfxr
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!        rflux(i,j,k)=fluxr(k,i,j)
!       enddo
!      enddo
!     enddo
!!
      CALL IDSDEF(1,IDS)
! UV-B scaling factor, if set up already, comment the next 2 lines out
      ids(IUVBF)  = 2
      ids(IUVBFC) = 2
! Ice conentration and thickness scaling factor
      ids(icemsk) = 3      ! ICE CONCENTRATION ()
      ids(isik)   = 2      ! ICE THICKNESS (M)
!
!wei added 10/24/2006
      ids(IZORL)  = 4
      ids(IHGT)   = 3
      ids(IVEG)   = 2
      ids(IUST)   = 3
      ids(ICHH)   = 4
      ids(ICMM)   = 4
      ids(ISRF)   = 5
      ids(ITEMP)  = 3
      ids(ISPHUM) = 6
      ids(IZNLW)  = 2
      ids(IMERW)  = 2
      ids(ISNC)   = 3
      ids(ISTC)   = 4
      ids(ISOILM) = 4
      ids(ISNOD)  = 6
      ids(ISNOWD) = 5
      ids(ICNP)   = 5
      ids(IPCPR)  = 6
      ids(ICPCPR) = 6
      ids(IRNOF)  = 5
!
      ILPDS = 28
      IF(ICEN2.EQ.2) ILPDS = 45
      IENS(1) = 1
      IENS(2) = IENST
      IENS(3) = IENSI
      IENS(4) = 1
      IENS(5) = 255
      IYR     = IDATE(4)
      IMO     = IDATE(2)
      IDA     = IDATE(3)
      IHR     = IDATE(1)
      IFHR    = NINT(ZHOUR)
      ITHR    = NINT(FHOUR)
      IF(FHOUR.GT.ZHOUR) THEN
        RTIME = 1./(3600.*(FHOUR-ZHOUR))
      ELSE
        RTIME = 0.
      ENDIF
      IF(SECSWR.GT.0.) THEN
        RTIMSW = 1./SECSWR
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
!     glolal=DUSFC*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
!
      if(me.eq.ioproc) then
!     print *,' ngrid for u flx=',ngrid
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IZWS,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IZWS),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
!         print *, ' called wryte unit noflx' ,noflx,ngrid
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '01)Zonal compt of momentum flux (N/m**2) land and sea surface '
        endif
      endif
 
!..........................................................
!     glolal=DVSFC*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IMWS,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IMWS),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '02)Merid compt of momentum flux (N/m**2) land and sea surface '
        endif
      endif
!..........................................................
!     glolal=DTSFC*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ISHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ISHFLX),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '03)Sensible heat flux (W/m**2) land and sea surface           '
        endif
      endif
!..........................................................
!     glolal=DQSFC*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ILHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ILHFLX),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '04)Latent heat flux (W/m**2) land and sea surface             '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ITEMP,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '05)Temperature (K) land and sea surface                       '
          stop
        endif
      endif
!..........................................................
!     glolal(:,:)=SMC(1,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        LBM=slmskful.EQ.1._kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ISOILM,I2DBLS,0,10,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '06)Volumetric soil moist content (frac) layer 10cm and 0cm    '
        endif
      endif
!..........................................................
!     glolal(:,:)=SMC(2,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        LBM=slmskful.EQ.1._kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
!lu  &            1,ISOILM,I2DBLS,10,200,IYR,IMO,IDA,IHR,
     +            1,ISOILM,I2DBLS,10,40,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
!lu  & '07)Volumetric soil moist content (frac) layer 200cm and 10cm  '
     & '07)Volumetric soil moist content (frac) layer 40cm and 10cm  '
        endif
      endif
!..........................................................
!     glolal(:,:)=STC(1,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        LBM=slmskful.EQ.1._kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ITEMP,I2DBLS,0,10,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '08)Temp (K) layer betw two depth below land sfc 10cm and 0cm  '
        endif
      endif
!..........................................................
!     glolal(:,:)=STC(2,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        LBM=slmskful.EQ.1._kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
!lu  &            1,ITEMP,I2DBLS,10,200,IYR,IMO,IDA,IHR,
     +            1,ITEMP,I2DBLS,10,40,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
!lu  & '09)Temp (K) layer betw two depth below land sfc 200cm and 10cm'
     & '09)Temp (K) layer betw two depth below land sfc 40cm and 10cm'
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ISNOWD,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISNOWD),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '10)Water equiv of accum snow depth (kg/m**2) land sea surface '
        endif
      endif
!..........................................................
!     glolal=DLWSFC*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IDLWF,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IDLWF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '11)Downward long wave radiation flux (W/m**2) land sea surface'
        endif
      endif
!..........................................................
!     glolal=ULWSFC*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IULWF,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IULWF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '12)Upward long wave radiation flux (W/m**2) land sea surface  '
        endif
      endif
!..........................................................
!.......  FIX FLUXES FOR APPROX DIURNAL CYCLE
      DO K=1,4
!       do j=1,LATS_NODE_R
!         do i=1,lonr
!            glolal(i,j)=rflux(i,j,k)*RTIMER(k)
!         enddo
!       enddo
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &              0,IPUR(K),ITLR(K),0,0,IYR,IMO,IDA,IHR,
     &              IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(K)),IENS,
     &              0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          if(k.eq.1)print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '13)Upward long wave radiation flux (W/m**2) top of atmosphere '
          if(k.eq.2)print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '14)Upward solar radiation flux (W/m**2) top of atmosphere     '
          if(k.eq.3)print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '15)Upward solar radiation flux (W/m**2) land and sea surface  '
          if(k.eq.4)print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '16)Downward solar radiation flux (W/m**2) land and sea surface'
        endif
      endif
      ENDDO
!..........................................................
!
!     For UV-B fluxes
!
!     do j=1,LATS_NODE_R
!       do i=1,lonr
!         glolal(i,j)=rflux(i,j,21)*rtimsw
!       enddo
!     enddo
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,129,ICEN,IGEN,
     &            0,IUVBF,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IUVBF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '17)UV-B Downward solar flux (W/m**2) land sea surface'
        endif
      endif
 
!     do j=1,LATS_NODE_R
!       do i=1,lonr
!         glolal(i,j)=rflux(i,j,22)*rtimsw
!       enddo
!     enddo
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,129,ICEN,IGEN,
     &            0,IUVBFC,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IUVBFC),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '18)clear sky UV-B Downward solar flux (W/m**2) land sea surface'
        endif
      endif
!
!     End UV-B fluxes
!
!..........................................................
!..........................................................
      DO K=5,7
!
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!        glolal(i,j)=rflux(i,j,k)*100.*rtimsw
!       enddo
!      enddo
!     where(glolal.ge.0.5)
!       kmskcv=1
!     elsewhere
!       kmskcv=0
!     endwhere
!!
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
!
        K4=4+(K-5)*4
        L=K4+1
        LBM=wrkga.Ge.0.5_kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &              0,IPUR(L),ITLR(L),0,0,IYR,IMO,IDA,IHR,
     &              IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(L)),IENS,
     &              0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          if(k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '19)Total cloud cover (percent) high cloud layer               '
          if(k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '23)Total cloud cover (percent) middle cloud layer             '
          if(k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '27)Total cloud cover (percent) low cloud layer                '
        endif
      endif
!        call baclose(noflx,ierr)
!
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!        if(rflux(i,j,k).gt.0.)then
!         glolal(i,j)=rflux(i,j,k+3)*1000./rflux(i,j,k)
!        else
!         glolal(i,j)=0.
!        endif
!       enddo
!      enddo
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
!
        L=K4+2
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &              1,IPUR(L),ITLR(L),0,0,IYR,IMO,IDA,IHR,
     &              IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(L)),IENS,
     &              0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          if(k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '20)Pressure (Pa) high cloud top level                         '
          if(k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '24)Pressure (Pa) middle cloud top level                       '
          if(k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '28)Pressure (Pa) low cloud top level                          '
        endif
      endif
!
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!        if(rflux(i,j,k).gt.0.)then
!         glolal(i,j)=rflux(i,j,k+6)*1000./rflux(i,j,k)
!        else
!         glolal(i,j)=0.
!        endif
!       enddo
!      enddo
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
!
        L=K4+3
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &              1,IPUR(L),ITLR(L),0,0,IYR,IMO,IDA,IHR,
     &              IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(L)),IENS,
     &              0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          if(k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '21)Pressure (Pa) high cloud bottom level                      '
          if(k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '25)Pressure (Pa) middle cloud bottom level                    '
          if(k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '29)Pressure (Pa) low cloud bottom level                       '
        endif
      endif
!
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!        if(rflux(i,j,k).gt.0.)then
!         glolal(i,j)=rflux(i,j,k+9)/rflux(i,j,k)
!        else
!         glolal(i,j)=0.
!        endif
!       enddo
!      enddo
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        L=K4+4
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &              1,IPUR(L),ITLR(L),0,0,IYR,IMO,IDA,IHR,
     &              IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPUR(L)),IENS,
     &              0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          if(k.eq.5)print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '22)Temperature (K) high cloud top level                       '
          if(k.eq.6)print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '26)Temperature (K) middle cloud top level                     '
          if(k.eq.7)print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '30)Temperature (K) low cloud top level                        '
        endif
      endif
!
      ENDDO
!!
!...................................................................
!     glolal=GESHEM*1.E3*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IPCPR,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IPCPR),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '31)Precipitation rate (kg/m**2/s) land and sea surface        '
        endif
      endif
!...................................................................
!     glolal=BENGSH*1.E3*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ICPCPR,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ICPCPR),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '32)Convective precipitation rate (kg/m**2/s) land sea surface '
        endif
      endif
!...................................................................
!     glolal=GFLUX*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        LBM=slmskful.NE.0._kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IGHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IGHFLX),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '33)Ground heat flux (W/m**2) land and sea surface             '
        endif
      endif
!...................................................................
!     buffo=MOD(slmskloc,2._kind_io8)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ISLMSK,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISLMSK),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '34)Land-sea mask (1=land; 0=sea) (integer) land sea surface   '
        endif
      endif
!...................................................................
!c-- XW: FOR SEA-ICE Nov04
!     buffo=MAX(slmskloc-1._kind_io8,0._kind_io8)
!     call unsplit2z(ioproc,wrkga,global_lats_r)
!     if(me.eq.ioproc) then
!       call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
!    &            0,ICEMSK,ISFC,0,0,IYR,IMO,IDA,IHR,
!    &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICEMSK),IENS,
!    &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
!       if(ierr.ne.0)print*,'wrtsfc gribit ierr=',ierr,'  ',
!    & '35)Ice concentration (ice=1; no ice=0) (1/0) land sea surface '
!     endif
!     IF(IERR.EQ.0 .and. me.eq.ioproc) CALL WRYTE(noflx,LG,G)
!
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ICEMSK,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICEMSK),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     x '35)Ice concentration (ice>0; no ice=0) (1/0) land sea surface '
        endif
      endif

!...................................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IZNLW,IELEV,0,10,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IZNLW),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '36)u wind (m/s) height above ground                           '
        endif
      endif
!...................................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IMERW,IELEV,0,10,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IMERW),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '37)v wind (m/s) height above ground                           '
        endif
      endif
!...................................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ITEMP,IELEV,0,2,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '38)Temperature (K) height above ground                        '
        endif
      endif
!...................................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ISPHUM,IELEV,0,2,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISPHUM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '39)Specific humidity (kg/kg) height above ground              '
        endif
      endif
!...................................................................
!     glolal=PSURF*1.E3
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IPRS,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IPRS),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '40)Pressure (Pa) land and sea surface                         '
        endif
      endif
!...................................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ITMX,IELEV,0,2,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IWIN,0,0,ICEN2,IDS(ITMX),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '41)Maximum temperature (K) height above ground                '
        endif
      endif
!...................................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ITMN,IELEV,0,2,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IWIN,0,0,ICEN2,IDS(ITMN),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '42)Minimum temperature (K) height above ground                '
        endif
      endif
!...................................................................
!     glolal=RUNOFF * 1.E3
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        LBM=slmskful.NE.0._kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IRNOF,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IACC,0,0,ICEN2,IDS(IRNOF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '43)Runoff (kg/m**2) land and sea surface                      '
        endif
      endif
!...................................................................
!     glolal=EP * RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.NE.0._kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IEP,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IEP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '44)Potential evaporation rate (w/m**/) land and sea surface   '
        endif
      endif
!...................................................................
!     glolal=CLDWRK * RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ICLDWK,ICOLMN,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ICLDWK),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '45)Cloud work function (J/Kg) total atmospheric column        '
        endif
      endif
!...................................................................
!     glolal=DUGWD*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IZGW,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IZGW),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '46)Zonal gravity wave stress (N/m**2) land and sea surface    '
        endif
      endif
!...................................................................
!     glolal=DVGWD*RTIME
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IMGW,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IMGW),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '47)Meridional gravity wave stress (N/m**2) land sea surface   '
        endif
      endif
!...................................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IHPBL,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IHPBL),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '48)Boundary layer height '
        endif
      endif
!...................................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IPWAT,ICOLMN,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IPWAT),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '49)Precipitable water (kg/m**2) total atmospheric column      '
        endif
      endif
!...................................................................
!
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!        if (rflux(i,j,4).GT.0.) then
!         glolal(i,j)=rflux(i,j,3)/rflux(i,j,4) * 100.
!         if (glolal(i,j).GT.100.) glolal(i,j)=100.
!        else
!         glolal(i,j)=0.
!        endif
!       enddo
!      enddo
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
!
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IALBDO,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IALBDO),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '50)Albedo (percent) land and sea surface                      '
        endif
      endif
!
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!        glolal(i,j)=rflux(i,j,26)*100.*rtimsw
!       enddo
!      enddo
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
!
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &              0,ICLDF,ICOLMN,0,0,IYR,IMO,IDA,IHR,
     &              IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ICLDF),IENS,
     &              0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '51)Total cloud cover (percent) total atmospheric column       '
        endif
      endif
!
! CONVECTIVE CLOUDS
! LABELED INSTANTANEOUS BUT ACTUALLY AVERAGED OVER FHSWR HOURS
!
!     glolal=CV*1.E2
!     where(glolal.ge.0.5)
!       kmskcv=1
!     elsewhere
!       kmskcv=0
!     endwhere
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=wrkga.Ge.0.5_kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ICLDF,ICVLYR,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICLDF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '52)Total cloud cover (percent) convective cloud layer         '
        endif
      endif
!.................................................
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!       glolal(i,j) = 0.
!       IF(CV(i,j).GT.0.) THEN
!!       ITOP=NINT(CVT(i,j))
!!       IF(ITOP.GE.1.AND.ITOP.LE.LEVS)
!!   &   glolal(i,j)=SI(ITOP+1)*PSURF(i,j)*1.E3
!...      cvt already a pressure (cb)...convert to Pa
!        glolal(i,j)=CVT(i,j)*1.E3
!       END IF
!      ENDDO
!     ENDDO
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IPRS,ICVTL,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IPRS),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '53)Pressure (Pa) convective cloud top level                   '
        endif
      endif
!.................................................
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!       glolal(i,j) = 0.
!       IF(CV(i,j).GT.0.) THEN
!!       Ibot=NINT(CVB(i,j))
!!       IF(Ibot.GE.1.AND.Ibot.LE.LEVS)
!!   &   glolal(i,j)=SI(IBOT)*PSURF(i,j)*1.E3
c...      cvb already a pressure (cb)...convert to Pa
!        glolal(i,j)=CVB(i,j)*1.E3
!       END IF
!      ENDDO
!     ENDDO
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IPRS,ICVBL,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IPRS),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '54)Pressure (Pa) convective cloud bottom level                '
        endif
      endif
!.................................................
!...   SAVE B.L. CLOUD AMOUNT
!
!      do j=1,LATS_NODE_R
!       do i=1,lonr
!        glolal(i,j)=rflux(i,j,27)*100.*rtimsw
!       enddo
!      enddo
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
!
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &              0,ICLDF,IBLLYR,0,0,IYR,IMO,IDA,IHR,
     &              IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ICLDF),IENS,
     &              0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
           print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '55)Total cloud cover (percent) boundary layer cloud layer     '
        endif
      endif
!
!c-- XW: FOR SEA-ICE Nov04
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.2._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ISIK,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISIK),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '56)Sea ice thickness (m) category 1'
        endif
      endif
!c-- XW: END SEA-ICE
!.................................................
!lu: add smc(3:4), stc(3:4), slc(1:4), snwdph, canopy
!lu: addition of 10 records starts here -------------------------------
      if(lsoil.gt.2)then
!     glolal(:,:)=SMC(3,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ISOILM,I2DBLS,40,100,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '57)Volumetric soil moist content (frac) layer 100cm and 40cm '
        endif
      endif
!..........................................................
!     glolal(:,:)=SMC(4,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ISOILM,I2DBLS,100,200,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '58)Volumetric soil moist content (frac) layer 200cm and 100cm '
        endif
      endif
!..........................................................
!     glolal(:,:)=STC(3,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ITEMP,I2DBLS,40,100,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '59)Temp (K) layer betw two depth below land sfc 100cm and 40cm'
        endif
      endif
!..........................................................
!     glolal(:,:)=STC(4,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ITEMP,I2DBLS,100,200,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '60)Temp (K) layer betw two depth below land sfc 200cm and 100cm'
        endif
      endif
      endif
!..........................................................
!     glolal(:,:)=SLC(1,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,130,ICEN,IGEN,
     &            1,ISLC,I2DBLS,0,10,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
         print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '61)Liquid soil moist content (frac) layer 10cm and 0cm  '
        endif
      endif
!..........................................................
!     glolal(:,:)=SLC(2,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,130,ICEN,IGEN,
     &            1,ISLC,I2DBLS,10,40,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '62)Liquid soil moist content (frac) layer 40cm and 10cm '
        endif
      endif
!..........................................................
      if(lsoil.gt.2)then
!     glolal(:,:)=SLC(3,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,130,ICEN,IGEN,
     &            1,ISLC,I2DBLS,40,100,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '63)Liquid soil moist content (frac) layer 100cm and 40cm'
        endif
      endif
!..........................................................
!     glolal(:,:)=SLC(4,:,:)
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,130,ICEN,IGEN,
     &            1,ISLC,I2DBLS,100,200,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISOILM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '64)Liquid soil moist content (frac) layer 200cm and 100cm'
        endif
      endif
      endif
!..........................................................
!     glolal=SNWDPH / 1.E3       !! convert from mm to m
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ISNOD,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISNOD),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
         print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '65)Snow depth (m) land surface                  '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ICNP,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICNP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '66)Canopy water content (kg/m^2) land surface      '
        endif
      endif
!lu: addition of 10 records ends here -------------------------------
!
!wei: addition of 30 records starts here -------------------------------
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IZORL,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IZORL),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '67)Surface roughness (m)       '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IVEG,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IVEG),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '68)Vegetation fraction (fractional) land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IVTP,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IVTP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '69)Vegetation type land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ISTP,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISTP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '70)Soil type land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,130,ICEN,IGEN,
     &            1,ISLO,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISLO),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '71)Slope type land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IUST,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IUST),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '72)Frictional velocity (m/s)      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IHGT,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IHGT),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '73)Surface height (m)      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IRST,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IRST),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '74)Freezing precip flag land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ICHH,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICHH),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '75)Exchange coefficient CH(m/s)      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,130,ICEN,IGEN,
     &            0,ICMM,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ICMM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '76)Exchange coefficient CM(m/s)      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IEP,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IEP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '77)Potential evaporation rate (w/m**2) land and sea surface   '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IDLWF,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IDLWF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '78)Downward long wave radiation flux (W/m**2)'
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IULWF,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IULWF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '79)Upward long wave radiation flux (W/m**2)  '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IUSWF,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IUSWF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '80)Upward short wave radiation flux (W/m**2)  '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IDSWF,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IDSWF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '81)Downward short wave radiation flux (W/m**2)   '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ISHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISHFLX),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '82)Sensible heat flux (W/m**2) land and sea surface           '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ILHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ILHFLX),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '83)Latent heat flux (W/m**2) land and sea surface             '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
        LBM=slmskful.NE.0._kind_io8
        call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IGHFLX,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IGHFLX),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '84)Ground heat flux (W/m**2) land and sea surface             '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ISRF,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IACC,0,0,ICEN2,IDS(ISRF),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '85)Surface runoff (kg/m^2) land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ITEMP,isglev,1,1,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ITEMP),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '86)Lowest model level Temp (K)       '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,ISPHUM,isglev,1,1,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISPHUM),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '87)Lowest model specific humidity (kg/kg)       '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IZNLW,isglev,1,1,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IZNLW),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '88)Lowest model u wind (m/s)      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            0,IMERW,isglev,1,1,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IMERW),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '89)Lowest model v wind (m/s)       '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IHGT,isglev,1,1,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(IHGT),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '90)Lowest model level height (m) land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IEVBS,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IEVBS),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '91)Direct evaporation from bare soil(W/m^2) land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,IEVCW,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(IEVBS),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '92)Canopy water evaporation(W/m^2) land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ITRAN,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ITRAN),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '93)Transpiration (W/m^2) land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,130,ICEN,IGEN,
     &            1,ISBS,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ISBS),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '94)Snow Sublimation (W/m^2) land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ISNC,ISFC,0,0,IYR,IMO,IDA,IHR,
     &            IFHOUR,IFHR,ITHR,IAVG,0,0,ICEN2,IDS(ISNC),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '95)Snow Cover (fraction) land surface      '
        endif
      endif
!..........................................................
      call unsplit2z(ioproc,wrkga,global_lats_r)
      if(me.eq.ioproc) then
      LBM=slmskful.EQ.1._kind_io8
      call gribit(wrkga,LBM,4,lonr,latr,16,CL1,ILPDS,2,ICEN,IGEN,
     &            1,ISTC,I2DBLS,0,200,IYR,IMO,IDA,IHR,
     &            IFHOUR,ITHR,0,INST,0,0,ICEN2,IDS(ISTC),IENS,
     &            0.,0.,0.,0.,0.,0.,G,LG,IERR)
        IF(IERR.EQ.0) then
          CALL WRYTE(noflx,LG,G)
        else
          print*,'wrtsfc gribit ierr=',ierr,'  ',
     & '96)Total column soil moisture (Kg/m^2) land surface      '
        endif
      endif
!wei: addition of 30 records ends here -------------------------------
!!
      if(me.eq.ioproc)
     &   PRINT *,'(wrtflx_w) GRIB FLUX FILE WRITTEN ',FHOUR,IDATE,noflx
!!
      RETURN
      END
       subroutine grids_move(ioproc)
!
!***********************************************************************
!
      use resol_def
      use mod_state
      use layout1
      use mpi_def
      implicit none
!
      integer ipt_lats_node_rl,nodesr
      integer lats_nodes_rl
!      integer lats_nodes_r(nodes),ipt,maxfld,ioproc,nproct
      integer ioproc
      integer proc,j,lat,msgtag,nproc,i,msgtag1,buff,startlat,ierr
      integer illen,ubound
       integer icount
       data icount/0/
         integer maxlats_comp
          save maxlats_comp
       integer kllen
!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
         if(icount .eq. 0) then
          allocate(ivar_global(10))
          allocate(ivar_global_a(10,nodes))
         ivar_global(1)=ipt_lats_node_r
         ivar_global(2)= lats_node_r
         ivar_global(3)=lats_node_r_max
         call mpi_gather(ivar_global,10,MPI_INTEGER,
     1 ivar_global_a,10,MPI_INTEGER,ioproc,MPI_COMM_ALL,ierr)
         icount=icount+1
         endif
!!
       if(allocated(buff_mult_pieces)) then
          continue
          else
            maxlats_comp=lats_node_r_max
            if(.not. liope .or. me .ne. ioproc) then
             continue
            else
!            maxlats_comp=ivar_global_a(3,ioproc)
            maxlats_comp=ivar_global_a(3,1)
            endif
!gwv watch this!!
          allocate
     1  (buff_mult_pieces(lonr,ngrids_sfcc,maxlats_comp,nodes))
          allocate
     1  (buff_mult_piecesf(lonr,0:ngrids_flx,maxlats_comp,nodes))
          allocate
     1  (buff_mult_piecesa(lonr,1:ngrids_flx+1+ngrids_sfcc,
     1  maxlats_comp,nodes))
          endif
 
!
!  big send
        IF (me.ne.ioproc) THEN
!
!         Sending the data
         msgtag=me
         illen=lats_node_r
              kllen=illen*lonr*(ngrids_flx+1+ngrids_sfcc)
! send the local grid domain
         CALL mpi_send
     &(buff_mult_piecea,kllen,MPI_R_IO,ioproc,
     &                  msgtag,MPI_COMM_ALL,info)
      ELSE
         if( MC_COMP .ne. MPI_COMM_NULL) then
! iotask is also a compute task.  send is replaced with direct
!  array copy
 
           buff_mult_piecesa(:,:,1:lats_node_r,ioproc+1)=
     1   buff_mult_piecea(:,:,1:lats_node_r)
!  END COMPUTE TASKS PORTION OF LOGIC
         endif
!  END COMPUTE TASKS PORTION OF LOGIC
!  receiving part of I/O task
 
!!
!!     for pes ioproc
        DO proc=1,nodes_comp
         if (proc.ne.ioproc+1) then
         msgtag=proc-1
         illen=ivar_global_a(2,proc)
          kllen=illen*lonr*(ngrids_flx+1+ngrids_sfcc)
          CALL mpi_recv
     1 (buff_mult_piecesa(1,1,1,proc),kllen
     1 ,MPI_R_IO,proc-1,
     &                msgtag,MPI_COMM_ALL,stat,info)
         endif
        enddo
      ENDIF
 
!!
      return
      end
 
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
 
 
