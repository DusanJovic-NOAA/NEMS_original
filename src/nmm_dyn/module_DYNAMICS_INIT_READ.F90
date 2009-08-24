                        module module_DYNAMICS_INIT_READ
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
use module_include
use module_dm_parallel,only : ids,ide,jds,jde &
                             ,ims,ime,jms,jme &
                             ,its,ite,jts,jte &
                             ,lm &
                             ,mype_share,npes,num_pts_max &
                             ,mpi_comm_comp &
                             ,dstrb
use module_exchange
use module_constants
use nemsio_module
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!---look-up tables------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint),private :: &
 mype
character(64):: &
 infile
real(kind=kfpt),allocatable,dimension(:,:):: &      !im,jm
 stdh                        ! standard deviation of topography height
integer(kind=kint):: &
ihrend &                    ! maximum forecast length, hours
,ntsd &
,ntstm_max &
,nfcst
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
       contains
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine dynamics_read_binary &
      (global &
      ,kss,kse &
      ,pdtop,pt,lpt2 &
      ,sg1,dsg1 &
      ,psg1,pdsg1 &
      ,sg2,dsg2,sgm &
      ,sgml1,psgml1,sgml2 &
      ,fis,sm,sice &
      ,pd,pdo,pint &
      ,u,v,q2,e2 &
      ,t,q,cw &
      ,tp,up,vp &
      ,rrw,dwdt,w &
      ,omgalf,div,z &
      ,rtop &
      ,tcu,tcv,tct &
      ,sp,indx_q,indx_cw,indx_rrw,indx_q2 &
      ,ntsti,ntstm &
      ,ihr,ihrst,idat &
      ,run,restart &
      ,num_water,water &
      ,num_tracers_total,tracers &
      ,p_qv,p_qc,p_qr &
      ,p_qi,p_qs,p_qg &
      ,dt,nhours_fcst)
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 kse &
,kss &
,num_water &
,num_tracers_total &
,p_qv &
,p_qc &
,p_qr &
,p_qi &
,p_qs &
,p_qg &
,indx_q &
,indx_cw &
,indx_rrw &
,indx_q2 &
,nhours_fcst

real(kind=kfpt),intent(in):: &
dt

integer(kind=kint),intent(out) :: &
 ihr &
,ihrst &
,lpt2 &
,ntsti &
,ntstm

integer(kind=kint),dimension(3),intent(out) :: &
 idat

real(kind=kfpt),intent(out) :: &
 pdtop &
,pt

real(kind=kfpt),dimension(1:lm),intent(out) :: &
 dsg1 &
,dsg2 &
,pdsg1 &
,psgml1 &
,sgml1 &
,sgml2

real(kind=kfpt),dimension(1:lm+1),intent(out) :: &
 psg1 &
,sg1 &
,sg2 &
,sgm

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out) :: &
 fis &
,pd &
,pdo &
,sice &
,sm

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out) :: &
 cw &
,dwdt &
,q &
,q2 &
,rrw &
,omgalf &
,div &
,z &
,rtop &
,tcu &
,tcv &
,tct &
,t &
,tp &
,u &
,up &
,v &
,vp &
,e2 &
,w

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(out) :: &
 pint

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:num_water),intent(out) :: &
 water

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:num_tracers_total),intent(out) :: &
 tracers


real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,kss:kse),intent(out) :: &
 sp

logical(kind=klog),intent(in) :: &
 global &
,restart

logical(kind=klog),intent(out) :: &
 run
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &                         ! index in x direction
,irtn &
,j &                         ! index in y direction
,k &                         ! index
,ks &                        ! tracer index
,l &                         ! index in p direction
,n &
,nrecs_skip_for_pt

integer(kind=kint) :: &      ! number of soil levels
 nsoil

integer(kind=kint) :: &
 iyear_fcst           &
,imonth_fcst          &
,iday_fcst            &
,ihour_fcst

real(kind=kfpt):: &
 tend,tend_max

real(kind=kfpt),allocatable,dimension(:,:) :: &
 temp1

logical(kind=klog) :: opened

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
      allocate(temp1(ids:ide,jds:jde),stat=i)
      allocate(stdh(ims:ime,jms:jme),stat=i)
!
!-----------------------------------------------------------------------

      select_unit: do n=51,59
        inquire(n,opened=opened)
        if(.not.opened)then
          nfcst=n
          exit select_unit
        endif
      enddo select_unit

!-----------------------------------------------------------------------
         if(.not.restart) then                     ! cold start
!-----------------------------------------------------------------------
!
      infile='main_input_filename'
      open(unit=nfcst,file=infile,status='old',form='unformatted')
!
!-----------------------------------------------------------------------
!***  First we need to extract the value of PT (pressure at top of domain)
!-----------------------------------------------------------------------
!
      if(mype==0)then
!        nrecs_skip_for_pt=6+5*lm+21 !<-- For current WPS input
       nrecs_skip_for_pt=6+5*lm+23 !zj 21
!
        do n=1,nrecs_skip_for_pt
          read(nfcst)
        enddo
!
        read(nfcst)pt
        write(0,*)' PT from input file equals ',pt
        rewind nfcst
      endif
      call mpi_bcast(pt,1,mpi_real,0,mpi_comm_comp,irtn)
!
!-----------------------------------------------------------------------
!
      read(nfcst) run,idat,ihrst,ihrend,ntsd
!     read(nfcst)

!     run=.true.

      if(mype==0)then
        write(0,*) 'run, idat,ntsd: ', run, idat, ntsd
      endif
!
!-----------------------------------------------------------------------
!***  Print the time information.
!-----------------------------------------------------------------------
!
      if(mype==0)then
        write(0,*)' Start year =',idat(3)
        write(0,*)' Start month=',idat(2)
        write(0,*)' Start day  =',idat(1)
        write(0,*)' Start hour =',ihrst
        write(0,*)' Timestep   =',dt
        write(0,*)' Steps/hour =',3600./dt
        if(.not.global)write(0,*)' Max fcst hours=',ihrend
      endif
!
!-----------------------------------------------------------------------
!
      read (nfcst) pt,pdtop,lpt2,sgm,sg1,dsg1,sgml1,sg2,dsg2,sgml2
!
      if(mype==0)then
        write(0,*) 'nmm_libutil reads of PT, PDTOP: ', PT, PDTOP
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        fis(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,fis,1,1,1,1,1)
      call halo_exch(fis,1,2,2)
!-----------------------------------------------------------------------
!
      if(mype==0)then
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        stdh(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,stdh,1,1,1,1,1)
      call halo_exch(stdh,1,2,2)
!-----------------------------------------------------------------------
!
      if(mype==0)then
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        sm(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,sm,1,1,1,1,1)
      call halo_exch(sm,1,2,2)
!-----------------------------------------------------------------------
!
      if(mype==0)then
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        pd(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,pd,1,1,1,1,1)
      call halo_exch(pd,1,2,2)
!-----------------------------------------------------------------------
!
      call mpi_barrier(mpi_comm_comp,irtn)
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          u(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,u,1,1,1,lm,l)
      enddo
      call halo_exch(u,lm,2,2)
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          v(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,v,1,1,1,lm,l)
      enddo
      call halo_exch(v,lm,2,2)
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
          write(0,*) 'L, T extremes: ', L, minval(temp1),maxval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          t(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,t,1,1,1,lm,l)
      enddo
      call halo_exch(t,lm,2,2)
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          q(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,q,1,1,1,lm,l)
      enddo
      call halo_exch(q,lm,2,2)
!
      do l=1,lm
      do j=jms,jme
      do i=ims,ime
        water(i,j,l,p_qv)=q(i,j,l)/(1.-q(i,j,l))    ! WRF water array uses mixing ratio for vapor
      enddo
      enddo
      enddo
!
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          cw(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,cw,1,1,1,lm,l)
      enddo
      call halo_exch(cw,lm,2,2)
!
      print*,'*** Read initial conditions in init from ',infile
!
!      close(unit=nfcst)
!-----------------------------------------------------------------------
      ntsti=ntsd+1
!
!
      tend_max=real(ihrend)
      ntstm_max=nint(tend_max*3600./dt)+1
      tend=real(nhours_fcst)
      ntstm=nint(tend*3600./dt)+1
      if(.not.global)then
         write(0,*)' Max runtime is ',tend_max,' hours'
      endif
      if(mype==0)then
        write(0,*)' Requested runtime is ',tend,' hours'
        write(0,*)' NTSTM=',ntstm
      endif
      if(ntstm>ntstm_max.and..not.global)then
        if(mype==0)then
          write(0,*)' Requested fcst length exceeds maximum'
          write(0,*)' Resetting to maximum'
        endif
        ntstm=min(ntstm,ntstm_max)
      endif
!
      ihr=nint(ntsd*dt/3600.)
!-----------------------------------------------------------------------
      do l=1,lm
        pdsg1(l)=dsg1(l)*pdtop
        psgml1(l)=sgml1(l)*pdtop+pt
!       write(0,*) 'L, pdsg1, psgml1: ', L, pdsg1(L), psgml1(L)
      enddo
!
      do l=1,lm+1
        psg1(l)=sg1(l)*pdtop+pt
      enddo
!-----------------------------------------------------------------------
      do j=jts,jte
        do i=its,ite
          pdo(i,j)=pd(i,j)
        enddo
      enddo
      call halo_exch(pdo,1,2,2)
!
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            up(i,j,l)=u(i,j,l)
            vp(i,j,l)=v(i,j,l)
            tp(i,j,l)=t(i,j,l)
          enddo
        enddo
      enddo
      call halo_exch(tp,lm,up,lm,vp,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            q2(i,j,l)=0.02
            rrw(i,j,l)=0.
            if(i.ge.ide  /2+1- 6.and.i.le.ide  /2+1+ 6.and. &
               j.ge.jde*3/4+1- 6.and.j.le.jde*3/4+1+ 6.) then !global
!               j.ge.jde  /2+1- 6.and.j.le.jde  /2+1+ 6.) then !regional
              rrw(i,j,l)=10.
            endif
            dwdt(i,j,l)=1.
            w(i,j,l)=0.
          enddo
        enddo
      enddo
      call halo_exch(dwdt,lm,2,2)
!
      do j=jts,jte
        do i=its,ite
          pint(i,j,1)=pt
        enddo
      enddo
!
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            pint(i,j,l+1)=pint(i,j,l)+dsg2(l)*pd(i,j)+pdsg1(l)
          enddo
        enddo
      enddo
      call halo_exch(pint,lm+1,2,2)
!
      call halo_exch(q2,lm,rrw,lm,2,2)
      do l=1,lm
        do j=jms,jme
          do i=ims,ime
            sp(i,j,l,indx_q  )=sqrt(max(q  (i,j,l),0.))
            sp(i,j,l,indx_cw )=sqrt(max(cw (i,j,l),0.))
            sp(i,j,l,indx_rrw)=sqrt(max(rrw(i,j,l),0.))
            sp(i,j,l,indx_q2 )=sqrt(max(q2 (i,j,l),0.))
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!---reading surface data------------------------------------------------
!-----------------------------------------------------------------------
!
       if(mype==0)then
         read(nfcst)temp1  ! ALBEDO
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! EPSR
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! MXSNAL
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! TSKIN
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! SST
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! SNO
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! SI
       endif
!
       if(mype==0)then
         read(nfcst)temp1  ! SICE
       endif
!      write(0,*) 'maxval(sice): ', maxval(temp1)
       do j=jms,jme
       do i=ims,ime
         sice(i,j)=0.
       enddo
       enddo
       call dstrb(temp1,sice,1,1,1,1,1)
       call halo_exch(sice,1,2,2)

!
!-----------------------------------------------------------------------
!
      if(mype==0)then
!        do n=1,13 !<-- For current WPS input
       do n=1,15 !zj 13
          read(nfcst)
        enddo
      endif
!
      if(mype==0)read(nfcst)pt
      write(0,*) 'pt read again here...pt: ', pt
      call mpi_bcast(pt,1,mpi_real,0,mpi_comm_comp,irtn)
!
!-----------------------------------------------------------------------
!
      if(mype==0)then
        close(nfcst)
      endif
!
!-----------------------------------------------------------------------
         else                                      ! restart
!-----------------------------------------------------------------------
!
      infile='restart_file'
      open(unit=nfcst,file=infile,status='old',form='unformatted')
!
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: INGEGER SCALARS
!-----------------------------------------------------------------------
      read(nfcst) iyear_fcst
      read(nfcst) imonth_fcst
      read(nfcst) iday_fcst
      read(nfcst) ihour_fcst
      read(nfcst) ! iminute_fcst
      read(nfcst) ! second_fcst
      read(nfcst) ! ntsd
      read(nfcst) ! im
      read(nfcst) ! jm
      read(nfcst) ! lm
      read(nfcst) ihrst
      read(nfcst) lpt2
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: INTEGER 1D ARRAYS
!-----------------------------------------------------------------------
      read(nfcst) idat
!-----------------------------------------------------------------------
!***  Print the time information.
!-----------------------------------------------------------------------
!
      if(mype==0)then
        write(0,*)'**** read in dynamics ***************'
        write(0,*)' Restart year =',iyear_fcst
        write(0,*)' Restart month=',imonth_fcst
        write(0,*)' Restart day  =',iday_fcst
        write(0,*)' Restart hour =',ihour_fcst
        write(0,*)' Original start year =',idat(3)
        write(0,*)' Original start month=',idat(2)
        write(0,*)' Original start day  =',idat(1)
        write(0,*)' Original start hour =',ihrst
        write(0,*)' Timestep   =',dt
        write(0,*)' Steps/hour =',3600./dt
        write(0,*)'*************************************'
      endif
!
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: INTEGER SCALARS
!-----------------------------------------------------------------------
      read(nfcst) ! mp_physics
      read(nfcst) ! sf_surface_physics
      read(nfcst) nsoil
      read(nfcst) ! nphs
      read(nfcst) ! nclod
      read(nfcst) ! nheat
      read(nfcst) ! nprec
      read(nfcst) ! nrdlw
      read(nfcst) ! nrdsw
      read(nfcst) ! nsrfc
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL SCALARS
!-----------------------------------------------------------------------
      read(nfcst) ! dt
      read(nfcst) ! dyh
      read(nfcst) pdtop
!-----------------------------------------------------------------------
      read(nfcst) pt
      read(nfcst) ! tlm0d
      read(nfcst) ! tph0d
      read(nfcst) ! tstart
      read(nfcst) ! dphd
      read(nfcst) ! dlmd
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 1D ARRAYS
!-----------------------------------------------------------------------
      read(nfcst) ! dxh
      read(nfcst) sg1
      read(nfcst) sg2
      read(nfcst) dsg1
      read(nfcst) dsg2
      read(nfcst) sgml1
      read(nfcst) sgml2
      read(nfcst) sgm
      read(nfcst) ! aphtim
      read(nfcst) ! ardlw
      read(nfcst) ! ardsw
      read(nfcst) ! asrfc
      read(nfcst) ! avcnvc
      read(nfcst) ! avrain
      read(nfcst) ! sldpth
      read(nfcst) ! mp_restart_state
      read(nfcst) ! tbpvs_state
      read(nfcst) ! tbpvs0_state
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: LOGICAL
!-----------------------------------------------------------------------
      read(nfcst) ! global
      read(nfcst) run
      if(mype==0)then
        write(0,*) 'run, idat,ntsd: ', run, idat, ntsd
      endif
      read(nfcst) ! adiabatic
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: INTEGER 2D ARRAYS
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst) ! isltyp
      endif
      if(mype==0)then
        read(nfcst) ! ivgtyp
      endif
      if(mype==0)then
        read(nfcst) ! ncfrcv
      endif
      if(mype==0)then
        read(nfcst) ! ncfrst
      endif
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 2D ARRAYS
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        fis(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,fis,1,1,1,1,1)
      call halo_exch(fis,1,2,2)
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst) !glat
        read(nfcst) !glon
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        pd(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,pd,1,1,1,1,1)
      call halo_exch(pd,1,2,2)
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst) !vlat
        read(nfcst) !vlon
        read(nfcst)temp1
      endif
      do j=jms,jme
      do i=ims,ime
        pdo(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,pdo,1,1,1,1,1)
      call halo_exch(pdo,1,2,2)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 3D ARRAYS (only DYN)
!-----------------------------------------------------------------------
!     call mpi_barrier(mpi_comm_comp,irtn)
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          w(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,w,1,1,1,lm,l)
      enddo
      call halo_exch(w,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          dwdt(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,dwdt,1,1,1,lm,l)
      enddo
      call halo_exch(dwdt,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm+1
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          pint(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,pint,1,1,1,lm+1,l)
      enddo
      call halo_exch(pint,lm+1,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          omgalf(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,omgalf,1,1,1,lm,l)
      enddo
      call halo_exch(omgalf,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          rrw(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,rrw,1,1,1,lm,l)
      enddo
      call halo_exch(rrw,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          div(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,div,1,1,1,lm,l)
      enddo
      call halo_exch(div,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          rtop(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,rtop,1,1,1,lm,l)
      enddo
      call halo_exch(rtop,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          tcu(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,tcu,1,1,1,lm,l)
      enddo
      call halo_exch(tcu,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          tcv(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,tcv,1,1,1,lm,l)
      enddo
      call halo_exch(tcv,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          tct(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,tct,1,1,1,lm,l)
      enddo
      call halo_exch(tct,lm,2,2)
!      write(0,*)'min,max tct=',minval(tct),maxval(tct)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          tp(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,tp,1,1,1,lm,l)
      enddo
      call halo_exch(tp,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          up(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,up,1,1,1,lm,l)
      enddo
      call halo_exch(up,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          vp(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,vp,1,1,1,lm,l)
      enddo
      call halo_exch(vp,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          e2(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,e2,1,1,1,lm,l)
      enddo
      call halo_exch(e2,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm-1
        if(mype==0)then
          read(nfcst)temp1
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          z(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,z,1,1,1,lm,l)
      enddo
      call halo_exch(z,lm,2,2)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 2D ARRAYS (contd.)
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ACFRCV
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ACFRST
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ACPREC
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ACSNOM
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ACSNOW
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! AKHS_OUT
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! AKMS_OUT
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ALBASE
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! albedo
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ALWIN
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ALWOUT
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ALWTOA
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ASWIN
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ASWOUT
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! ASWTOA
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! BGROFF
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CFRACH
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CFRACL
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CFRACM
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! cldefi
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! cmc
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CNVBOT
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CNVTOP
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CPRATE
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CUPPT
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CUPREC
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CZEN
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! CZMEAN
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! epsr
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! GRNFLX
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! HBOTD
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! HBOTS
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! HTOPD
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! HTOPS
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! mxsnal
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! PBLH
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! POTEVP
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! PREC
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! PSHLTR
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! Q10
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! QSH
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! QSHLTR
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! QWBS
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! QZ0
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! RADOT
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! RLWIN
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! RLWTOA
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! RSWIN
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! RSWINC
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! RSWOUT
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SFCEVP
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SFCEXC
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SFCLHX
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SFCSHX
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! si
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! sice
      endif
      do j=jms,jme
      do i=ims,ime
        sice(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,sice,1,1,1,1,1)
      call halo_exch(sice,1,2,2)
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SIGT4
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! sm
      endif
      do j=jms,jme
      do i=ims,ime
        sm(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,sm,1,1,1,1,1)
      call halo_exch(sm,1,2,2)
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SMSTAV
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SMSTOT
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! sno
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SNOPCX
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SOILTB
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! sr
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SSROFF
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! sst
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! SUBSHX
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! tg
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! TH10
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! THS
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! THZ0
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! TSHLTR
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! TWBS
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! U10
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  !ustar
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! UZ0
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! V10
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  !vegfrc
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! VZ0
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  !z0
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! tskin
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! akhs
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! akms
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! hbot
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! htop
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! rswtoa
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! potflx
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! rmol
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! t2
      endif
!-----------------------------------------------------------------------
      if(mype==0)then
        read(nfcst)temp1  ! z0base
      endif
      if(mype==0)then
        read(nfcst)temp1  ! tlmin
      endif
      if(mype==0)then
        read(nfcst)temp1  ! tlmax
      endif
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 3D ARRAYS
!-----------------------------------------------------------------------
      call mpi_barrier(mpi_comm_comp,irtn)
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1 ! cldfra
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          cw(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,cw,1,1,1,lm,l)
      enddo
      call halo_exch(cw,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          q(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,q,1,1,1,lm,l)
      enddo
      call halo_exch(q,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          q2(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,q2,1,1,1,lm,l)
      enddo
      call halo_exch(q2,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1 ! rlwtt
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1 ! rswtt
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
          write(0,*) 'L, T extremes: ', L, minval(temp1),maxval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          t(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,t,1,1,1,lm,l)
      enddo
      call halo_exch(t,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1 ! tcucn
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1 ! train
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          u(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,u,1,1,1,lm,l)
      enddo
      call halo_exch(u,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          v(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,v,1,1,1,lm,l)
      enddo
      call halo_exch(v,lm,2,2)
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1 ! xlen_mix
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1 ! f_ice
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1 ! f_rimef
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1 ! f_rain
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm+1
        if(mype==0)then
          read(nfcst)temp1 ! rthblten
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,lm+1
        if(mype==0)then
          read(nfcst)temp1 ! rqvblten
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,nsoil
        if(mype==0)then
          read(nfcst)temp1 ! sh2o
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,nsoil
        if(mype==0)then
          read(nfcst)temp1 ! smc
        endif
      enddo
!-----------------------------------------------------------------------
      do l=1,nsoil
        if(mype==0)then
          read(nfcst)temp1 ! stc
        endif
      enddo
!-----------------------------------------------------------------------
      do n=indx_rrw+1,num_tracers_total
      do l=1,lm
        if(mype==0)then
          read(nfcst)temp1
        endif
        do j=jms,jme
        do i=ims,ime
          tracers(i,j,l,n)=0.
        enddo
        enddo
        call dstrb(temp1,tracers(:,:,:,n),1,1,1,lm,l)
      enddo
      enddo
!
      call halo_exch(tracers,lm,num_tracers_total,1,2,2)
!
      do n=1,num_water
      do l=1,lm
        do j=jms,jme
        do i=ims,ime
          water(i,j,l,n)=tracers(i,j,l,n+num_tracers_total-num_water)
        enddo
        enddo
      enddo
      enddo
!
!-----------------------------------------------------------------------
!
      if(mype==0)then
        close(nfcst)
      endif
!
!-----------------------------------------------------------------------
      ntsti=ntsd+1
!
      tend_max=real(ihrend)
      ntstm_max=nint(tend_max*3600./dt)+1
      tend=real(nhours_fcst)
      ntstm=nint(tend*3600./dt)+1
      if(.not.global)then
         write(0,*)' Max runtime is ',tend_max,' hours'
      endif
      if(mype==0)then
        write(0,*)' Requested runtime is ',tend,' hours'
        write(0,*)' NTSTM=',ntstm
      endif
      if(ntstm>ntstm_max.and..not.global)then
        if(mype==0)then
          write(0,*)' Requested fcst length exceeds maximum'
          write(0,*)' Resetting to maximum'
        endif
        ntstm=min(ntstm,ntstm_max)
      endif
!
      ihr=nint(ntsd*dt/3600.)
!-----------------------------------------------------------------------
      do l=1,lm
        pdsg1(l)=dsg1(l)*pdtop
        psgml1(l)=sgml1(l)*pdtop+pt
      enddo
!
      do l=1,lm+1
        psg1(l)=sg1(l)*pdtop+pt
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jms,jme
          do i=ims,ime
            sp(i,j,l,indx_q  )=sqrt(max(q  (i,j,l),0.))
            sp(i,j,l,indx_cw )=sqrt(max(cw (i,j,l),0.))
            sp(i,j,l,indx_rrw)=sqrt(max(rrw(i,j,l),0.))
            sp(i,j,l,indx_q2 )=sqrt(max(q2 (i,j,l),0.))
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
         endif                                     ! cold start /restart
!-----------------------------------------------------------------------
!
      deallocate(temp1)
      deallocate(stdh)
      if(mype==0)then
        write(0,*)' EXIT SUBROUTINE INIT pt=',pt
      endif
!-----------------------------------------------------------------------
!
                        end subroutine dynamics_read_binary
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
                        subroutine dynamics_read_nemsio &
      (global &
      ,kss,kse &
      ,pdtop,pt,lpt2 &
      ,sg1,dsg1 &
      ,psg1,pdsg1 &
      ,sg2,dsg2,sgm &
      ,sgml1,psgml1,sgml2 &
      ,fis,sm,sice &
      ,pd,pdo,pint &
      ,u,v,q2,e2 &
      ,t,q,cw &
      ,tp,up,vp &
      ,rrw,dwdt,w &
      ,omgalf,div,z &
      ,rtop &
      ,tcu,tcv,tct &
      ,sp,indx_q,indx_cw,indx_rrw,indx_q2 &
      ,ntsti,ntstm &
      ,ihr,ihrst,idat &
      ,run,restart &
      ,num_water,water &
      ,num_tracers_total,tracers &
      ,p_qv,p_qc,p_qr &
      ,p_qi,p_qs,p_qg &
      ,dt,nhours_fcst)
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none

!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 kse &
,kss &
,num_water &
,num_tracers_total &
,p_qv &
,p_qc &
,p_qr &
,p_qi &
,p_qs &
,p_qg &
,indx_q &
,indx_cw &
,indx_rrw &
,indx_q2 &
,nhours_fcst

real(kind=kfpt),intent(in):: &
dt


integer(kind=kint),intent(out) :: &
 ihr &
,ihrst &
,lpt2 &
,ntsti &
,ntstm

integer(kind=kint),dimension(3),intent(out) :: &
 idat

real(kind=kfpt),intent(out) :: &
 pdtop &
,pt

real(kind=kfpt),dimension(1:lm),intent(out) :: &
 dsg1 &
,dsg2 &
,pdsg1 &
,psgml1 &
,sgml1 &
,sgml2

real(kind=kfpt),dimension(1:lm+1),intent(out) :: &
 psg1 &
,sg1 &
,sg2 &
,sgm

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out) :: &
 fis &
,pd &
,pdo &
,sice &
,sm

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(out) :: &
 cw &
,dwdt &
,q &
,q2 &
,rrw &
,omgalf &
,div &
,z &
,rtop &
,tcu &
,tcv &
,tct &
,t &
,tp &
,u &
,up &
,v &
,vp &
,e2 &
,w

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(out) :: &
 pint

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:num_water),intent(out) :: &
 water

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,1:num_tracers_total),intent(out) :: &
 tracers


real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm,kss:kse),intent(out) :: &
 sp

logical(kind=klog),intent(in) :: &
 global &
,restart

logical(kind=klog),intent(out) :: &
 run
!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &                         ! index in x direction
,irtn &
,j &                         ! index in y direction
,k &                         ! index
,ks &                        ! tracer index
,l &                         ! index in p direction
,n &
,nrecs_skip_for_pt &
,ierr

integer(kind=kint) :: &      ! number of soil levels
 nsoil

integer(kind=kint) :: &
 iyear_fcst           &
,imonth_fcst          &
,iday_fcst            &
,ihour_fcst
!
!jws
integer(kind=kint) :: idate(7)
integer(kind=kint) :: fcstdate(7)
character(2)       ::tn
!jwe
real(kind=kfpt):: &
 tend,tend_max

real(kind=kfpt),allocatable,dimension(:) :: &
 temp1


logical(kind=klog) :: opened

type(nemsio_gfile) :: gfile
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
      allocate(temp1((ide-ids+1)*(jde-jds+1)),stat=i)
      allocate(stdh(ims:ime,jms:jme),stat=i)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!***  initiaizel nemsio_module
!-----------------------------------------------------------------------
!
      call nemsio_init()
!-----------------------------------------------------------------------
!
      if(.not.restart) then                     ! cold start
!-----------------------------------------------------------------------
!
      infile='main_input_filename_nemsio'
      call nemsio_open(gfile,infile,'read',iret=ierr)
      if(ierr/=0) write(0,*)'ERROR: open file ',trim(infile),' has failed'
!
!-----------------------------------------------------------------------
!***  First we need to extract the value of PT (pressure at top of domain)
!-----------------------------------------------------------------------
!
      call nemsio_getheadvar(gfile,'PT',pt,ierr)
!
!-----------------------------------------------------------------------
!*** get run,idat,ihrst,ihrend,ntsd
!
      call nemsio_getheadvar(gfile,'run',run,ierr)
      call nemsio_getheadvar(gfile,'IDAT',idat,ierr)
      call nemsio_getheadvar(gfile,'ihrst',ihrst,ierr)
      call nemsio_getheadvar(gfile,'ihrend',ihrend,ierr)
      call nemsio_getheadvar(gfile,'ntsd',ntsd,ierr)

      if(mype==0)then
        write(0,*) 'run, idat,ntsd: ', run, idat, ntsd
      endif
!
!-----------------------------------------------------------------------
!***  Print the time information.
!-----------------------------------------------------------------------
!
      if(mype==0)then
        write(0,*)' Start year =',idat(3)
        write(0,*)' Start month=',idat(2)
        write(0,*)' Start day  =',idat(1)
        write(0,*)' Start hour =',ihrst
        write(0,*)' Timestep   =',dt
        write(0,*)' Steps/hour =',3600./dt
        if(.not.global)write(0,*)' Max fcst hours=',ihrend
      endif
!
!-----------------------------------------------------------------------
!
      call nemsio_getheadvar(gfile,'pt',pt,ierr)
      call nemsio_getheadvar(gfile,'pdtop',pdtop,ierr)
      call nemsio_getheadvar(gfile,'lpt2',lpt2,ierr)
      call nemsio_getheadvar(gfile,'sgm',sgm,ierr)
      call nemsio_getheadvar(gfile,'sg1',sg1,ierr)
      call nemsio_getheadvar(gfile,'dsg1',dsg1,ierr)
      call nemsio_getheadvar(gfile,'sgml1',sgml1,ierr)
      call nemsio_getheadvar(gfile,'sg2',sg2,ierr)
      call nemsio_getheadvar(gfile,'dsg2',dsg2,ierr)
      call nemsio_getheadvar(gfile,'sgml2',sgml2,ierr)
!
      if(mype==0)then
        call nemsio_readrecv(gfile,'fis','sfc',1,temp1,iret=ierr)
!        write(0,*)'in init_nemsio,fis=',maxval(temp1),minval(temp1)
      endif
      do j=jms,jme
      do i=ims,ime
        fis(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,fis,1,1,1,1,1)
      call halo_exch(fis,1,2,2)
!-----------------------------------------------------------------------
!
      if(mype==0)then
        call nemsio_readrecv(gfile,'stdh','sfc',1,temp1,iret=ierr)
!        write(0,*)'in init_nemsio,stdh=',maxval(temp1),minval(temp1)
      endif
      do j=jms,jme
      do i=ims,ime
        stdh(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,stdh,1,1,1,1,1)
      call halo_exch(stdh,1,2,2)

!-----------------------------------------------------------------------
!
      if(mype==0)then
        call nemsio_readrecv(gfile,'sm','sfc',1,temp1,iret=ierr)
!        write(0,*)'in init_nemsio,sm=',maxval(temp1),minval(temp1)
      endif
      do j=jms,jme
      do i=ims,ime
        sm(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,sm,1,1,1,1,1)
      call halo_exch(sm,1,2,2)

!-----------------------------------------------------------------------
!
      if(mype==0)then
        call nemsio_readrecv(gfile,'dpres','hybrid sig lev',1,temp1,iret=ierr)
        write(0,*)'in init_nemsio,pd=',maxval(temp1),minval(temp1)
      endif
      do j=jms,jme
      do i=ims,ime
        pd(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,pd,1,1,1,1,1)
      call halo_exch(pd,1,2,2)

!-----------------------------------------------------------------------
!
      call mpi_barrier(mpi_comm_comp,ierr)
      do l=1,lm
        if(mype==0)then
        call nemsio_readrecv(gfile,'ugrd','mid layer',l,temp1,iret=ierr)
!        write(0,*)'in init_nemsio,ugrd=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          u(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,u,1,1,1,lm,l)
      enddo
      call halo_exch(u,lm,2,2)

!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
        call nemsio_readrecv(gfile,'vgrd','mid layer',l,temp1,iret=ierr)
!        write(0,*)'in init_nemsio,vgrd=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          v(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,v,1,1,1,lm,l)
      enddo
      call halo_exch(v,lm,2,2)

!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
        call nemsio_readrecv(gfile,'tmp','mid layer',l,temp1,iret=ierr)
          write(0,*) 'L, T extremes: ', L, minval(temp1),maxval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          t(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,t,1,1,1,lm,l)
      enddo
      call halo_exch(t,lm,2,2)

!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
        call nemsio_readrecv(gfile,'spfh','mid layer',l,temp1,iret=ierr)
!        write(0,*)'in init_nemsio,spfh=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          q(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,q,1,1,1,lm,l)
      enddo
      call halo_exch(q,lm,2,2)

!
      do l=1,lm
      do j=jms,jme
      do i=ims,ime
        water(i,j,l,p_qv)=q(i,j,l)/(1.-q(i,j,l))    ! WRF water array uses mixing ratio for vapor
      enddo
      enddo
      enddo
!
!-----------------------------------------------------------------------
!
      do l=1,lm
        if(mype==0)then
        call nemsio_readrecv(gfile,'clwmr','mid layer',l,temp1,iret=ierr)
!        write(0,*)'in init_nemsio,cldmr=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          cw(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,cw,1,1,1,lm,l)
      enddo
      call halo_exch(cw,lm,2,2)

!
      print*,'*** Read initial conditions in init_nemsio from ',infile
!
!      close(unit=nfcst)
!-----------------------------------------------------------------------
      ntsti=ntsd+1
!
      tend_max=real(ihrend)
      ntstm_max=nint(tend_max*3600./dt)+1
      tend=real(nhours_fcst)
      ntstm=nint(tend*3600./dt)+1
      if(ntstm>ntstm_max.and..not.global)then
        ntstm=min(ntstm,ntstm_max)
      endif
!
      ihr=nint(ntsd*dt/3600.)
!-----------------------------------------------------------------------
      do l=1,lm
        pdsg1(l)=dsg1(l)*pdtop
        psgml1(l)=sgml1(l)*pdtop+pt
!       write(0,*) 'L, pdsg1, psgml1: ', L, pdsg1(L), psgml1(L)
      enddo
!
      do l=1,lm+1
        psg1(l)=sg1(l)*pdtop+pt
      enddo
!-----------------------------------------------------------------------
      do j=jts,jte
        do i=its,ite
          pdo(i,j)=pd(i,j)
        enddo
      enddo
      call halo_exch(pdo,1,2,2)
!
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            up(i,j,l)=u(i,j,l)
            vp(i,j,l)=v(i,j,l)
            tp(i,j,l)=t(i,j,l)
          enddo
        enddo
      enddo
      call halo_exch(tp,lm,up,lm,vp,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            q2(i,j,l)=0.02
            rrw(i,j,l)=0.
            if(i.ge.ide  /2+1- 6.and.i.le.ide  /2+1+ 6.and. &
               j.ge.jde*3/4+1- 6.and.j.le.jde*3/4+1+ 6.) then !global
!               j.ge.jde  /2+1- 6.and.j.le.jde  /2+1+ 6.) then !regional
              rrw(i,j,l)=10.
            endif
            dwdt(i,j,l)=1.
            w(i,j,l)=0.
          enddo
        enddo
      enddo
      call halo_exch(dwdt,lm,2,2)
!
      do j=jts,jte
        do i=its,ite
          pint(i,j,1)=pt
        enddo
      enddo
!
      do l=1,lm
        do j=jts,jte
          do i=its,ite
            pint(i,j,l+1)=pint(i,j,l)+dsg2(l)*pd(i,j)+pdsg1(l)
          enddo
        enddo
      enddo
      call halo_exch(pint,lm+1,2,2)
!
      call halo_exch(q2,lm,rrw,lm,2,2)
      do l=1,lm
        do j=jms,jme
          do i=ims,ime
            sp(i,j,l,indx_q  )=sqrt(max(q  (i,j,l),0.))
            sp(i,j,l,indx_cw )=sqrt(max(cw (i,j,l),0.))
            sp(i,j,l,indx_rrw)=sqrt(max(rrw(i,j,l),0.))
            sp(i,j,l,indx_q2 )=sqrt(max(q2 (i,j,l),0.))
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!---reading surface data------------------------------------------------
!-----------------------------------------------------------------------
!
       if(mype==0)then
        call nemsio_readrecv(gfile,'sice','sfc',1,temp1,iret=ierr)
       endif
       do j=jms,jme
       do i=ims,ime
         sice(i,j)=0.
       enddo
       enddo
       call dstrb(temp1,sice,1,1,1,1,1)
       call halo_exch(sice,1,2,2)
!
!-----------------------------------------------------------------------
!
       call nemsio_close(gfile,iret=ierr)
!
!-----------------------------------------------------------------------
         else                                      ! restart
!-----------------------------------------------------------------------
!
      infile='restart_file_nemsio'
      call nemsio_open(gfile,infile,'read',iret=ierr)
      if(ierr/=0) write(0,*)'ERROR: open file ',trim(infile),' has failed'
!
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: INGEGER SCALARS
!-----------------------------------------------------------------------
      call nemsio_getheadvar(gfile,'FCSTDATE',FCSTDATE,ierr)
      iyear_fcst=FCSTDATE(1)
      imonth_fcst=FCSTDATE(2)
      iday_fcst=FCSTDATE(3)
      ihour_fcst=FCSTDATE(4)
      call nemsio_getheadvar(gfile,'IHRST',ihrst,ierr)
      call nemsio_getheadvar(gfile,'LPT2',lpt2,ierr)
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: INTEGER 1D ARRAYS
!-----------------------------------------------------------------------
!      read(nfcst) idat
      call nemsio_getheadvar(gfile,'IDAT',idat,ierr)
!-----------------------------------------------------------------------
!***  Print the time information.
!-----------------------------------------------------------------------
!
      if(mype==0)then
        write(0,*)'**** read in dynamics ***************'
        write(0,*)' Restart year =',iyear_fcst
        write(0,*)' Restart month=',imonth_fcst
        write(0,*)' Restart day  =',iday_fcst
        write(0,*)' Restart hour =',ihour_fcst
        write(0,*)' Original start year =',idat(3)
        write(0,*)' Original start month=',idat(2)
        write(0,*)' Original start day  =',idat(1)
        write(0,*)' Original start hour =',ihrst
        write(0,*)' Timestep   =',dt
        write(0,*)' Steps/hour =',3600./dt
        write(0,*)'*************************************'
      endif
!
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: INTEGER SCALARS
!-----------------------------------------------------------------------
      call nemsio_getheadvar(gfile,'NSOIL',nsoil,ierr)
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL SCALARS
!-----------------------------------------------------------------------
      call nemsio_getheadvar(gfile,'PDTOP',pdtop,ierr)
!-----------------------------------------------------------------------
      call nemsio_getheadvar(gfile,'PT',pt,ierr)
!      write(0,*)' in rst,pdtop=',pdtop,'pt=',pt,'nsoil=',nsoil,'idat=', &
!       idat,'ntsd=',ntsd,'fcstdate=',fcstdate,'ihrst=',ihrst,'lpt2=',lpt2
!      read(nfcst) pt
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 1D ARRAYS
!-----------------------------------------------------------------------
      call nemsio_getheadvar(gfile,'SG1',sg1,ierr)
      call nemsio_getheadvar(gfile,'SG2',sg2,ierr)
      call nemsio_getheadvar(gfile,'DSG1',dsg1,ierr)
      call nemsio_getheadvar(gfile,'DSG2',dsg2,ierr)
      call nemsio_getheadvar(gfile,'SGML1',sgml1,ierr)
      call nemsio_getheadvar(gfile,'SGML2',sgml2,ierr)
      call nemsio_getheadvar(gfile,'SGM',sgm,ierr)
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: LOGICAL
!-----------------------------------------------------------------------
      call nemsio_getheadvar(gfile,'RUN',run,ierr)
      if(mype==0)then
        write(0,*) 'run, idat,ntsd: ', run, idat, ntsd
      endif
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 2D ARRAYS
!-----------------------------------------------------------------------
      if(mype==0)then
        call nemsio_readrecv(gfile,'fis','sfc',1,temp1,iret=ierr)
        write(0,*)'in init restart,fis=',maxval(temp1),minval(temp1)
      endif
      do j=jms,jme
      do i=ims,ime
        fis(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,fis,1,1,1,1,1)
      call halo_exch(fis,1,2,2)
!-----------------------------------------------------------------------
      if(mype==0)then
        call nemsio_readrecv(gfile,'dpres','hybrid sig lev',1,temp1,iret=ierr)
        write(0,*)'in init restart,pd=',maxval(temp1),minval(temp1)
      endif
      do j=jms,jme
      do i=ims,ime
        pd(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,pd,1,1,1,1,1)
      call halo_exch(pd,1,2,2)
!-----------------------------------------------------------------------
      if(mype==0)then
        call nemsio_readrecv(gfile,'pdo','sfc',1,temp1,iret=ierr)
        write(0,*)'in init restart,pdo=',maxval(temp1),minval(temp1)
      endif
      do j=jms,jme
      do i=ims,ime
        pdo(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,pdo,1,1,1,1,1)
      call halo_exch(pdo,1,2,2)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 3D ARRAYS (only DYN)
!-----------------------------------------------------------------------
!     call mpi_barrier(mpi_comm_comp,ierr)
      do l=1,lm
        if(mype==0)then
        call nemsio_readrecv(gfile,'vvel','mid layer',l,temp1,iret=ierr)
        write(0,*)'in init restart,vvel=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          w(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,w,1,1,1,lm,l)
      enddo
      call halo_exch(w,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'dwdt','mid layer',l,temp1,iret=ierr)
        endif
        do j=jms,jme
        do i=ims,ime
          dwdt(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,dwdt,1,1,1,lm,l)
      enddo
      call halo_exch(dwdt,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm+1
        if(mype==0)then
          call nemsio_readrecv(gfile,'pres','layer',l,temp1,iret=ierr)
        write(0,*)'in init restart,pres layer=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          pint(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,pint,1,1,1,lm+1,l)
      enddo
      call halo_exch(pint,lm+1,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'omgalf','mid layer',l,temp1,iret=ierr)
        endif
        do j=jms,jme
        do i=ims,ime
          omgalf(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,omgalf,1,1,1,lm,l)
      enddo
      call halo_exch(omgalf,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'rrw','mid layer',l,temp1,iret=ierr)
          write(0,*)'in init restart,rrw layer=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          rrw(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,rrw,1,1,1,lm,l)
      enddo
      call halo_exch(rrw,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'div','mid layer',l,temp1,iret=ierr)
        endif
        do j=jms,jme
        do i=ims,ime
          div(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,div,1,1,1,lm,l)
      enddo
      call halo_exch(div,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'rtop','mid layer',l,temp1,iret=ierr)
        endif
        do j=jms,jme
        do i=ims,ime
          rtop(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,rtop,1,1,1,lm,l)
      enddo
      call halo_exch(rtop,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'tcu','mid layer',l,temp1,iret=ierr)
        endif
        do j=jms,jme
        do i=ims,ime
          tcu(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,tcu,1,1,1,lm,l)
      enddo
      call halo_exch(tcu,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'tcv','mid layer',l,temp1,iret=ierr)
        endif
        do j=jms,jme
        do i=ims,ime
          tcv(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,tcv,1,1,1,lm,l)
      enddo
      call halo_exch(tcv,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'tct','mid layer',l,temp1,iret=ierr)
          write(0,*)'min,max tct=',minval(tct),maxval(tct)
        endif
        do j=jms,jme
        do i=ims,ime
          tct(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,tct,1,1,1,lm,l)
      enddo
      call halo_exch(tct,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'tp','mid layer',l,temp1,iret=ierr)
        write(0,*)'in init restart,tp =',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          tp(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,tp,1,1,1,lm,l)
      enddo
      call halo_exch(tp,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'up','mid layer',l,temp1,iret=ierr)
        write(0,*)'in init restart,up =',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          up(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,up,1,1,1,lm,l)
      enddo
      call halo_exch(up,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'vp','mid layer',l,temp1,iret=ierr)
        write(0,*)'in init restart,vp =',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          vp(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,vp,1,1,1,lm,l)
      enddo
      call halo_exch(vp,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'e2','mid layer',l,temp1,iret=ierr)
        write(0,*)'in init restart,e2 =',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          e2(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,e2,1,1,1,lm,l)
      enddo
      call halo_exch(e2,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'z','mid layer',l,temp1,iret=ierr)
        write(0,*)'in init restart,z =',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          z(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,z,1,1,1,lm,l)
      enddo
      call halo_exch(z,lm,2,2)
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 2D ARRAYS (contd.)
!-----------------------------------------------------------------------
      if(mype==0)then
        call nemsio_readrecv(gfile,'sice','sfc',1,temp1,iret=ierr)
        write(0,*)'in init restart,sice =',maxval(temp1),minval(temp1)
      endif
      do j=jms,jme
      do i=ims,ime
        sice(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,sice,1,1,1,1,1)
      call halo_exch(sice,1,2,2)
!-----------------------------------------------------------------------
      if(mype==0)then
        call nemsio_readrecv(gfile,'sm','sfc',1,temp1,iret=ierr)
        write(0,*)'in init restart,sm =',maxval(temp1),minval(temp1)
      endif
      do j=jms,jme
      do i=ims,ime
        sm(i,j)=0.
      enddo
      enddo
      call dstrb(temp1,sm,1,1,1,1,1)
      call halo_exch(sm,1,2,2)
!-----------------------------------------------------------------------
!   READ FROM RESTART FILE: REAL 3D ARRAYS
!-----------------------------------------------------------------------
      call mpi_barrier(mpi_comm_comp,ierr)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'clwmr','mid layer',l,temp1,iret=ierr)
!          write(0,*)'in init restart,clwmr =',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          cw(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,cw,1,1,1,lm,l)
      enddo
!      write(0,*)'in init restart after1,clwmr =',maxval(cw),minval(cw)
      call halo_exch(cw,lm,2,2)
!      write(0,*)'in init restart after2,clwmr =',maxval(cw),minval(cw)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'spfh','mid layer',l,temp1,iret=ierr)
        write(0,*)'in init restart,spfh =',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          q(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,q,1,1,1,lm,l)
      enddo
      call halo_exch(q,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'q2','mid layer',l,temp1,iret=ierr)
        endif
        do j=jms,jme
        do i=ims,ime
          q2(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,q2,1,1,1,lm,l)
      enddo
      call halo_exch(q2,lm,2,2)
!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'tmp','mid layer',l,temp1,iret=ierr)
!        write(0,*)'in init restart,tmp =',maxval(temp1),minval(temp1)
          write(0,*) 'L, T extremes: ', L, minval(temp1),maxval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          t(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,t,1,1,1,lm,l)
      enddo
      call halo_exch(t,lm,2,2)

!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'ugrd','mid layer',l,temp1,iret=ierr)
        endif
        do j=jms,jme
        do i=ims,ime
          u(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,u,1,1,1,lm,l)
      enddo
      call halo_exch(u,lm,2,2)

!-----------------------------------------------------------------------
      do l=1,lm
        if(mype==0)then
          call nemsio_readrecv(gfile,'vgrd','mid layer',l,temp1,iret=ierr)
        endif
        do j=jms,jme
        do i=ims,ime
          v(i,j,l)=0.
        enddo
        enddo
        call dstrb(temp1,v,1,1,1,lm,l)
      enddo
      call halo_exch(v,lm,2,2)

!-----------------------------------------------------------------------
      do n=indx_rrw+1,num_tracers_total
      do l=1,lm
        if(mype==0)then
          write(tn,'(I2.2)')n
          call nemsio_readrecv(gfile,'tracers_'//tn,'mid layer',l,temp1,iret=ierr)
!          write(0,*)'get tracer, name=','tracers_'//tn,maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          tracers(i,j,l,n)=0.
        enddo
        enddo
        call dstrb(temp1,tracers(:,:,:,n),1,1,1,lm,l)
      enddo
      enddo
      call halo_exch(tracers,lm,num_tracers_total,1,2,2)
!
      do n=1,num_water
      do l=1,lm
        do j=jms,jme
        do i=ims,ime
          water(i,j,l,n)=tracers(i,j,l,n+num_tracers_total-num_water)
        enddo
        enddo
      enddo
      enddo

!
!-----------------------------------------------------------------------
!
      call nemsio_close(gfile,iret=ierr)
!
!-----------------------------------------------------------------------
      ntsti=ntsd+1
!
      tend_max=real(ihrend)
      ntstm_max=nint(tend_max*3600./dt)+1
      tend=real(nhours_fcst)
      ntstm=nint(tend*3600./dt)+1
      if(.not.global)then
         write(0,*)' Max runtime is ',tend_max,' hours'
      endif
      if(mype==0)then
        write(0,*)' Requested runtime is ',tend,' hours'
        write(0,*)' NTSTM=',ntstm
      endif
      if(ntstm>ntstm_max.and..not.global)then
        if(mype==0)then
          write(0,*)' Requested fcst length exceeds maximum'
          write(0,*)' Resetting to maximum'
        endif
        ntstm=min(ntstm,ntstm_max)
      endif
!
      ihr=nint(ntsd*dt/3600.)
!-----------------------------------------------------------------------
      do l=1,lm
        pdsg1(l)=dsg1(l)*pdtop
        psgml1(l)=sgml1(l)*pdtop+pt
      enddo
!
      do l=1,lm+1
        psg1(l)=sg1(l)*pdtop+pt
      enddo
!-----------------------------------------------------------------------
      do l=1,lm
        do j=jms,jme
          do i=ims,ime
            sp(i,j,l,indx_q  )=sqrt(max(q  (i,j,l),0.))
            sp(i,j,l,indx_cw )=sqrt(max(cw (i,j,l),0.))
            sp(i,j,l,indx_rrw)=sqrt(max(rrw(i,j,l),0.))
            sp(i,j,l,indx_q2 )=sqrt(max(q2 (i,j,l),0.))
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
         endif                                     ! cold start /restart
!-----------------------------------------------------------------------
!
      deallocate(temp1)
      deallocate(stdh)
      if(mype==0)then
        write(0,*)' EXIT SUBROUTINE INIT pt=',pt
      endif
!-----------------------------------------------------------------------
!
                        end subroutine dynamics_read_nemsio
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
                       end module module_DYNAMICS_INIT_READ
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------


