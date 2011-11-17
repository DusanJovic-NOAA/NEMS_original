!-----------------------------------------------------------------------
                        module module_DYNAMICS_INIT_READ_BIN
!-----------------------------------------------------------------------
use esmf_mod
use module_include
use module_dm_parallel,only : ids,ide,jds,jde &
                             ,ims,ime,jms,jme &
                             ,its,ite,jts,jte &
                             ,its_h2,ite_h2,jts_h2,jte_h2 &
                             ,lm &
                             ,mype_share,npes,num_pts_max &
                             ,mpi_comm_comp &
                             ,dstrb
use module_exchange
use module_constants
!-----------------------------------------------------------------------
!
      implicit none
!
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
      ,sbd,wbd &              
      ,dphd,dlmd &            
      ,tph0d,tlm0d &          
      ,fis,sm,sice &
      ,pd,pdo,pint &
      ,u,v,q2,e2 &
      ,t,q,cw &
      ,tp,up,vp &
      ,o3,dwdt,w &
      ,omgalf,div,z &
      ,rtop &
      ,tcu,tcv,tct &
      ,sp,indx_q,indx_cw,indx_o3,indx_q2 &
      ,ntsti,ntstm &
      ,ihr,ihrst,idat &
      ,run,restart &
      ,num_water,water &
      ,num_tracers_total,tracers &
      ,p_qv,p_qc,p_qr &
      ,p_qi,p_qs,p_qg &
      ,dt,nhours_fcst &
      ,lnsv &
      ,ubs,ubn,ubw,ube &
      ,vbs,vbn,vbw,vbe &
      ,my_domain_id &
      ,i_parent_start,j_parent_start &
       )
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
integer(kind=kint),intent(in) :: &
 kse &
,kss &
,lnsv &
,my_domain_id &
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
,indx_o3 &
,indx_q2 &
,nhours_fcst

real(kind=kfpt),intent(in):: &
dt

integer(kind=kint),intent(out) :: &
 ihr &
,ihrst &
,i_parent_start &
,j_parent_start &
,lpt2 &
,ntsti &
,ntstm

integer(kind=kint),dimension(3),intent(out) :: &
 idat

real(kind=kfpt),intent(out) :: &
 dlmd &
,dphd &
,pdtop &
,pt &
,tlm0d &
,tph0d &
,sbd &
,wbd 

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
,o3 &
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

real(kind=kfpt),dimension(ims:ime,1:lnsv,1:lm,1:2),intent(out) :: &
 ubs &
,ubn &
,vbs &
,vbn

real(kind=kfpt),dimension(1:lnsv,jms:jme,1:lm,1:2),intent(out) :: &
 ubw &
,ube &
,vbw &
,vbe

logical(kind=klog),intent(in) :: &
 global &
,restart 

logical(kind=klog),intent(out) :: &
 run
!
!---------------------
!***  Local variables
!---------------------
!
integer(kind=kint) :: &
 i &                         ! index in x direction
,iend &
,irtn &
,j &                         ! index in y direction
,jend &
,k &                         ! index
,kount &
,ks &                        ! tracer index
,l &                         ! index in p direction
,length &
,n

integer(kind=kint) :: &      ! number of soil levels
 nsoil

integer(kind=kint) :: &      ! dimensions from input file
 im &
,jm &
,lmm &
,lnsh

integer(kind=kint) :: &
 iyear_fcst           &
,imonth_fcst          &
,iday_fcst            &
,ihour_fcst

real(kind=kfpt):: &
 tend,tend_max

real(kind=kfpt),dimension(:),allocatable :: &
 all_bc_data

real(kind=kfpt),allocatable,dimension(:,:) :: &
 temp1

logical(kind=klog) :: opened


integer(kind=kint) :: &
 ierr &
,mype &
,rc
character(64):: &
 infile
real(kind=kfpt),allocatable,dimension(:,:):: &      !im,jm
 stdh                        ! standard deviation of topography height
integer(kind=kint):: &
ihrend &                    ! maximum forecast length, hours
,ntsd &
,ntstm_max &
,nfcst

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
      read_blocks: if(.not.restart) then                     ! cold start
!-----------------------------------------------------------------------
!
        write(infile,'(a,i2.2)')'input_domain_',my_domain_id
        open(unit=nfcst,file=infile,status='old',form='unformatted')
!
!-----------------------------------------------------------------------
!
        read (nfcst) run,idat,ihrst,ihrend,ntsd
        read (nfcst) pt,pdtop,lpt2,sgm,sg1,dsg1,sgml1,sg2,dsg2,sgml2
        read (nfcst) i_parent_start,j_parent_start
        read (nfcst) dlmd,dphd,wbd,sbd,tlm0d,tph0d
        read (nfcst) im,jm,lmm,lnsh
!
!-----------------------------------------------------------------------
!***  Print the time & domain information.
!-----------------------------------------------------------------------
!
        if(mype==0)then
          write(0,*) 'run, idat,ntsd: ', run, idat, ntsd
          write(0,*)' Start year =',idat(3)
          write(0,*)' Start month=',idat(2)
          write(0,*)' Start day  =',idat(1)
          write(0,*)' Start hour =',ihrst
          write(0,*)' Timestep   =',dt
          write(0,*)' Steps/hour =',3600./dt
          if(.not.global)write(0,*)' Max fcst hours=',ihrend
          write(0,*) 'nmm_dyn reads of PT, PDTOP: ',PT,PDTOP
          write(0,*) 'nmm_dyn reads of I_PARENT_START, J_PARENT_START: ',i_parent_start,j_parent_start
          write(0,*) 'nmm_dyn reads of TLM0D, TPH0D: ',tlm0d,tph0d
          write(0,*) 'nmm_dyn reads of DLMD, DPHD: ',dlmd,dphd
          write(0,*) 'nmm_dyn reads of WBD, SBD: ',wbd,sbd
          write(0,*) 'nmm_dyn reads of IM, JM, LM, LNSH: ',im,jm,lmm,lnsh
        endif
!
!-----------------------------------------------------------------------
!
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
!-----------------------------------------------------------------------
!
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1   ! O3
          endif
!         do j=jms,jme
!         do i=ims,ime
!           o3(i,j,l)=0.
!         enddo
!         enddo
!         call dstrb(temp1,o3,1,1,1,lm,l)
        enddo
!       call halo_exch(o3,lm,2,2)
!
!       print*,'*** Read initial conditions in init from ',infile
!
!-----------------------------------------------------------------------
        ntsti=ntsd+1
!
        tend_max=real(ihrend)
        ntstm_max=nint(tend_max*3600./dt)+1
        tend=real(nhours_fcst)
        ntstm=nint(tend*3600./dt)+1
        if(.not.global)then
          if(mype==0)then
            write(0,*)' Max runtime is ',tend_max,' hours'
          endif
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
          do j=jms,jme
            do i=ims,ime
              q2(i,j,l)=0.02
              o3(i,j,l)=0.
              if(i.ge.ide  /2+1- 6.and.i.le.ide  /2+1+ 6.and. &
                 j.ge.jde*3/4+1- 6.and.j.le.jde*3/4+1+ 6.) then !global
!                 j.ge.jde  /2+1- 6.and.j.le.jde  /2+1+ 6.) then !regional
                o3(i,j,l)=10.
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
        call halo_exch(q2,lm,o3,lm,2,2)
        do l=1,lm
          do j=jms,jme
            do i=ims,ime
              sp(i,j,l,indx_q )=sqrt(max(q (i,j,l),0.))
              sp(i,j,l,indx_cw)=sqrt(max(cw(i,j,l),0.))
              sp(i,j,l,indx_o3)=sqrt(max(o3(i,j,l),0.))
              sp(i,j,l,indx_q2)=sqrt(max(q2(i,j,l),0.))
            enddo
          enddo
        enddo
!
!-----------------------------------------------------------------------
!---reading surface data------------------------------------------------
!-----------------------------------------------------------------------
!
        if(mype==0)then
          read(nfcst)temp1  ! ALBEDO
        endif
!
        if(mype==0)then
          read(nfcst)temp1  ! ALBASE
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
          read(nfcst) temp1  ! TG
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! CMC
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! SR
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! USTAR
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! Z0
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! Z0BASE
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! STC
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! SMC
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! SH2O
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! ISLTYP
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! IVGTYP
        endif
!
        if(mype==0)then
          read(nfcst) temp1  ! VEGFRC
        endif
!
        if(mype==0)then
          read(nfcst) ! temp1  ! DZSOIL
        endif
!
        if(mype==0)then
          read(nfcst) ! temp1  ! SLDPTH
        endif
!
!       if(mype==0)then               ! here will be 14 orography fields for GWD
!         do n=1,14
!           read(nfcst) temp1
!         enddo
!       endif
!
!-----------------------------------------------------------------------
!
        close(nfcst)
!
!-----------------------------------------------------------------------
      else  read_blocks                         ! restart
!-----------------------------------------------------------------------
!
        write(infile,'(a,i2.2)')'restart_file_',my_domain_id
        open(unit=nfcst,file=infile,status='old',form='unformatted'     &
            ,iostat=ierr)
        if(ierr/=0)then
          write(0,*)' Unable to open ',trim(infile)                     &
                   ,' in DYNAMICS_READ_BINARY'
          write(0,*)' ABORTING!'
          call esmf_finalize(terminationflag=esmf_abort                 &
                            ,rc             =rc)
        endif
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
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
!
        read(nfcst) i_parent_start
        read(nfcst) j_parent_start
        write(0,*)' INIT restart i_parent_start=',i_parent_start,' j_parent_start=',j_parent_start
        read(nfcst) lpt2
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 1D arrays
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
!***  Read from restart file: Integer scalars
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
!***  Read from restart file: Real scalars
!-----------------------------------------------------------------------
        read(nfcst) ! dt
        read(nfcst) ! dyh
        read(nfcst) pdtop
!-----------------------------------------------------------------------
        read(nfcst) pt
        read(nfcst) tlm0d
        read(nfcst) tph0d
        read(nfcst) ! tstart
        read(nfcst) dphd
        read(nfcst) dlmd
        read(nfcst) sbd
        read(nfcst) wbd
        write(0,*)' INIT restart tlm0d=',tlm0d,' tph0d=',tph0d,' dphd=',dphd,' dlmd=',dlmd,' sbd=',sbd,' wbd=',wbd
!-----------------------------------------------------------------------
!***  Read from restart file: Real 1D arrays
!-----------------------------------------------------------------------
        read(nfcst) ! dxh
        read(nfcst) sg1
        read(nfcst) sg2
        read(nfcst) dsg1
        read(nfcst) dsg2
        read(nfcst) sgml1
        read(nfcst) sgml2
        read(nfcst) sgm
        read(nfcst) ! sldpth
        read(nfcst) ! mp_restart_state
        read(nfcst) ! tbpvs_state
        read(nfcst) ! tbpvs0_state
!
!-----------------------------------------------------------------------
!***  Read in the full-domain 1-D datastring of boundary winds.
!***  Each task isolates its own piece of that data.
!-----------------------------------------------------------------------
!
        length=2*2*2*lnsv*lm*((ide-ids)+(jde-jds))
        allocate(all_bc_data(1:length))
        write(0,*)' INIT allocated all_bc_data to length=',length
!
        read(nfcst) all_bc_data
!
!-----------------------------------------------------------------------
!
        kount=0
!
!-----------------------------------------------------------------------
!
        iend=min(ite_h2,ide-1)
        do n=1,2
        do l=1,lm
        do j=1,lnsv
        do i=ids,ide-1
          if(jts==jds.and.i>=its_h2.and.i<=iend)then                         !<-- South boundary tasks extract their BC winds
            ubs(i,j,l,n)=all_bc_data(kount+1)
            vbs(i,j,l,n)=all_bc_data(kount+2)
          endif
          kount=kount+2
        enddo
        enddo
        enddo
        enddo
!
        do n=1,2
        do l=1,lm
        do j=1,lnsv
        do i=ids,ide-1
          if(jte==jde.and.i>=its_h2.and.i<=iend)then                         !<-- North boundary tasks extract their BC winds
            ubn(i,j,l,n)=all_bc_data(kount+1)
            vbn(i,j,l,n)=all_bc_data(kount+2)
          endif
          kount=kount+2
        enddo
        enddo
        enddo
        enddo
!
        jend=min(jte_h2,jde-1)
        do n=1,2
        do l=1,lm
        do j=jds,jde-1
        do i=1,lnsv
          if(its==ids.and.j>=jts_h2.and.j<=jend)then                         !<-- West boundary tasks extract their BC winds
            ubw(i,j,l,n)=all_bc_data(kount+1)
            vbw(i,j,l,n)=all_bc_data(kount+2)
          endif
          kount=kount+2
        enddo
        enddo
        enddo
        enddo
!
        do n=1,2
        do l=1,lm
        do j=jds,jde-1
        do i=1,lnsv
          if(ite==ide.and.j>=jts_h2.and.j<=jend)then                         !<-- West boundary tasks extract their BC winds
            ube(i,j,l,n)=all_bc_data(kount+1)
            vbe(i,j,l,n)=all_bc_data(kount+2)
          endif
          kount=kount+2
        enddo
        enddo
        enddo
        enddo
!
        deallocate(all_bc_data)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Logical
!-----------------------------------------------------------------------
        read(nfcst) ! global
        read(nfcst) run
        if(mype==0)then
          write(0,*) 'run, idat,ntsd: ', run, idat, ntsd
        endif
        read(nfcst) ! adiabatic
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 2D arrays
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
!***  Read from restart file: Real 2D arrays
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
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays (only DYN)
!-----------------------------------------------------------------------
!       call mpi_barrier(mpi_comm_comp,irtn)
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
            o3(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,o3,1,1,1,lm,l)
        enddo
        call halo_exch(o3,lm,2,2)
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
!        write(0,*)'min,max tct=',minval(tct),maxval(tct)
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
!----------- psgdt -----------------------------------------------------
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
        do n=1,indx_o3
          do l=1,lm
            if(mype==0)then
              read(nfcst)temp1
            endif
            do j=jms,jme
            do i=ims,ime
              sp(i,j,l,n)=0.
            enddo
            enddo
            call dstrb(temp1,sp(:,:,:,n),1,1,1,lm,l)
          enddo
        enddo
        call halo_exch(sp,lm,indx_o3,1,2,2)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays (contd.)
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
!-----------------------------------------------------------------------
        if(mype==0)then
          read(nfcst)temp1  ! tlmin
        endif
!-----------------------------------------------------------------------
        if(mype==0)then
          read(nfcst)temp1  ! tlmax
        endif
!-----------------------------------------------------------------------
        if(mype==0)then
          read(nfcst)temp1  ! acutim
        endif
!-----------------------------------------------------------------------
        if(mype==0)then
          read(nfcst)temp1  ! aphtim
        endif
!-----------------------------------------------------------------------
        if(mype==0)then
          read(nfcst)temp1  ! ardlw
        endif
!-----------------------------------------------------------------------
        if(mype==0)then
          read(nfcst)temp1  ! ardsw
        endif
!-----------------------------------------------------------------------
        if(mype==0)then
          read(nfcst)temp1  ! asrfc
        endif
!-----------------------------------------------------------------------
        if(mype==0)then
          read(nfcst)temp1  ! avrain
        endif
!-----------------------------------------------------------------------
        if(mype==0)then
          read(nfcst)temp1  ! avcnvc
        endif
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays
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
            dwdt(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,dwdt,1,1,1,lm,l)
        enddo
        call halo_exch(dwdt,lm,2,2)
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
        do n=indx_o3+1,num_tracers_total                                   !<-- The first 'indx_o3' arrays are unallocated pointers
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
        close(nfcst)
!
!-----------------------------------------------------------------------
        ntsti=ntsd+1
!
        tend_max=real(ihrend)
        ntstm_max=nint(tend_max*3600./dt)+1
        tend=real(nhours_fcst)
        ntstm=nint(tend*3600./dt)+1
        if(.not.global)then
          if(mype==0)then
            write(0,*)' Max runtime is ',tend_max,' hours'
          endif
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
      endif  read_blocks                        ! cold start /restart
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
!
      end module module_DYNAMICS_INIT_READ_BIN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------


