!-----------------------------------------------------------------------
                        module module_INIT_READ_BIN
!-----------------------------------------------------------------------
use module_include
use module_dm_parallel,only : ids,ide,jds,jde &
                             ,ims,ime,jms,jme &
                             ,its,ite,jts,jte &
                             ,its_h2,ite_h2,jts_h2,jte_h2 &
                             ,lm &
                             ,mype_share,npes,num_pts_max &
                             ,mpi_comm_comp &
                             ,dstrb,idstrb
use module_exchange
use module_constants
use module_solver_internal_state,only: solver_internal_state
use module_microphysics_nmm
!-----------------------------------------------------------------------
!
      implicit none
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
      PRIVATE
      PUBLIC :: read_binary,physics_read_gwd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
       contains
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine read_binary &
      (INT_STATE,global &
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
      ,NSOIL,rc )
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
,nhours_fcst &
,nsoil

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

integer(kind=kint),intent(out) :: &
 rc
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

real(kind=kfpt),dimension(NSOIL) :: &
 soil1din

real(kind=kfpt),dimension(:),allocatable :: &
 all_bc_data

real(kind=kfpt),allocatable,dimension(:,:) :: &
 temp1

integer(kind=kfpt),allocatable,dimension(:,:) :: &
 itemp

real(kind=kfpt),allocatable,dimension(:,:,:) :: &
 tempsoil

logical(kind=klog) :: opened

integer(kind=kint) :: &
 ierr &
,mype

character(64):: &
 infile

real(kind=kfpt),allocatable,dimension(:,:):: &      !im,jm
 stdh                        ! standard deviation of topography height

integer(kind=kint):: &
 ihrend &                    ! maximum forecast length, hours
,ntsd &
,ntstm_max &
,nfcst

      TYPE(SOLVER_INTERNAL_STATE),POINTER :: INT_STATE                    !<-- The physics internal state

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
      DO L=1,LM+1
        PSG1(L)=SG1(L)*PDTOP+PT
      ENDDO
      DO L=1,LM
        PDSG1(L)=DSG1(L)*PDTOP
        PSGML1(L)=SGML1(L)*PDTOP+PT
      ENDDO
!
!-----------------------------------------------------------------------
!***  Before moving on, transfer values to the internal state.
!-----------------------------------------------------------------------
!
      int_state%PT=PT
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
!-----------------------------------------------------------------------
!***  Proceed with getting fields from input file.
!***  NOTE: Two records were already read at the top of this routine.
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
          do j=jms,jme
          do i=ims,ime
            o3(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,o3,1,1,1,lm,l)
        enddo
        call halo_exch(o3,lm,2,2)
!
!-----------------------------------------------------------------------
!
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
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst)temp1  ! ALBASE
        endif
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst)temp1  ! EPSR
        endif
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst)temp1  ! MXSNAL
        endif
      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst)temp1  ! TSKIN
        endif
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst)temp1  ! SST
        endif
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst)temp1  ! SNO
        endif
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst)temp1  ! SI
        endif
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst)temp1  ! SICE
        endif
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst) temp1  ! TG
        endif
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst) temp1  ! CMC
        endif
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst) temp1  ! SR
        endif
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst) temp1  ! USTAR
        endif
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst) temp1  ! Z0
        endif
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!
        if(mype==0)then
          read(nfcst) temp1  ! Z0BASE
        endif
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1)
      CALL HALO_EXCH(int_state%Z0BASE,1,3,3)
!
      ALLOCATE(TEMPSOIL(1:NSOIL,IDS:IDE,JDS:JDE),STAT=I)
!
        if(mype==0)then
          read(nfcst) TEMPSOIL  ! STC
        endif
      CALL DSTRB(TEMPSOIL(1,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,1),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(2,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,2),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(3,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,3),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(4,IDS:IDE,JDS:JDE),int_state%STC(IMS:IME,JMS:JME,4),1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst) TEMPSOIL  ! SMC
        endif
      CALL DSTRB(TEMPSOIL(1,:,:),int_state%SMC(:,:,1),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(2,:,:),int_state%SMC(:,:,2),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(3,:,:),int_state%SMC(:,:,3),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(4,:,:),int_state%SMC(:,:,4),1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst) TEMPSOIL  ! SH2O
        endif
      CALL DSTRB(TEMPSOIL(1,:,:),int_state%SH2O(:,:,1),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(2,:,:),int_state%SH2O(:,:,2),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(3,:,:),int_state%SH2O(:,:,3),1,1,1,1,1)
      CALL DSTRB(TEMPSOIL(4,:,:),int_state%SH2O(:,:,4),1,1,1,1,1)
!
      DEALLOCATE(TEMPSOIL)
      ALLOCATE(ITEMP(IDS:IDE,JDS:JDE),STAT=I)
!
        if(mype==0)then
          read(nfcst) ITEMP  ! ISLTYP
        endif
      CALL IDSTRB(ITEMP,int_state%ISLTYP)
!
        if(mype==0)then
          read(nfcst) ITEMP  ! IVGTYP
        endif
      CALL IDSTRB(ITEMP,int_state%IVGTYP)
!
      DEALLOCATE(ITEMP)
!
        if(mype==0)then
          read(nfcst) temp1  ! VEGFRC
        endif
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1)
!
        if(mype==0)then
          read(nfcst) SOIL1DIN  ! DZSOIL
        endif
!
        if(mype==0)then
          read(nfcst) SOIL1DIN  ! SLDPTH
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
                   ,' in CORE_READ_BINARY'
          rc = ierr
          return
        endif
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
        read(nfcst) iyear_fcst
        read(nfcst) imonth_fcst
        read(nfcst) iday_fcst
        read(nfcst) ihour_fcst
        read(nfcst) !iminute_fcst
        read(nfcst) ! second_fcst
        read(nfcst) ! ntsd
        read(nfcst) ! im
        read(nfcst) ! jm
        read(nfcst) ! lm
        read(nfcst) ihrst
        read(nfcst) i_parent_start
        read(nfcst) j_parent_start
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
          write(0,*)'**** read in core *******************'
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
        read(nfcst) ! nsoil
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
        read(nfcst) pt
        read(nfcst) tlm0d
        read(nfcst) tph0d
        read(nfcst) ! tstart
        read(nfcst) dphd
        read(nfcst) dlmd
        read(nfcst) sbd
        read(nfcst) wbd
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
        read(nfcst) int_state%sldpth
        read(nfcst) int_state%mp_restart_state
        read(nfcst) int_state%tbpvs_state
        read(nfcst) int_state%tbpvs0_state
!
!-----------------------------------------------------------------------
!***  Read in the full-domain 1-D datastring of boundary winds.
!***  Each task isolates its own piece of that data.
!-----------------------------------------------------------------------
!
        length=2*2*2*lnsv*lm*((ide-ids)+(jde-jds))
        allocate(all_bc_data(1:length))
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
!***  Before moving on, transfer values to the internal state.
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

!-----------------------------------------------------------------------
!***  Read from restart file: Integer 2D arrays
!-----------------------------------------------------------------------
!
      ALLOCATE(ITEMP(IDS:IDE,JDS:JDE),STAT=I)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%ISLTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%IVGTYP)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRCV)
!
      IF(MYPE==0)THEN
        READ(NFCST) ITEMP
      ENDIF
      CALL IDSTRB(ITEMP,int_state%NCFRST)
!
      DEALLOCATE(ITEMP)
!
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
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays (contd.)
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  ACFRCV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRCV(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACFRCV,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACFRST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACFRST(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACFRST,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACPREC(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACPREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACSNOM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACSNOM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ACSNOW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ACSNOW(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ACSNOW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKHS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKHS_OUT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%AKHS_OUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKMS_OUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%AKMS_OUT(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%AKMS_OUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALBASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%ALBASE(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%ALBASE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALBEDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALBEDO,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ALWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ALWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  BGROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%BGROFF,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACL,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CFRACM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CFRACM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CLDEFI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CLDEFI,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CMC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CMC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CNVBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVBOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CNVTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CNVTOP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CPRATE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CPRATE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CUPPT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPPT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CUPREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CUPREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CZEN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZEN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  CZMEAN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%CZMEAN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  EPSR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%EPSR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  GRNFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%GRNFLX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOTD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOTS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOTS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOPD
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPD,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOPS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOPS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNOW ALBEDO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%MXSNAL,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PBLH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PBLH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  POTEVP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTEVP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PREC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PREC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  PSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%PSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  Q10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Q10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QSH
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSH,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QWBS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  QZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%QZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RADOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RADOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RLWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RLWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RLWTOA,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWINC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWINC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWOUT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWOUT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCEVP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEVP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCEXC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCEXC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCLHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCLHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SFCSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SFCSHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SI
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SI,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SICE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SICE,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SIGT4
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SIGT4,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SM (Seamask)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
!
      DO J=JMS,JME
      DO I=IMS,IME
        int_state%SM(I,J)=0.
      ENDDO
      ENDDO
      CALL DSTRB(TEMP1,int_state%SM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SMSTAV
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTAV,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SMSTOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SMSTOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNO
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNO,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SNOPCX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SNOPCX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SOILTB
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SOILTB,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SSROFF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SSROFF,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SST
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SST,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  SUBSHX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%SUBSHX,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TG
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TG,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TH10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TH10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  THS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%THS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  THZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%THZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TSHLTR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSHLTR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  TWBS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TWBS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  U10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%U10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  USTAR
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%USTAR,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  UZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%UZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  V10
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%V10,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  VEGFRC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%VEGFRC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  VZ0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST) TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%VZ0,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  Z0
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0,1,1,1,1,1)
      CALL HALO_EXCH(int_state%Z0,1,3,3)
!-----------------------------------------------------------------------
!***  TSKIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TSKIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKHS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKHS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AKMS
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AKMS,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HBOT
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HBOT,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  HTOP
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%HTOP,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  RSWTOA
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RSWTOA,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  POTFLX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%POTFLX,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  RMOL
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%RMOL,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  T2
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%T2,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  Z0BASE
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%Z0BASE,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  TLMIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TLMIN,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  TLMAX
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%TLMAX,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  ACUTIM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ACUTIM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  APHTIM
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%APHTIM,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ARDLW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ARDLW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ARDSW
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ARDSW,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  ASRFC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%ASRFC,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AVRAIN
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AVRAIN,1,1,1,1,1)
!-----------------------------------------------------------------------
!***  AVCNVC
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NFCST)TEMP1
      ENDIF
      CALL DSTRB(TEMP1,int_state%AVCNVC,1,1,1,1,1)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays
!-----------------------------------------------------------------------
        call mpi_barrier(mpi_comm_comp,irtn)
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1 ! cldfra
          endif
          call dstrb(temp1,int_state%CLDFRA,1,1,1,lm,l)
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
          do j=jms,jme
          do i=ims,ime
            int_state%RLWTT(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%RLWTT,1,1,1,lm,l)
        enddo
!-----------------------------------------------------------------------
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1 ! rswtt
          endif
          do j=jms,jme
          do i=ims,ime
            int_state%RSWTT(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%RSWTT,1,1,1,lm,l)
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
          do j=jms,jme
          do i=ims,ime
            int_state%TCUCN(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%TCUCN,1,1,1,lm,l)
        enddo
!-----------------------------------------------------------------------
        do l=1,lm
          if(mype==0)then
            read(nfcst)temp1 ! train
          endif
          do j=jms,jme
          do i=ims,ime
            int_state%TRAIN(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,int_state%TRAIN,1,1,1,lm,l)
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
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! XLEN_MIX
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%XLEN_MIX(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%XLEN_MIX,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! F_ICE
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_ICE(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_ICE,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! F_RIMEF
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RIMEF(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_RIMEF,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!
      DO K=1,LM
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1   ! F_RAIN
        ENDIF
!
        DO J=JMS,JME
        DO I=IMS,IME
          int_state%F_RAIN(I,J,K)=0.
        ENDDO
        ENDDO
!
        CALL DSTRB(TEMP1,int_state%F_RAIN,1,1,1,LM,K)
      ENDDO
!-----------------------------------------------------------------------
!***  SH2O, SMC, STC
!-----------------------------------------------------------------------
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
!          write(0,*) 'lev, min, max for SH2O: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%SH2O,1,1,1,NSOIL,K)
!
      ENDDO
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
!          write(0,*) 'lev, min, max for SMC: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%SMC,1,1,1,NSOIL,K)
!
      ENDDO
!
      DO K=1,NSOIL
!
        IF(MYPE==0)THEN
          READ(NFCST)TEMP1
!          write(0,*) 'lev, min, max for STC: ', k,minval(TEMP1),maxval(TEMP1)
        ENDIF
!
        CALL DSTRB(TEMP1,int_state%STC,1,1,1,NSOIL,K)
!
      ENDDO
!
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
      end subroutine read_binary
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      SUBROUTINE PHYSICS_READ_GWD(INFILE,NGWD,INT_STATE                &
                                 ,MYPE,MPI_COMM_COMP                   &
                                 ,IDS,IDE,JDS,JDE,RC)
!----------------------------------------------------------------------
!
!------------------------
!***  Argument variables
!------------------------
!
      INTEGER,INTENT(IN) :: NGWD,MYPE,MPI_COMM_COMP
      INTEGER,INTENT(IN) :: IDS,IDE,JDS,JDE
!
      CHARACTER(LEN=*),INTENT(IN) :: INFILE
!
      TYPE(SOLVER_INTERNAL_STATE),POINTER,INTENT(INOUT) :: INT_STATE     !<-- The physics internal state
!
      INTEGER,INTENT(OUT) :: RC
!
!---------------------
!***  Local variables
!---------------------
!
      INTEGER :: IERR
!
      REAL,DIMENSION(:,:),ALLOCATABLE :: TEMP_GWD
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      RC = 0
!
      ALLOCATE(TEMP_GWD(IDS:IDE,JDS:JDE))      
!
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        OPEN(unit=NGWD,file=INFILE,status='old',form='unformatted'      &
            ,iostat=IERR)
        IF(IERR/=0)THEN
          WRITE(0,*)' Unable to open file ',TRIM(INFILE)                &
                   ,' in PHYSICS_READ_GWD'
          RC = IERR
          RETURN
        ENDIF
      ENDIF
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HSTDV,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HCNVX,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYS,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYSW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HASYNW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENS,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENSW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HLENNW,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HANGL,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HANIS,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HSLOP,1,1,1,1,1)
!-----------------------------------------------------------------------
      IF(MYPE==0)THEN
        READ(NGWD)TEMP_GWD
      ENDIF
!
      CALL DSTRB(TEMP_GWD,int_state%HZMAX,1,1,1,1,1)
!-----------------------------------------------------------------------
!
      IF(MYPE==0)THEN
        CLOSE(NGWD)
      ENDIF
!
      DEALLOCATE(TEMP_GWD)
!-----------------------------------------------------------------------
!
      END SUBROUTINE PHYSICS_READ_GWD 
!
!-----------------------------------------------------------------------
!#######################################################################
!-----------------------------------------------------------------------
!
      end module module_INIT_READ_BIN
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------

