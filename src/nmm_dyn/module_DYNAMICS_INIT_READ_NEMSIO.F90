                        module module_DYNAMICS_INIT_READ_NEMSIO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
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
use module_nemsio
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
                        subroutine dynamics_read_nemsio &
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
      ,dt,nhours_fcst &
      ,lnsv &
      ,ubs,ubn,ubw,ube &
      ,vbs,vbn,vbw,vbe &
      ,my_domain_id &
      ,i_parent_start,j_parent_start &
       )
!
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
,indx_rrw &
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
,n &
,ierr

integer(kind=kint) :: &      ! number of soil levels
 nsoil

integer(kind=kint) :: &
 iyear_fcst           &
,imonth_fcst          &
,iday_fcst            &
,ihour_fcst
 
integer(kind=kint) :: idate(7)
integer(kind=kint) :: fcstdate(7)
character(2)       ::tn
 
real(kind=kfpt):: &
 tend,tend_max

real(kind=kfpt),dimension(:),allocatable :: &
 all_bc_data &
,temp1

logical(kind=klog) :: opened

type(nemsio_gfile) :: gfile

integer(kind=kint) :: &
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

!
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
!***  Initialize nemsio_module
!-----------------------------------------------------------------------
!
      call nemsio_init()
!-----------------------------------------------------------------------
!
      read_blocks: if(.not.restart) then                 ! cold start
!-----------------------------------------------------------------------
!
!!!     infile='main_input_filename_nemsio'
        write(infile,'(a,i2.2,a)')'input_domain_',my_domain_id,'_nemsio'
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
!***  Get run,idat,ihrst,ihrend,ntsd
!-----------------------------------------------------------------------
!
        call nemsio_getheadvar(gfile,'run',run,ierr)
        call nemsio_getheadvar(gfile,'IDAT',idat,ierr)
        call nemsio_getheadvar(gfile,'ihrst',ihrst,ierr)
        call nemsio_getheadvar(gfile,'ihrend',ihrend,ierr)
        call nemsio_getheadvar(gfile,'ntsd',ntsd,ierr)
!
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
!          write(0,*)'in init_nemsio,fis=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          fis(i,j)=0.
        enddo
        enddo
        call dstrb(temp1,fis,1,1,1,1,1)
        call halo_exch(fis,1,2,2)
!
!-----------------------------------------------------------------------
!
        if(mype==0)then
          call nemsio_readrecv(gfile,'stdh','sfc',1,temp1,iret=ierr)
!          write(0,*)'in init_nemsio,stdh=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          stdh(i,j)=0.
        enddo
        enddo
        call dstrb(temp1,stdh,1,1,1,1,1)
        call halo_exch(stdh,1,2,2)
!
!-----------------------------------------------------------------------
!
        if(mype==0)then
          call nemsio_readrecv(gfile,'sm','sfc',1,temp1,iret=ierr)
!          write(0,*)'in init_nemsio,sm=',maxval(temp1),minval(temp1)
        endif
        do j=jms,jme
        do i=ims,ime
          sm(i,j)=0.
        enddo
        enddo
        call dstrb(temp1,sm,1,1,1,1,1)
        call halo_exch(sm,1,2,2)
!
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
!
!-----------------------------------------------------------------------
!
        call mpi_barrier(mpi_comm_comp,ierr)
        do l=1,lm
          if(mype==0)then
          call nemsio_readrecv(gfile,'ugrd','mid layer',l,temp1,iret=ierr)
!          write(0,*)'in init_nemsio,ugrd=',maxval(temp1),minval(temp1)
          endif
          do j=jms,jme
          do i=ims,ime
            u(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,u,1,1,1,lm,l)
        enddo
        call halo_exch(u,lm,2,2)
!
!-----------------------------------------------------------------------
!
        do l=1,lm
          if(mype==0)then
          call nemsio_readrecv(gfile,'vgrd','mid layer',l,temp1,iret=ierr)
!          write(0,*)'in init_nemsio,vgrd=',maxval(temp1),minval(temp1)
          endif
          do j=jms,jme
          do i=ims,ime
            v(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,v,1,1,1,lm,l)
        enddo
        call halo_exch(v,lm,2,2)
!
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
!
!-----------------------------------------------------------------------
!
        do l=1,lm
          if(mype==0)then
          call nemsio_readrecv(gfile,'spfh','mid layer',l,temp1,iret=ierr)
!          write(0,*)'in init_nemsio,spfh=',maxval(temp1),minval(temp1)
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
!           write(0,*)'in init_nemsio,cldmr=',maxval(temp1),minval(temp1)
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
!
!-----------------------------------------------------------------------
        do l=1,lm
          pdsg1(l)=dsg1(l)*pdtop
          psgml1(l)=sgml1(l)*pdtop+pt
!         write(0,*) 'L, pdsg1, psgml1: ', L, pdsg1(L), psgml1(L)
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
!                 j.ge.jde  /2+1- 6.and.j.le.jde  /2+1+ 6.) then !regional
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
!
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
        call nemsio_getheadvar(gfile,'i_par_sta',i_parent_start,ierr)
        call nemsio_getheadvar(gfile,'j_par_sta',j_parent_start,ierr)
        call nemsio_getheadvar(gfile,'dlmd',dlmd,ierr)
        call nemsio_getheadvar(gfile,'dphd',dphd,ierr)
        call nemsio_getheadvar(gfile,'wbd',wbd,ierr)
        call nemsio_getheadvar(gfile,'sbd',sbd,ierr)
        call nemsio_getheadvar(gfile,'tlm0d',tlm0d,ierr)
        call nemsio_getheadvar(gfile,'tph0d',tph0d,ierr)
!
        call mpi_bcast(pt,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(i_parent_start,1,mpi_integer,0,mpi_comm_comp,irtn)
        call mpi_bcast(j_parent_start,1,mpi_integer,0,mpi_comm_comp,irtn)
        call mpi_bcast(dlmd,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(dphd,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(wbd,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(sbd,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(tlm0d,1,mpi_real,0,mpi_comm_comp,irtn)
        call mpi_bcast(tph0d,1,mpi_real,0,mpi_comm_comp,irtn)
!
!-----------------------------------------------------------------------
!
        call nemsio_close(gfile,iret=ierr)
!
!-----------------------------------------------------------------------
      else read_blocks                          ! restart
!-----------------------------------------------------------------------
!
        write(infile,'(a,i2.2,a)')'restart_file_',my_domain_id,'_nemsio'
        call nemsio_open(gfile,infile,'read',iret=ierr)
        if(ierr/=0) write(0,*)'ERROR: open file ',trim(infile),' has failed'
!
!-----------------------------------------------------------------------
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'FCSTDATE',FCSTDATE,ierr)
        iyear_fcst=FCSTDATE(1)
        imonth_fcst=FCSTDATE(2)
        iday_fcst=FCSTDATE(3)
        ihour_fcst=FCSTDATE(4)
        call nemsio_getheadvar(gfile,'IHRST',ihrst,ierr)
        call nemsio_getheadvar(gfile,'i_par_sta',i_parent_start,ierr)
        call nemsio_getheadvar(gfile,'j_par_sta',j_parent_start,ierr)
        write(0,*)' INIT restart i_parent_start=',i_parent_start,' j_parent_start=',j_parent_start
        call nemsio_getheadvar(gfile,'LPT2',lpt2,ierr)
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 1D arrays
!-----------------------------------------------------------------------
!        read(nfcst) idat
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
!***  Read from restart file: Integer scalars
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'NSOIL',nsoil,ierr)
!-----------------------------------------------------------------------
!***  Read from restart file: Real scalars
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'PDTOP',pdtop,ierr)
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'PT',pt,ierr)
        call nemsio_getheadvar(gfile,'TLM0D',tlm0d,ierr)
        call nemsio_getheadvar(gfile,'TPH0D',tph0d,ierr)
        call nemsio_getheadvar(gfile,'DPHD',dphd,ierr)
        call nemsio_getheadvar(gfile,'DLMD',dlmd,ierr)
        call nemsio_getheadvar(gfile,'SBD',sbd,ierr)
        call nemsio_getheadvar(gfile,'WBD',wbd,ierr)
!       write(0,*)' in rst,pdtop=',pdtop,'pt=',pt,'nsoil=',nsoil,'idat=', &
!        idat,'ntsd=',ntsd,'fcstdate=',fcstdate,'ihrst=',ihrst,'lpt2=',lpt2
!       read(nfcst) pt
!-----------------------------------------------------------------------
!***  Read from restart file: Real 1D arrays
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'SG1',sg1,ierr)
        call nemsio_getheadvar(gfile,'SG2',sg2,ierr)
        call nemsio_getheadvar(gfile,'DSG1',dsg1,ierr)
        call nemsio_getheadvar(gfile,'DSG2',dsg2,ierr)
        call nemsio_getheadvar(gfile,'SGML1',sgml1,ierr)
        call nemsio_getheadvar(gfile,'SGML2',sgml2,ierr)
        call nemsio_getheadvar(gfile,'SGM',sgm,ierr)
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
        if(mype==0)then
        call nemsio_getheadvar(gfile,'ALL_BC_DATA',all_bc_data,ierr)
          write(0,*)' nemsio reads all_bc_data ierr=',ierr
        endif
!
        call mpi_bcast(all_bc_data,length,mpi_real,0,mpi_comm_comp,irtn)
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
        call nemsio_getheadvar(gfile,'RUN',run,ierr)
        if(mype==0)then
          write(0,*) 'run, idat,ntsd: ', run, idat, ntsd
        endif
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays
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
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays (only DYN)
!-----------------------------------------------------------------------
!       call mpi_barrier(mpi_comm_comp,ierr)
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
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays (contd.)
!-----------------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays
!-----------------------------------------------------------------------
!
        call mpi_barrier(mpi_comm_comp,ierr)
!-----------------------------------------------------------------------
        do l=1,lm
          if(mype==0)then
            call nemsio_readrecv(gfile,'clwmr','mid layer',l,temp1,iret=ierr)
!            write(0,*)'in init restart,clwmr =',maxval(temp1),minval(temp1)
          endif
          do j=jms,jme
          do i=ims,ime
            cw(i,j,l)=0.
          enddo
          enddo
          call dstrb(temp1,cw,1,1,1,lm,l)
        enddo
!        write(0,*)'in init restart after1,clwmr =',maxval(cw),minval(cw)
        call halo_exch(cw,lm,2,2)
!        write(0,*)'in init restart after2,clwmr =',maxval(cw),minval(cw)
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
!          write(0,*)'in init restart,tmp =',maxval(temp1),minval(temp1)
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
!
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
!
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
!
!-----------------------------------------------------------------------
        do n=indx_rrw+1,num_tracers_total                                    !<-- The first 'indx_rrw' arrays are unallocated pointers
        do l=1,lm
          if(mype==0)then
            write(tn,'(I2.2)')n
            call nemsio_readrecv(gfile,'tracers_'//tn,'mid layer',l,temp1,iret=ierr)
!            write(0,*)'get tracer, name=','tracers_'//tn,maxval(temp1),minval(temp1)
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
      end subroutine dynamics_read_nemsio
!
!-----------------------------------------------------------------------
!
      end module module_DYNAMICS_INIT_READ_NEMSIO
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------


