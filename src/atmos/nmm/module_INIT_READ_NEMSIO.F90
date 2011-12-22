                        module module_INIT_READ_NEMSIO
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
use nemsio_module_mpi
use module_control,only : timef
USE MODULE_SOLVER_INTERNAL_STATE,ONLY: SOLVER_INTERNAL_STATE
!
!-----------------------------------------------------------------------
!
      implicit none
!
      private
!
      public read_nemsio
      public physics_read_input_nemsio
      public physics_read_restt_nemsio
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
                        subroutine read_nemsio &
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

!---------------------
!***  Local variables
!---------------------
!
integer(kind=kint) :: &
 i &                         ! index in x direction
,iend &
,ierr &
,irtn &
,j &                         ! index in y direction
,jend &
,k &                         ! index
,kount &
,ks &                        ! tracer index
,l &                         ! index in p direction
,length &
,n &
,rc

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
,tmp

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

integer nrec,mysize,myrank
integer fldsize,fldst,js,recn
character(16),allocatable :: recname(:), reclevtyp(:)
integer,allocatable       :: reclev(:)
!
real(8) :: stime,etime,stime1
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      stime=timef()
!
      mype=mype_share
!
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
        call nemsio_open(gfile,infile,'read',mpi_comm_comp,iret=ierr)
        if(ierr/=0)then
          write(0,*)'ERROR: Unable to open file ',trim(infile)          &
                   ,' in READ_NEMSIO'
          write(0,*)' ABORTING!'
          call esmf_finalize(terminationflag=esmf_abort                 &
                            ,rc             =rc)
        endif
!!!     write(0,*)'after nemsio_open, t=',timef()-stime
!
        call nemsio_getfilehead(gfile,nrec=nrec,iret=ierr)
!
        allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
        call nemsio_getfilehead(gfile,recname=recname,reclevtyp=reclevtyp,   &
                                reclev=reclev,iret=ierr)
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
        fldsize=(jte-jts+1)*(ite-its+1)
        allocate(tmp((ite-its+1)*(jte-jts+1)*nrec),stat=i)
!
        stime1=timef()
        call nemsio_denseread(gfile,its,ite,jts,jte,tmp,iret=ierr)
        if(ierr/=0) then
          print *,'WRONG: Could not read all the fields in the file!'
        endif
!
!-- fis
        fis=0.
        call getrecn(recname,reclevtyp,reclev,nrec,'fis','sfc',1,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              fis(i,j)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
        call halo_exch(fis,1,2,2)
!        write(0,*)'in init_nemsio,fis=',maxval(fis),minval(fis),'fis(its:its+2,1:3)=', &
!           fis(its:its+2,jts:jts+2),'its=',its,'ite=',ite,'jts=',jts,'jte=',jte, &
!          'tmp(1:3)=',tmp(1:3)

!        if(mype==0)then
!          call nemsio_readrecv(gfile,'fis','sfc',1,temp1,iret=ierr)
!      if(mype==0) write(0,*)'dyn_init_read,after get fis, t=',timef()-stime
!!          write(0,*)'in init_nemsio,fis=',maxval(temp1),minval(temp1)
!        endif
!        do j=jms,jme
!        do i=ims,ime
!          fis(i,j)=0.
!        enddo
!        enddo
!      if(mype==0) write(0,*)'dyn_init_read,bf dstrb fis, t=',timef()-stime
!        call dstrb(temp1,fis,1,1,1,1,1)
!      if(mype==0) write(0,*)'dyn_init_read,after dstrb fis, t=',timef()-stime
!        call halo_exch(fis,1,2,2)
!      if(mype==0) write(0,*)'dyn_init_read,after exch fis, t=',timef()-stime
!
!-----------------------------------------------------------------------
!
!-- sm
        sm=0.
        call getrecn(recname,reclevtyp,reclev,nrec,'sm','sfc',1,recn)
!        write(0,*)'after sm, recn=',recn,'nrec=',nrec 
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              sm(i,j)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
        call halo_exch(sm,1,2,2)
!          call nemsio_readrecv(gfile,'sm','sfc',1,temp1,iret=ierr)
!         write(0,*)'in init_nemsio,sm=',maxval(sm),minval(sm)
!
!-----------------------------------------------------------------------
!
!-- dpres
        pd=0.
        call getrecn(recname,reclevtyp,reclev,nrec,'dpres','hybrid sig lev',1,recn)
!        write(0,*)'after pd, recn=',recn,'nrec=',nrec 
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              pd(i,j)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
        call halo_exch(pd,1,2,2)
!          call nemsio_readrecv(gfile,'dpres','hybrid sig lev',1,temp1,iret=ierr)
!          write(0,*)'in init_nemsio,pd=',maxval(pd),minval(pd)
!
!-----------------------------------------------------------------------
!
!-- ugrd
      u=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'ugrd','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              u(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(u,lm,2,2)
!      if(mype==0) write(0,*)'dyn_init_read,bf get ugrd2 , t=',timef()-stime
!          call nemsio_readrecv(gfile,'ugrd','mid layer',l,temp1,iret=ierr)
!          write(0,*)'in init_nemsio,ugrd=',maxval(tmp1),minval(temp1)
!
!-----------------------------------------------------------------------
!
!--vgrd
        v=0.
        do l=1,lm
         call getrecn(recname,reclevtyp,reclev,nrec,'vgrd','mid layer',l,recn)
         if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              v(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
         endif
        enddo
        call halo_exch(v,lm,2,2)
!
!-----------------------------------------------------------------------
!--tmp
        t=0.
        do l=1,lm
         call getrecn(recname,reclevtyp,reclev,nrec,'tmp','mid layer',l,recn)
         if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              t(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
         endif
        enddo
        call halo_exch(t,lm,2,2)
!
!-----------------------------------------------------------------------
!
!--spfh
        q=0.
        do l=1,lm
         call getrecn(recname,reclevtyp,reclev,nrec,'spfh','mid layer',l,recn)
         if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              q(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
         endif
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
!-- clwmr
        cw=0.
        do l=1,lm
         call getrecn(recname,reclevtyp,reclev,nrec,'clwmr','mid layer',l,recn)
         if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              cw(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
         endif
        enddo
        call halo_exch(cw,lm,2,2)
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
!-- sice
        sice=0.
        call getrecn(recname,reclevtyp,reclev,nrec,'sice','sfc',1,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              sice(i,j)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
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
        call nemsio_open(gfile,infile,'read',mpi_comm_comp,iret=ierr)
        if(ierr/=0)then
          write(0,*)'ERROR: open file ',trim(infile),' has failed'
          write(0,*)' ABORTING!'
          call esmf_finalize(terminationflag=esmf_abort                 &
                            ,rc             =rc)
        endif
!!!     write(0,*)'dyn_init_read,aft open nemsio,',trim(infile),'time=',timef()-stime
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
        if(mype == 0) then
         write(0,*)' INIT restart i_parent_start=',i_parent_start,' j_parent_start=',j_parent_start,'ierr=',ierr
        endif
        call nemsio_getheadvar(gfile,'LPT2',lpt2,ierr)
!-----------------------------------------------------------------------
!***  Read from restart file: Integer 1D arrays
!-----------------------------------------------------------------------
        call nemsio_getheadvar(gfile,'IDAT',idat,ierr)
!-----------------------------------------------------------------------
!***  Print the time information.
!-----------------------------------------------------------------------
!
        if(mype==0)then
          write(0,*)'**** in read_nemsio *****************'
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
!
        call nemsio_getheadvar(gfile,'ALL_BC_DATA',all_bc_data,ierr)
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
!
      fldsize=(jte-jts+1)*(ite-its+1)
!      write(0,*)'its=',its,'ite=',ite,'jts=',jts,'jte=',jte,'fldsize=',fldsize
!
      call nemsio_getfilehead(gfile,nrec=nrec,iret=ierr)
      allocate(recname(nrec),reclevtyp(nrec),reclev(nrec))
      call nemsio_getfilehead(gfile,recname=recname,reclevtyp=reclevtyp,   &
                              reclev=reclev,iret=ierr)
!
      allocate(tmp(fldsize*nrec))
      tmp=0.
      stime=timef()
      call nemsio_denseread(gfile,its,ite,jts,jte,tmp,iret=ierr)
!!!   write(0,*)'aft nemsio_denseread, time=',timef()-stime,'tmp=',maxval(tmp), &
!!!     minval(tmp)
!
!-- fis
      fis=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'fis','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            fis(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(fis,1,2,2)

!          call nemsio_readrecv(gfile,'fis','sfc',1,temp1,iret=ierr)
!      write(0,*)'in init restart dyn,fis=',maxval(fis),minval(fis),'fiss=', &
!           fis(its:its+2,jts:jts+2),'fisend=',fis(ite-2:ite,jte-2:jte)
!-----------------------------------------------------------------------
!--pd
      pd=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'dpres','hybrid sig lev',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            pd(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(pd,1,2,2)
!      write(0,*)'in init restart dyn,pd=',maxval(pd),minval(pd),'pd=', &
!           pd(its:its+2,jts:jts+2),'pde=',pd(ite-2:ite,jte-2:jte)
!-----------------------------------------------------------------------
!-- pdo
      pdo=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'pdo','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            pdo(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(pdo,1,2,2)
!      write(0,*)'in init restart dyn,pdo=',maxval(pdo),minval(pdo),'pdo=', &
!           pdo(its:its+2,jts:jts+2),'pdoe=',pdo(ite-2:ite,jte-2:jte)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays (only DYN)
!-----------------------------------------------------------------------
!w
      w=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'vvel','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              w(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(w,lm,2,2)
!
!-----------------------------------------------------------------------
!-- dwdt
      dwdt=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'dwdt','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              dwdt(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(dwdt,lm,2,2)
!-----------------------------------------------------------------------
!-- pres
      pint=0.
      do l=1,lm+1
        call getrecn(recname,reclevtyp,reclev,nrec,'pres','layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              pint(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(pint,lm+1,2,2)
!-----------------------------------------------------------------------
!-- omgalf
      omgalf=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'omgalf','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              omgalf(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(omgalf,lm,2,2)
!-----------------------------------------------------------------------
!-- o3mr
      o3=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'o3mr','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              o3(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(o3,lm,2,2)
!-----------------------------------------------------------------------
!-- div
      div=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'div','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              div(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(div,lm,2,2)
!-----------------------------------------------------------------------
!-- rtop
      rtop=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'rtop','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              rtop(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(rtop,lm,2,2)
!-----------------------------------------------------------------------
!-- tcu
      tcu=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tcu','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              tcu(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(tcu,lm,2,2)
!-----------------------------------------------------------------------
!-- tcv
      tcv=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tcv','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              tcv(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(tcv,lm,2,2)
!-----------------------------------------------------------------------
!-- tct
      tct=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tct','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              tct(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(tct,lm,2,2)
!-----------------------------------------------------------------------
!--tp
      tp=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tp','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              tp(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(tp,lm,2,2)
!-----------------------------------------------------------------------
!-- up
      up=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'up','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              up(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(up,lm,2,2)
!-----------------------------------------------------------------------
!-- vp
      vp=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'vp','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              vp(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(vp,lm,2,2)
!-----------------------------------------------------------------------
!-- e2
      e2=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'e2','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              e2(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(e2,lm,2,2)
!-----------------------------------------------------------------------
!-- z
      z=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'z','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              z(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(z,lm,2,2)
!-----------------------------------------------------------------------
!-- tracers
      do n=1,indx_o3
        sp(:,:,:,n)=0.
        write(tn,'(I2.2)')n
        do l=1,lm
          call getrecn(recname,reclevtyp,reclev,nrec,'tracers_prev_'//tn,  &
               'mid layer',l,recn)
          if(recn/=0) then
            fldst=(recn-1)*fldsize
            do j=jts,jte
              js=(j-jts)*(ite-its+1)
              do i=its,ite
                sp(i,j,l,n)=tmp(i-its+1+js+fldst)
              enddo
            enddo
          endif
        enddo
      enddo
      call halo_exch(sp,lm,indx_o3,1,2,2)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 2D arrays (contd.)
!-----------------------------------------------------------------------
!
!-- sice
      sice=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'sice','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            sice(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(sice,1,2,2)
!-----------------------------------------------------------------------
!-- sm
      sm=0.
      call getrecn(recname,reclevtyp,reclev,nrec,'sm','sfc',1,recn)
      if(recn/=0) then
        fldst=(recn-1)*fldsize
        do j=jts,jte
          js=(j-jts)*(ite-its+1)
          do i=its,ite
            sm(i,j)=tmp(i-its+1+js+fldst)
          enddo
        enddo
      endif
      call halo_exch(sm,1,2,2)
!
!-----------------------------------------------------------------------
!***  Read from restart file: Real 3D arrays
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!-- cw
      cw=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'clwmr','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              cw(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(cw,lm,2,2)
!      write(0,*)'in init restart after2,clwmr =',maxval(cw),minval(cw)
!-----------------------------------------------------------------------
!-- spfh
      q=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'spfh','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              q(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(q,lm,2,2)
!-----------------------------------------------------------------------
!-- q2
      q2=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'q2','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              q2(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(q2,lm,2,2)
!-----------------------------------------------------------------------
!-- t
      t=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'tmp','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              t(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(t,lm,2,2)
!-----------------------------------------------------------------------
!-- u
      u=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'ugrd','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              u(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(u,lm,2,2)
!-----------------------------------------------------------------------
!-- v
      v=0.
      do l=1,lm
        call getrecn(recname,reclevtyp,reclev,nrec,'vgrd','mid layer',l,recn)
        if(recn/=0) then
          fldst=(recn-1)*fldsize
          do j=jts,jte
            js=(j-jts)*(ite-its+1)
            do i=its,ite
              v(i,j,l)=tmp(i-its+1+js+fldst)
            enddo
          enddo
        endif
      enddo
      call halo_exch(v,lm,2,2)
!-----------------------------------------------------------------------
!-- rest of tracers
      tracers(:,:,:,indx_o3+1:num_tracers_total)=0.
      do n=indx_o3+1,num_tracers_total                                     !<-- The first 'indx_o3' arrays are unallocated pointers
        write(tn,'(I2.2)')n
        do l=1,lm
          call getrecn(recname,reclevtyp,reclev,nrec,'tracers_'//tn,     &
                       'mid layer',l,recn)
          if(recn/=0) then
            fldst=(recn-1)*fldsize
            do j=jts,jte
              js=(j-jts)*(ite-its+1)
              do i=its,ite
                tracers(i,j,l,n)=tmp(i-its+1+js+fldst)
              enddo
            enddo
          endif
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
       call mpi_barrier(mpi_comm_comp,ierr)
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
      etime=timef()
!!!   write(0,*)'restart=',restart,' READ_NEMSIO=',etime-stime
!
      deallocate(tmp)
      deallocate(stdh)
      if(mype==0)then
        write(0,*)' EXIT SUBROUTINE INIT pt=',pt
      endif
!-----------------------------------------------------------------------
!
      end subroutine read_nemsio
!
!-----------------------------------------------------------------------
!
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
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: INT_STATE                 
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
      TYPE(SOLVER_INTERNAL_STATE),POINTER :: INT_STATE
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
          if(trim(recname(i))==trim(fldname).and.                        &
            trim(reclevtyp(i))==trim(fldlevtyp) .and.                    &
            reclev(i)==fldlev) then
            recn=i
            return
          endif
        enddo
!
        if(recn==0) print *,'WARNING: field ',trim(fldname),' ',         &
          trim(fldlevtyp),' ',fldlev,' is not in the nemsio file!'
!
!-----------------------------------------------------------------------
!
      END SUBROUTINE getrecn
!
!-----------------------------------------------------------------------
!
      end module module_INIT_READ_NEMSIO
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------


