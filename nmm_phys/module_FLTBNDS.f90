!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
                        module module_fltbnds
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
use module_include
use module_control,only : im,jm,klog,kint,kfpt &
                         ,s_bdy,n_bdy,w_bdy,e_bdy &
                         ,read_global_sums,write_global_sums
!
use module_dm_parallel,only : ids,ide,jds,jde &
                             ,ims,ime,jms,jme &
                             ,its,ite,jts,jte &
                             ,its_h1,ite_h1,jts_h1,jte_h1 &
                             ,its_h2,ite_h2,jts_h2,jte_h2 &
                             ,its_b1,ite_b2,jts_b2,jte_b2 &
                             ,local_istart,local_iend &
                             ,local_jstart,local_jend &
                             ,gather_layers,scatter_layers &
                             ,mpi_comm_comp,mpi_intra &
                             ,mype_share
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 pi=3.141592653589793238462643383279502884197169399375105820 &
,pih=pi/2. &
,tpi=pi+pi &
,dtr=pi/180.
!
integer(kind=kint) :: &
 ipe_start_north &
,ipe_start_south &
,ipe_end_north &
,ipe_end_south &
,iunit_pole_sums &
,jh_start_fft_north &
,jh_end_fft_north &
,jh_start_fft_south &
,jh_end_fft_south &
,jv_start_fft_north &
,jv_end_fft_north &
,jv_start_fft_south &
,jv_end_fft_south &
,lm_fft &   ! Max number of model layers per task for FFTs
,msize_dummy_fft 
!
integer(kind=kint),dimension(mpi_status_size),private :: &
 istatw
!
integer(kind=kint),private :: &
 mype
integer(kind=kint),allocatable,dimension(:) :: &
 k1_fft &
,k2_fft
!
real(kind=kfpt),allocatable,dimension(:,:,:) :: &
 hn &
,un &
,vn &
,wn
!
logical(kind=klog) :: &
 fft_north &
,fft_south 
!
logical(kind=klog),allocatable,dimension(:) :: &
 my_domain_has_fft_lats_h &
,my_domain_has_fft_lats_v
!
!-----------------------------------------------------------------------
      contains
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine prefft &
(dlmd,dphd,sbd,lm &
,khfilt,kvfilt &
,hfilt,vfilt &
,craux1,craux2,craux3,rcaux1,rcaux2,rcaux3 &
,inpes,jnpes)
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------
real(kind=kfpt),parameter :: &
 cxnc=0.0 &
,rwind=1./3.0 &
,cfilt=1.

integer(kind=kint),intent(in) :: &
 inpes &
,jnpes &
,lm

integer(kind=kint),dimension(jds:jde),intent(out) :: &
 khfilt &
,kvfilt

real(kind=kfpt),intent(in) :: &
 dlmd &
,dphd &
,sbd
 
real(kind=kfpt),dimension(ids:(ide-3)/2+1,jds:jde),intent(out) :: &
 hfilt &
,vfilt
 
real(kind=kdbl),dimension(1:25000),intent(out) :: &
 craux1 &
,rcaux1
 
real(kind=kdbl),dimension(1:20000),intent(out) :: &
 craux2 &
,rcaux2
 
real(kind=kdbl),dimension(1:1),intent(out) :: &
 craux3 &
,rcaux3

!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------

logical(kind=klog) :: &
 first

integer(kind=kint) :: &
 icycle &
,icyclc &
,imax &
,istat &
,j &
,jmax &
,jmax_north &
,jmax_south &
,k &
,k2 &
,ks &
,l_remain &
,lyr_frac_north &
,lyr_frac_south &
,n &
,nnew &
,npe &
,npes &
,npes_north &
,npes_south &
,nsmud

real(kind=kfpt) :: &
 cpf &
,cph &
,cx &
,cxn &
,dlm &
,dph &
,flt &
,rcph &   
,rcycle &
,sb &
,sub_j &
,sxl &
,tph &
,x &
,xl

real(kind=kfpt),dimension(ids:ide-3) :: &
 rarg
 
complex(kind=kfpt),dimension(ids:(ide-3)/2+1) :: &
 carg

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      icycle=ide-3       ! icycle is preferrably an even number
      icyclc=icycle/2+1
      rcycle=1./icycle
!                                                      
      dlm=dlmd*dtr                                                      
      dph=dphd*dtr
      sb=sbd*dtr
!
      cpf=dph/dlm
!-----------------------------------------------------------------------
!***  Maximum number of layers used per task for FFTs.                 
!***  Since the FFT regions cap the poles, there are 2*LM layers
!***  in play.
!-----------------------------------------------------------------------
!
      npes=inpes*jnpes     
!
!-----------------------------------------------------------------------
!-------------filters on h and v rows-----------------------------------
!-----------------------------------------------------------------------
!
      jh_start_fft_south=-10
      jh_end_fft_south=-10
      jh_start_fft_north=-10
      jh_end_fft_north=-10
!
      sub_j=real(jds+1)
!
!-----------------------------------------------------------------------
!
      h_filters: do j=jds+2,jde-2
!
!-----------------------------------------------------------------------
        tph=sb+(j-sub_j)*dph
        cph=cos(tph)
        rcph=(cph/cpf)
        ks=icycle
        nsmud=0
        khfilt(j)=icyclc+1
!
!-----------------------------------------------------------------------
!***  Perform Fourier filtering where cos(phi) < dph/dlm .
!-----------------------------------------------------------------------
!
        if(rcph.lt.cfilt) then
!
          if(tph<0.)then
            if(jh_start_fft_south<-5)jh_start_fft_south=j
            jh_end_fft_south=j
          else
            if(jh_start_fft_north<-5)jh_start_fft_north=j
            jh_end_fft_north=j
          endif
!
          khfilt(j)=icyclc
          hfilt(1,j)=1.
          first=.true.
!
          do k=2,icyclc
            x=(k-1)*dlm*0.5
            xl=min(x/rcph*cfilt,pih)
            sxl=sin(xl)+sin(2.*xl)*0.5*rwind
!
            if(rcph/sxl.gt.cfilt) then
              flt=1.
              ks=k
            else
!
!--filter definition----------------------------------------------------
              cx=cos(x)
              if(abs(cx).lt.1.e-7) first=.false.
              if(first) then
                nnew=(log(sxl/(xl*(1.+rwind)))/log(cx)+0.5)
                if(nnew.gt.nsmud) nsmud=nnew
                if(xl.eq.pih) first=.false.
              endif
!-----------------------------------------------------------------------
              cxn=cx**nsmud
!
              if(cxn.lt.cxnc) then
                khfilt(j)=k-1
!
                do k2=k,icyclc
                  hfilt(k2,j)=0.
                enddo
!
!!!!!           if(mype.eq.0) print*,j,khfilt(j)
                cycle h_filters
              endif
!
              flt=min(cxn,1.)
            endif
!
            hfilt(k,j)=flt
!
!!!!   if(mype.eq.0.and.j.lt.10) &
!!!!   print*,'hfilt ',j,k,nsmud,hfilt(k,j)
!
          enddo
!
        endif
!-----------------------------------------------------------------------
!
      enddo h_filters
!
!-----------------------------------------------------------------------
!
      if(ks.gt.khfilt(j)) ks=0
!
      khfilt(jds)=khfilt(jds+2)
      khfilt(jds+1)=1
      khfilt(jde-1)=1
      khfilt(jde  )=khfilt(jde-2)
!
      do k=1,icycle/2+1
        hfilt(k,jds)=hfilt(k,jds+2)
        hfilt(k,jds+1)=0.
        hfilt(k,jde-1)=0.
        hfilt(k,jde  )=hfilt(k,jde-2)
      enddo
!
      hfilt(1,jds+1)=1.
      hfilt(1,jde-1)=1.
!-----------------------------------------------------------------------
!
      jv_start_fft_south=-10
      jv_end_fft_south=-10
      jv_start_fft_north=-10
      jv_end_fft_north=-10
!
      sub_j=jds+0.5
!
!-----------------------------------------------------------------------
!
      v_filters: do j=jds+1,jde-2
!
!-----------------------------------------------------------------------
        tph=sb+(j-sub_j)*dph
        cph=cos(tph)
        rcph=(cph/cpf)
        ks=icycle
        nsmud=0
        kvfilt(j)=icyclc+1
!
!-----------------------------------------------------------------------
!***  Perform Fourier filtering where cos(phi) < dph/dlm .
!-----------------------------------------------------------------------
        if(rcph.lt.cfilt) then
!
          if(tph<0.)then
            if(jv_start_fft_south<-5)jv_start_fft_south=j
            jv_end_fft_south=j
          else
            if(jv_start_fft_north<-5)jv_start_fft_north=j
            jv_end_fft_north=j
          endif
!
          kvfilt(j)=icyclc
          vfilt(1,j)=1.
          first=.true.
!
          do k=2,icyclc
            x=(k-1)*dlm*0.5
            xl=min(x/rcph*cfilt,pih)
            sxl=sin(xl)+sin(2.*xl)*0.5*rwind
!
            if(rcph/sxl.gt.cfilt) then
              flt=1.
              ks=k
            else
!
!--filter definition----------------------------------------------------
              cx=cos(x)
              if(abs(cx).lt.1.e-7) first=.false.
              if(first) then
                nnew=(log(sxl/(xl*(1.+rwind)))/log(cx)+0.5)
                if(nnew.gt.nsmud) nsmud=nnew
                if(xl.eq.pih) first=.false.
              endif
!-----------------------------------------------------------------------
              cxn=cx**nsmud
!
              if(cxn.lt.cxnc) then
                kvfilt(j)=k-1
!
                do k2=k,icyclc
                  vfilt(k2,j)=0.
                enddo
!
!!!!            if(mype.eq.0.and.j.lt.10) &
!!!!            print*,'vfilt ',j,k,nsmud,vfilt(k,j),vfilt(k+1,j)
                cycle v_filters
              endif
!
              flt=min(cxn,1.)
            endif
!
            vfilt(k  ,j)=flt
!
!!!!   if(mype.eq.0.and.j.lt.10) &
!!!!   print*,'vfilt ',j,k,nsmud,vfilt(k,j)
!
          enddo
!
        endif
!
!-----------------------------------------------------------------------
!
      enddo v_filters
!
!-----------------------------------------------------------------------
!
      if(ks.gt.kvfilt(j)) ks=0
!
      kvfilt(jds)=kvfilt(jds+1)
      kvfilt(jde-1)=kvfilt(jde-2)
!
      do k=1,icycle/2+1
        vfilt(k,jds)=vfilt(k,jds+1)
        vfilt(k,jde-1)=vfilt(k,jde-2)
      enddo
!
!-----------------------------------------------------------------------
!***  Initialize the 1 FFT algorithms
!-----------------------------------------------------------------------
!
      call srcft(1,rarg,icycle,carg,icyclc,icycle,1,+1,1.     &
                ,rcaux1,25000,rcaux2,20000,rcaux3,1)
      call scrft(1,carg,icyclc,rarg,icycle,icycle,1,-1,rcycle &
                ,craux1,25000,craux2,20000,craux3,1)
!
!----------------------------------------------------------------------
!
      if(jh_start_fft_south==jds+2)jh_start_fft_south=jds+1
      if(jh_end_fft_north==jde-2)jh_end_fft_north=jde-1
!
!-----------------------------------------------------------------------
!***  Identify tasks as handling Northern or Southern Hemipshere
!***  model layers for FFTs.
!-----------------------------------------------------------------------
!
      fft_south=.false.
      fft_north=.false.
      ipe_start_south=0
      ipe_end_north=npes-1
!
      if(jts<jde/2)then
        fft_south=.true.
      else
        fft_north=.true.
      endif
!
      south: do n=0,npes-1
        if(local_jstart(n)<jde/2)then
          ipe_end_south=n
        else
          exit south
        endif
      enddo south
!
      north: do n=npes-1,0,-1
        if(local_jstart(n)<jde/2)then
          ipe_start_north=n+1
          exit north
        endif
      enddo north
!
!-----------------------------------------------------------------------
!***  Does this task's subdomain contain any FFT latitudes?
!-----------------------------------------------------------------------
!
      allocate(my_domain_has_fft_lats_h(0:npes-1),stat=istat)
      allocate(my_domain_has_fft_lats_v(0:npes-1),stat=istat)
!
      do n=0,npes-1
        my_domain_has_fft_lats_h(n)=.false.
        my_domain_has_fft_lats_v(n)=.false.
      enddo
!
      do n=0,npes-1
!
        if(fft_south.and.local_jstart(n)<=jh_end_fft_south.or. &
           fft_north.and.local_jend(n)>=jh_start_fft_north)then
          my_domain_has_fft_lats_h(n)=.true.
        endif
!
        if(fft_south.and.local_jstart(n)<=jv_end_fft_south.or. &
           fft_north.and.local_jend(n)>=jv_start_fft_north)then
          my_domain_has_fft_lats_v(n)=.true.
        endif
!
      enddo
!
!-----------------------------------------------------------------------
!***  Assign layers to each task.
!***  k1_fft and k2_fft are the first and last model layers in
!***  a task's group of layers over which it will apply FFTs.
!***  Groups of model layers will be assigned from top down
!***  in the southern or northern hemisphere and then divided
!***  if there are more than 2*LM MPI tasks being used.
!***  When there are "remainder" layers, give them one at a time
!***  to each task in the row until they are used up.
!-----------------------------------------------------------------------
!
      allocate(k1_fft(0:npes-1),stat=istat)
      allocate(k2_fft(0:npes-1),stat=istat)
!
      npes_south=ipe_end_south-ipe_start_south+1
      lyr_frac_south=lm/npes_south
      l_remain=lm-npes_south*lyr_frac_south
!
      k2=0
      do npe=0,ipe_end_south
        k1_fft(npe)=k2+1
        k2=k1_fft(npe)+lyr_frac_south-1
        if(l_remain>0)then
          k2=k2+1
          l_remain=l_remain-1
        endif
        k2_fft(npe)=k2
      enddo
!
      npes_north=ipe_end_north-ipe_start_north+1
      lyr_frac_north=lm/npes_north
      l_remain=lm-npes_north*lyr_frac_north
!
      k2=0
      do npe=ipe_start_north,ipe_end_north
        k1_fft(npe)=k2+1
        k2=k1_fft(npe)+lyr_frac_north-1
        if(l_remain>0)then
          k2=k2+1
          l_remain=l_remain-1
        endif
        k2_fft(npe)=k2
      enddo
!
      lm_fft=max(lyr_frac_south,lyr_frac_north)+1
!
!-----------------------------------------------------------------------
!***  Allocate the working arrays for the FFTs depending on which
!***  hemisphere this task is in.
!-----------------------------------------------------------------------
!
      if(fft_south)then
        allocate(hn(ids:ide,jh_start_fft_south:jh_end_fft_south,1:lm_fft) &
                ,stat=istat)
        allocate(wn(ids:ide,jv_start_fft_south:jv_end_fft_south,1:lm_fft) &
                ,stat=istat)
        allocate(un(ids:ide,jv_start_fft_south:jv_end_fft_south,1:lm_fft) &
                ,stat=istat)
        allocate(vn(ids:ide,jv_start_fft_south:jv_end_fft_south,1:lm_fft) &
                ,stat=istat)
      elseif(fft_north)then
        allocate(hn(ids:ide,jh_start_fft_north:jh_end_fft_north,1:lm_fft) &
                ,stat=istat)
        allocate(wn(ids:ide,jv_start_fft_north:jv_end_fft_north,1:lm_fft) &
                ,stat=istat)
        allocate(un(ids:ide,jv_start_fft_north:jv_end_fft_north,1:lm_fft) &
                ,stat=istat)
        allocate(vn(ids:ide,jv_start_fft_north:jv_end_fft_north,1:lm_fft) &
                ,stat=istat)
      endif
!
!-----------------------------------------------------------------------
!***  Maximum size of dummy space for FFTs
!-----------------------------------------------------------------------
!
      imax=(ide-ids+1)/inpes+1
      jmax_south=max(jh_end_fft_south-jh_start_fft_south &
                    ,jv_end_fft_south-jv_start_fft_south)+1
      jmax_north=max(jh_end_fft_north-jh_start_fft_north &
                    ,jv_end_fft_north-jv_start_fft_north)+1
      jmax=max(jmax_south,jmax_north)
      msize_dummy_fft=imax*jmax*lm_fft
!
!-----------------------------------------------------------------------
!
                        endsubroutine prefft
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine fftfhn &
(km &
,khfilt &
,hfilt &
,field_h &
,craux1,craux2,craux3,rcaux1,rcaux2,rcaux3 &
,npes)

!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 km &
,npes        ! Number of compute tasks

integer(kind=kint),dimension(jds:jde),intent(in):: &
 khfilt

real(kind=kfpt),dimension(ids:(ide-3)/2+1,jds:jde),intent(in):: &
 hfilt

real(kind=kdbl),dimension(1:25000),intent(in) :: &
 craux1 &
,rcaux1

real(kind=kdbl),dimension(1:20000),intent(in) :: &
 craux2 &
,rcaux2
 
real(kind=kdbl),dimension(1:1),intent(in) :: &
 craux3 &
,rcaux3

real(kind=kfpt),dimension(ims:ime,jms:jme,1:km),intent(inout):: &
 field_h

!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------

integer(kind=kint) :: &
 i &
,icycle &
,icyclc &
,ipe_end &
,ipe_start &
,j &
,jend &
,jh_end_fft &
,jh_start_fft &
,jstart &
,k1 &
,k2 &
,l &
,n

real(kind=kfpt) :: &
 an &
,as &
,rcycle

!real(kind=kfpt),dimension(ids:ide,jts:jte,1:lm_fft) :: &
!real(kind=kfpt),dimension(ids:ide,jh_start_fft:jh_end_fft,1:lm_fft) :: &
! hn

real(kind=kfpt),dimension(ids:ide-3) :: &
 rarg
                                                                                                                                              
complex(kind=kfpt),dimension(ids:(ide-3)/2+1) :: &
 carg
integer(kind=kint) :: &
 nend
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
      icycle=ide-3
      icyclc=icycle/2+1
      rcycle=1./icycle
!-----------------------------------------------------------------------
!***  k1 and k2 are starting/ending vertical indices of model layers
!***  that this task will use in applying FFTs.
!-----------------------------------------------------------------------
!
      k1=k1_fft(mype)
      k2=k2_fft(mype)
!
!-----------------------------------------------------------------------
!***  Select hemisphere-dependent variables.
!-----------------------------------------------------------------------
!
      if(fft_south)then
        jh_start_fft=jh_start_fft_south
        jh_end_fft=jh_end_fft_south
        ipe_start=ipe_start_south
        ipe_end=ipe_end_south
      elseif(fft_north)then
        jh_start_fft=jh_start_fft_north
        jh_end_fft=jh_end_fft_north
        ipe_start=ipe_start_north
        ipe_end=ipe_end_north
      endif
!
!-----------------------------------------------------------------------
!***  Gather the model layer data from full latitude circles
!***  onto appropriate tasks for the FFTs.
!-----------------------------------------------------------------------
!
      call gather_layers(field_h,km,npes,msize_dummy_fft &
                        ,lm_fft,k1_fft,k2_fft &
                        ,local_istart,local_iend &
                        ,local_jstart,local_jend &
                        ,jh_start_fft,jh_end_fft &
                        ,ipe_start,ipe_end &
                        ,my_domain_has_fft_lats_h &
                        ,hn)
!
!-----------------------------------------------------------------------
!
!jaa      n=0
      nend=k2-k1+1
!jaa      kloop1: do l=k1,k2
      kloop1: do n=1,nend
!
!-----------------------------------------------------------------------
!jaa        n=n+1
!-----------------------------------------------------------------------
        if(fft_south)then
        as=0.
!-----------------------------------------------------------------------
          do i=ids+1,ide-2
            as=hn(i,jds+1,n)+as
          enddo
!
          as=as*rcycle
!
          do i=ids,ide
            hn(i,jds+1,n)=as
          enddo
!
!-----------------------------------------------------------------------
        elseif(fft_north)then
!-----------------------------------------------------------------------
        an=0.
          do i=ids+1,ide-2
            an=hn(i,jde-1,n)+an
          enddo
!
          an=an*rcycle
!
          do i=ids,ide
            hn(i,jde-1,n)=an
          enddo
!
!-----------------------------------------------------------------------
        endif   
!-----------------------------------------------------------------------
!
      enddo kloop1
!jaa      nend = n
!
!-----------------------------------------------------------------------
!***  jstart and jend are the starting/ending rows on which this task
!***  will apply FFTs
!-----------------------------------------------------------------------
!
      jstart=max(jh_start_fft,jds+2)
      jend=min(jh_end_fft,jde-2)
!
!-----------------------------------------------------------------------
!
!jaa      n=0
!jaa      kloop2: do l=k1,k2
!
!.......................................................................
!$omp parallel do private (n,j,i,rarg,carg,rcaux2,craux2)
!.......................................................................
!
      kloop2: do n=1,nend
!
!-----------------------------------------------------------------------
!jaa        n=n+1

        do j=jstart,jend
!-----------------------------------------------------------------------
          if(khfilt(j)<=icyclc) then
!-----------------------------------------------------------------------
            do i=ids+1,ide-2
              rarg(i-1)=hn(i,j,n)
            enddo
!
            call srcft(0,rarg,icycle,carg,icyclc,icycle,1,+1,1.     &
                      ,rcaux1,25000,rcaux2,20000,rcaux3,1)
!
            do i=1,khfilt(j)
              carg(i)=carg(i)*hfilt(i,j)
            enddo
!
            do i=khfilt(j)+1,icyclc
              carg(i)=0.
            enddo
!
            call scrft(0,carg,icyclc,rarg,icycle,icycle,1,-1,rcycle &
     &                ,craux1,25000,craux2,20000,craux3,1)
!
            do i=ids+1,ide-2
              hn(i,j,n)=rarg(i-1)
            enddo
!
            hn(ide-1,j,n)=rarg(1)
!-----------------------------------------------------------------------
          endif
!-----------------------------------------------------------------------
        enddo
!-----------------------------------------------------------------------
!
      enddo kloop2
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  Now scatter the model layer data from full latitude circles
!***  back to the appropriate tasks.
!-----------------------------------------------------------------------
!
      call scatter_layers(hn,km,npes,msize_dummy_fft &
                         ,lm_fft,k1_fft,k2_fft &
                         ,local_istart,local_iend &
                         ,local_jstart,local_jend &
                         ,jh_start_fft,jh_end_fft &
                         ,ipe_start,ipe_end &
                         ,my_domain_has_fft_lats_h &
                         ,field_h)
!!!   call mpi_abort(mpi_intra,1,ierr)
!
!-----------------------------------------------------------------------
!
                        endsubroutine fftfhn
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine fftfwn &
(km &
,kvfilt &
,vfilt &
,field_w &
,craux1,craux2,craux3,rcaux1,rcaux2,rcaux3 &
,npes)
!
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------

integer(kind=kint),intent(in) :: &
 km &
,npes      ! Number of compute tasks

integer(kind=kint),dimension(jds:jde),intent(in):: &
 kvfilt

real(kind=kfpt),dimension(ids:(ide-3)/2+1,jds:jde),intent(in):: &
 vfilt

real(kind=kdbl),dimension(1:25000),intent(in) :: &
 craux1 &
,rcaux1
 
real(kind=kdbl),dimension(1:20000),intent(in) :: &
 craux2 &
,rcaux2
 
real(kind=kdbl),dimension(1:1),intent(in) :: &
 craux3 &
,rcaux3

real(kind=kfpt),dimension(ims:ime,jms:jme,1:km),intent(inout):: &
 field_w

!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------

integer(kind=kint) :: &
 i &
,ipe_start &
,ipe_end &
,icycle &
,icyclc &
,iloc_mype &
,j &
,jend &
,jstart &
,jv_end_fft &
,jv_start_fft &
,k1 &
,k2 &
,l &
,n

real(kind=kfpt) :: &
 rcycle

!real(kind=kfpt),dimension(ids:ide,jts:jte,1:lm_fft) :: &
!real(kind=kfpt),dimension(ids:ide,jv_start_fft:jv_end_fft,1:lm_fft) :: &
! wn

real(kind=kfpt),dimension(ids:ide-3) :: &
 rarg
                                                                                                                                              
complex(kind=kfpt),dimension(ids:(ide-3)/2+1) :: &
 carg
integer(kind=kint) :: &
 nend 
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
      icycle=ide-3
      icyclc=icycle/2+1
      rcycle=1./icycle
!-----------------------------------------------------------------------
!***  k1 and k2 are starting/ending vertical indices of model layers
!***  that this task will use in applying FFTs.
!-----------------------------------------------------------------------
!
      k1=k1_fft(mype)
      k2=k2_fft(mype)
!
!-----------------------------------------------------------------------
!***  Select hemisphere-dependent variables.
!-----------------------------------------------------------------------
!
      if(fft_south)then
        jv_start_fft=jv_start_fft_south
        jv_end_fft=jv_end_fft_south
        ipe_start=ipe_start_south
        ipe_end=ipe_end_south
      elseif(fft_north)then
        jv_start_fft=jv_start_fft_north
        jv_end_fft=jv_end_fft_north
        ipe_start=ipe_start_north
        ipe_end=ipe_end_north
      endif
!
!-----------------------------------------------------------------------
!***  Gather the model layer data from full latitude circles
!***  onto appropriate tasks for the FFTs.
!-----------------------------------------------------------------------
!
      call gather_layers(field_w,km,npes,msize_dummy_fft &
                        ,lm_fft,k1_fft,k2_fft &
                        ,local_istart,local_iend &
                        ,local_jstart,local_jend &
                        ,jv_start_fft,jv_end_fft &
                        ,ipe_start,ipe_end &
                        ,my_domain_has_fft_lats_v &
                        ,wn)
!
!-----------------------------------------------------------------------
!***  jstart and jend are the starting/ending rows on which this task
!***  will apply FFTs
!-----------------------------------------------------------------------
!
      jstart=max(jv_start_fft,jds+1)
      jend=min(jv_end_fft,jde-2)
!
!-----------------------------------------------------------------------
!
!jaa      n=0
      nend=k2-k1+1
!jaa      kloop: do l=k1,k2
!.......................................................................
!$omp parallel do private (n,j,i,rarg,carg,rcaux2,craux2)
!.......................................................................
      kloop1: do n=1,nend
!
!jaa        n=n+1
!
        do j=jstart,jend
!-----------------------------------------------------------------------
          if(kvfilt(j)<=icycle) then
!-----------------------------------------------------------------------
            do i=ids+1,ide-2
              rarg(i-1)=wn(i,j,n)
            enddo
!
            call srcft(0,rarg,icycle,carg,icyclc,icycle,1,+1,1.     &
                      ,rcaux1,25000,rcaux2,20000,rcaux3,1)
!
            do i=1,kvfilt(j)
              carg(i)=carg(i)*vfilt(i,j)
            enddo
!
            do i=kvfilt(j)+1,icyclc
              carg(i)=0.
            enddo
!
            call scrft(0,carg,icyclc,rarg,icycle,icycle,1,-1,rcycle &
                      ,craux1,25000,craux2,20000,craux3,1)
!
            do i=ids+1,ide-2
              wn(i,j,n)=rarg(i-1)
            enddo
!
            wn(ide-1,j,n)=rarg(1)
!-----------------------------------------------------------------------
          endif
!-----------------------------------------------------------------------
        enddo
!-----------------------------------------------------------------------
!
      enddo kloop1
!.......................................................................
!$omp end parallel do
!.......................................................................
!
!-----------------------------------------------------------------------
!***  Now scatter the model layer data from full latitude circles
!***  back to the appropriate tasks.
!-----------------------------------------------------------------------
!
      call scatter_layers(wn,km,npes,msize_dummy_fft &
                         ,lm_fft,k1_fft,k2_fft &
                         ,local_istart,local_iend &
                         ,local_jstart,local_jend &
                         ,jv_start_fft,jv_end_fft &
                         ,ipe_start,ipe_end &
                         ,my_domain_has_fft_lats_v &
                         ,field_w)
!
!-----------------------------------------------------------------------
!
                        endsubroutine fftfwn
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine fftfuvn &
(km &
,kvfilt &
,vfilt &
,u,v &
,craux1,craux2,craux3,rcaux1,rcaux2,rcaux3 &
,npes)
!
!-----------------------------------------------------------------------
      implicit none
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------
                                                                                                                                              
integer(kind=kint),intent(in) :: &
 km &
,npes      ! Number of compute tasks
                                                                                                                                              
integer(kind=kint),dimension(jds:jde),intent(in):: &
 kvfilt
                                    
real(kind=kfpt),dimension(ids:(ide-3)/2+1,jds:jde),intent(in):: &
 vfilt
                                    
real(kind=kdbl),dimension(1:25000),intent(in) :: &
 craux1 &
,rcaux1
 
real(kind=kdbl),dimension(1:20000),intent(in) :: &
 craux2 &
,rcaux2
 
real(kind=kdbl),dimension(1:1),intent(in) :: &
 craux3 &
,rcaux3

real(kind=kfpt),dimension(ims:ime,jms:jme,1:km),intent(inout):: &
 u &
,v
                                    
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
                                    
integer(kind=kint) :: &
 i &
,ipe_start &
,ipe_end &
,icycle &
,icyclc &
,iloc_mype &
,j &
,jend &
,jstart &
,jv_end_fft &
,jv_start_fft &
,k1 &
,k2 &
,l &
,n
                                    
real(kind=kfpt) :: &
 anu &
,anv &
,asu &
,asv &
,rcycle
     
real(kind=kfpt),dimension(1:ide-3):: &
 rargu &
,rargv
       
complex(kind=kfpt),dimension(ids:(ide-3)/2+1) :: &
 cargu &
,cargv

!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      icycle=ide-3
      icyclc=icycle/2+1
      rcycle=1./icycle
!-----------------------------------------------------------------------
!***  k1 and k2 are starting/ending vertical indices of model layers
!***  that this task will use in applying FFTs.
!-----------------------------------------------------------------------
!
      k1=k1_fft(mype)
      k2=k2_fft(mype)
!
!-----------------------------------------------------------------------
!***  Select hemisphere-dependent variables.
!-----------------------------------------------------------------------
!
      if(fft_south)then
        jv_start_fft=jv_start_fft_south
        jv_end_fft=jv_end_fft_south
        ipe_start=ipe_start_south
        ipe_end=ipe_end_south
!!!   endif
!!!   if(fft_north)then
      elseif(fft_north)then
        jv_start_fft=jv_start_fft_north
        jv_end_fft=jv_end_fft_north
        ipe_start=ipe_start_north
        ipe_end=ipe_end_north
      endif
!
!-----------------------------------------------------------------------
!***  Gather the model layer data from full latitude circles
!***  onto appropriate tasks for the FFTs.
!-----------------------------------------------------------------------
!
      call gather_layers(u,km,npes,msize_dummy_fft &
                        ,lm_fft,k1_fft,k2_fft &
                        ,local_istart,local_iend &
                        ,local_jstart,local_jend &
                        ,jv_start_fft,jv_end_fft &
                        ,ipe_start,ipe_end &
                        ,my_domain_has_fft_lats_v &
                        ,un)
      call gather_layers(v,km,npes,msize_dummy_fft &
                        ,lm_fft,k1_fft,k2_fft &
                        ,local_istart,local_iend &
                        ,local_jstart,local_jend &
                        ,jv_start_fft,jv_end_fft &
                        ,ipe_start,ipe_end &
                        ,my_domain_has_fft_lats_v &
                        ,vn)
!
!-----------------------------------------------------------------------
!
!      n=0
!      kloop1: do l=k1,k2
!
!-----------------------------------------------------------------------
!        asu=0.
!        asv=0.
!        anu=0.
!        anv=0.
!        n=n+1
!-----------------------------------------------------------------------
!        if(fft_south)then
!-----------------------------------------------------------------------
!          do i=ids+1,ide-2
!            asu=un(i,jds+1,n)+asu
!            asv=vn(i,jds+1,n)+asv
!          enddo
!
!          asu=asu*rcycle
!          asv=asv*rcycle
!
!          do i=ids,ide
!            un(i,jds+1,n)=un(i,jds+1,n)-asu
!            vn(i,jds+1,n)=vn(i,jds+1,n)-asv
!          enddo
!-----------------------------------------------------------------------
!        elseif(fft_north)then
!-----------------------------------------------------------------------
!          do i=ids+1,ide-2
!            anu=un(i,jde-2,n)+anu
!            anv=vn(i,jde-2,n)+anv
!          enddo
!
!          anu=anu*rcycle
!          anv=anv*rcycle
!
!          do i=ids,ide
!            un(i,jde-2,n)=un(i,jde-2,n)-anu
!            vn(i,jde-2,n)=vn(i,jde-2,n)-anv
!          enddo
!-----------------------------------------------------------------------
!        endif
!-----------------------------------------------------------------------
!
!      enddo kloop1
!
!-----------------------------------------------------------------------
!***  jstart and jend are the starting/ending rows on which this task
!***  will apply FFTs
!-----------------------------------------------------------------------
!
      jstart=max(jv_start_fft,jds+1)
      jend=min(jv_end_fft,jde-2)
!
!-----------------------------------------------------------------------
!
      n=0
      kloop: do l=k1,k2
!
        n=n+1
!
        do j=jstart,jend
!-----------------------------------------------------------------------
          if(kvfilt(j).le.icyclc) then
!-----------------------------------------------------------------------
            do i=ids+1,ide-2
              rargu(i-1)=un(i,j,n)
              rargv(i-1)=vn(i,j,n)
            enddo
!
            call srcft(0,rargu,icycle,cargu,icyclc,icycle,1,+1,1.     &
                      ,rcaux1,25000,rcaux2,20000,rcaux3,1)
            call srcft(0,rargv,icycle,cargv,icyclc,icycle,1,+1,1.     &
                      ,rcaux1,25000,rcaux2,20000,rcaux3,1)
!
            do i=1,kvfilt(j)
              cargu(i)=cargu(i)*vfilt(i,j)
              cargv(i)=cargv(i)*vfilt(i,j)
            enddo
!
            do i=kvfilt(j)+1,icyclc
              cargu(i)=0.
              cargv(i)=0.
            enddo
!
            call scrft(0,cargu,icyclc,rargu,icycle,icycle,1,-1,rcycle &
                      ,craux1,25000,craux2,20000,craux3,1)
            call scrft(0,cargv,icyclc,rargv,icycle,icycle,1,-1,rcycle &
                      ,craux1,25000,craux2,20000,craux3,1)
!
            do i=ids+1,ide-2
              un(i,j,n)=rargu(i-1)
              vn(i,j,n)=rargv(i-1)
            enddo
!
            un(ide-1,j,n)=rargu(1)
            vn(ide-1,j,n)=rargv(1)
!-----------------------------------------------------------------------
          endif
!-----------------------------------------------------------------------
        enddo
!-----------------------------------------------------------------------
      enddo kloop
!-----------------------------------------------------------------------
!***  Now scatter the model layer data from full latitude circles
!***  back to the appropriate tasks.
!-----------------------------------------------------------------------
!
      call scatter_layers(un,km,npes,msize_dummy_fft &
                         ,lm_fft,k1_fft,k2_fft &
                         ,local_istart,local_iend &
                         ,local_jstart,local_jend &
                         ,jv_start_fft,jv_end_fft &
                         ,ipe_start,ipe_end &
                         ,my_domain_has_fft_lats_v &
                         ,u)
      call scatter_layers(vn,km,npes,msize_dummy_fft &
                         ,lm_fft,k1_fft,k2_fft &
                         ,local_istart,local_iend &
                         ,local_jstart,local_jend &
                         ,jv_start_fft,jv_end_fft &
                         ,ipe_start,ipe_end &
                         ,my_domain_has_fft_lats_v &
                         ,v)
!
!-----------------------------------------------------------------------
!
                        endsubroutine fftfuvn
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine poavhn &
(i_start,i_end,j_start,j_end,km,hn,ntsd)
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 i_start &
,i_end &
,j_start &
,j_end &
,km &
,ntsd
!
real(kind=kfpt),dimension(i_start:i_end,j_start:j_end,km),intent(inout):: &
 hn
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,ierr &
,irecv &
,l
!
real(kind=kfpt) :: &
 rcycle
!
real(kind=kfpt),dimension(1:km) :: &
 an &
,an_g &
,as &
,as_g
 
integer(kind=kint) :: istat
logical(kind=klog) :: opened
logical(kind=klog),save :: sum_file_is_open=.false.
character(10) :: fstatus
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!
      rcycle=1./(ide-3)
      do l=1,km
        as(l)=0.
        an(l)=0.
      enddo
!
!-----------------------------------------------------------------------
!***  Sum the values along the row north of the southern boundary.
!-----------------------------------------------------------------------
!
      if(s_bdy)then
        do l=1,km
          do i=its_b1,ite_b2
            as(l)=hn(i,jts+1,l)+as(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Generate the global sum of the previous sums from each task
!***  if they are not being read in.
!-----------------------------------------------------------------------
!
      if(.not.read_global_sums)then
!
        call mpi_allreduce(as,as_g,km,mpi_real,mpi_sum,mpi_comm_comp &
                          ,irecv)
      endif
!-----------------------------------------------------------------------
!***  For bit reproducibility, read/write global sums.
!-----------------------------------------------------------------------
!
      bits_1: if(read_global_sums.or.write_global_sums)then
        if(.not.sum_file_is_open.and.mype==0)then
          open_unit: do l=51,59
            inquire(l,opened=opened)
            if(.not.opened)then
              iunit_pole_sums=l
              if(read_global_sums)fstatus='OLD'
              if(write_global_sums)fstatus='REPLACE'
              open(unit=iunit_pole_sums,file='global_pole_sums',status=fstatus &
                  ,form='UNFORMATTED',iostat=istat)
              sum_file_is_open=.true.
              exit open_unit
            endif
          enddo open_unit
        endif
!
!***  Read in south/north global sums.
!
        if(read_global_sums)then
          if(mype==0)then
            do l=1,km
              read(iunit_pole_sums)as_g(l),an_g(l)
            enddo
          endif
!
          call mpi_bcast(as_g,km,mpi_real,0,mpi_comm_comp,ierr)
          call mpi_bcast(an_g,km,mpi_real,0,mpi_comm_comp,ierr)
!
        endif
!
      endif bits_1
!
!-----------------------------------------------------------------------
!***  Reset the array values in that same row to the global sum.
!-----------------------------------------------------------------------
!
      if(s_bdy)then
        do l=1,km
          as(l)=as_g(l)*rcycle
          do i=its,ite
            hn(i,jts+1,l)=as(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Now sum the values along the row south of the northern boundary.
!-----------------------------------------------------------------------
!
      if(n_bdy)then
        do l=1,km
          do i=its_b1,ite_b2
            an(l)=hn(i,jte-1,l)+an(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Generate the global sum of the previous sums from each task
!***  if they are not being read in.
!-----------------------------------------------------------------------
!
      if(.not.read_global_sums)then
!
        call mpi_allreduce(an,an_g,km,mpi_real,mpi_sum,mpi_comm_comp &
                          ,irecv)
      endif
!-----------------------------------------------------------------------
!***  For bit reproducibility, write global sums.
!-----------------------------------------------------------------------
!
      bits_2: if(write_global_sums)then
!
        if(mype==0)then
          do l=1,km
            write(iunit_pole_sums)as_g(l),an_g(l)
          enddo
        endif
!
      endif bits_2
!
!-----------------------------------------------------------------------
!***  Reset the array values in that same row to the global sum.
!-----------------------------------------------------------------------
!
      if(n_bdy)then
        do l=1,km
          an(l)=an_g(l)*rcycle
          do i=its,ite
            hn(i,jte-1,l)=an(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!
                        endsubroutine poavhn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine nhpoav &
(i_start,i_end,j_start,j_end,km,hn,ntsd)
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 i_start &
,i_end &
,j_start &
,j_end &
,km &
,ntsd
!
real(kind=kfpt),dimension(i_start:i_end,j_start:j_end,km),intent(inout):: &
 hn
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,ierr &
,irecv &
,l
!
real(kind=kfpt) :: &
 rcycle
!
real(kind=kfpt),dimension(1:km) :: &
 an &
,an_g &
,as &
,as_g
 
integer(kind=kint) :: istat
logical(kind=klog) :: opened
character(10) :: fstatus
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!
      rcycle=1./(ide-3)
      do l=1,km
        as(l)=0.
        an(l)=0.
      enddo
!
!-----------------------------------------------------------------------
!***  Sum the values two rows north of the southern boundary.
!-----------------------------------------------------------------------
!
      if(s_bdy)then
        do l=1,km
          do i=its_b1,ite_b2
            as(l)=hn(i,jts+2,l)+as(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Generate the global sum of the previous sums from each task
!***  if they are not being read in.
!-----------------------------------------------------------------------
!
      if(.not.read_global_sums)then
!
        call mpi_allreduce(as,as_g,km,mpi_real,mpi_sum,mpi_comm_comp &
                          ,irecv)
      endif
!-----------------------------------------------------------------------
!***  For bit reproducibility, read/write global sums.
!-----------------------------------------------------------------------
!
      bits_1: if(read_global_sums.or.write_global_sums)then
!
!***  Read in south/north global sums.
!
        if(read_global_sums)then
          if(mype==0)then
            do l=1,km
              read(iunit_pole_sums)as_g(l),an_g(l)
            enddo
          endif
!
          call mpi_bcast(as_g,km,mpi_real,0,mpi_comm_comp,ierr)
          call mpi_bcast(an_g,km,mpi_real,0,mpi_comm_comp,ierr)
!
        endif
!
      endif bits_1
!
!-----------------------------------------------------------------------
!***  Reset the array values one row north of the southern boundary row
!***  to the global sum of values two rows north of the boundary.
!-----------------------------------------------------------------------
!
      if(s_bdy)then
        do l=1,km
          as(l)=as_g(l)*rcycle
          do i=its,ite
            hn(i,jts+1,l)=as(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Now sum the values two rows south of the northern boundary.
!-----------------------------------------------------------------------
!
      if(n_bdy)then
        do l=1,km
          do i=its_b1,ite_b2
            an(l)=hn(i,jte-2,l)+an(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  Generate the global sum of the previous sums from each task
!***  if they are not being read in.
!-----------------------------------------------------------------------
!
      if(.not.read_global_sums)then
!
        call mpi_allreduce(an,an_g,km,mpi_real,mpi_sum,mpi_comm_comp &
                          ,irecv)
      endif
!-----------------------------------------------------------------------
!***  For bit reproducibility, write global sums.
!-----------------------------------------------------------------------
!
      bits_2: if(write_global_sums)then
!
        if(mype==0)then
          do l=1,km
            write(iunit_pole_sums)as_g(l),an_g(l)
          enddo
        endif
!
      endif bits_2
!
!-----------------------------------------------------------------------
!***  Reset the array values one row south of the northern boundary row
!***  to the global sum of values two rows south of the boundary.
!-----------------------------------------------------------------------
!
      if(n_bdy)then
        do l=1,km
          an(l)=an_g(l)*rcycle
          do i=its,ite
            hn(i,jte-1,l)=an(l)
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!
                        endsubroutine nhpoav
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine swaphn &
(hn,i_start,i_end,j_start,j_end,km,inpes)
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------
integer(kind=kint),intent(in) :: &
 i_start &
,i_end &
,j_start &
,j_end &
,km &
,inpes
!
real(kind=kfpt),dimension(i_start:i_end,j_start:j_end,km),intent(inout):: &
 hn
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 irecv &
,isend &
,j &
,l &
,length &
,ntask
!
integer(kind=kint),dimension(mpi_status_size) :: &
 jstat
!
real(kind=kfpt) :: &
 ave 

real(kind=kfpt),dimension(2,jts_h1:jte_h1,km) :: &
 eastx &
,westx
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!***  Each task along the eastern and western "boundaries" needs to
!***  exchange with its counterpart on the opposite side.
!
!***  What is being said here is that i=1 (west) and i=ite-2 (east)
!***  are coincident, i=3 (west) and i=ite (east) are coincident,
!***  and since i=2 (west) and i=ite-1 (east) are coincident and
!***  in the true integration domain then they will equal their 
!***  mean.
!-----------------------------------------------------------------------
!
      length=2*(jte_h1-jts_h1+1)*km
!
      if(e_bdy)then
        do l=1,km
        do j=jts_h1,jte_h1
          eastx(1,j,l)=hn(ite-2,j,l)
          eastx(2,j,l)=hn(ite-1,j,l)
        enddo
        enddo
!
        ntask=mype-inpes+1
        call mpi_send(eastx,length,mpi_real,ntask,mype &
                     ,mpi_comm_comp,isend)
!
        call mpi_recv(westx,length,mpi_real,ntask,ntask &
                     ,mpi_comm_comp,jstat,irecv)
!
        do l=1,km
        do j=jts_h1,jte_h1
          ave=(westx(1,j,l)+hn(ite-1,j,l))*0.5
          hn(ite-1,j,l)=ave
          hn(ite,j,l)=westx(2,j,l) 
        enddo
        enddo
!-----------------------------------------------------------------------
      elseif(w_bdy)then
!
        ntask=mype+inpes-1  
        call mpi_recv(eastx,length,mpi_real,ntask,ntask &
                     ,mpi_comm_comp,jstat,irecv)
!
        do l=1,km
        do j=jts_h1,jte_h1
          westx(1,j,l)=hn(its+1,j,l)
          westx(2,j,l)=hn(its+2,j,l)
        enddo
        enddo
!
        call mpi_send(westx,length,mpi_real,ntask,mype &
                     ,mpi_comm_comp,isend)
!
        do l=1,km
        do j=jts_h1,jte_h1
          hn(its,j,l)=eastx(1,j,l)
          ave=(hn(its+1,j,l)+eastx(2,j,l))*0.5
          hn(its+1,j,l)=ave
        enddo
        enddo
!
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
!
                        endsubroutine swaphn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine polehn &
(hn,i_start,i_end,j_start,j_end,km,inpes,jnpes)
!
!-----------------------------------------------------------------------
!***  Create polar mass point arrays holding all values around the
!***  southernmost and northernmost latitude circles in the
!***  true integration regions.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 i_end &
,i_start &
,inpes &
,j_end &
,j_start &
,jnpes &
,km
!
real(kind=kfpt),dimension(ims:ime,jms:jme,km),intent(inout):: &
 hn
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,iadd &
,ierr &
,ind &
,ipe &
,irecv &
,isend &
,istat &
,jpe &
,kount &
,l &
,npe_end &
,npe_start &
,num_pes &
,nwords_max &
,nwords_mine &
,nwords_remote
!
integer(kind=kint),dimension(2) :: &
 i_index
integer(kind=kint),dimension(4) :: &
 ihandle
!
real(kind=kfpt),dimension(ids-3:ide+3,km) :: &
 h_northpole &
,h_southpole
!
real(kind=kfpt),allocatable,dimension(:) :: &
 h_recv &
,h_send
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  First do the tasks at the south pole.
!-----------------------------------------------------------------------
!
      south_tasks: if(s_bdy)then
!
!----------------------------------------------------------------
!***  Allocate the arrays that will hold all local mass points in
!***  the latitude circle of interest.
!----------------------------------------------------------------
!
        nwords_max=(ide-ids+1)*km
        allocate(h_recv(1:nwords_max),stat=istat)
        allocate(h_send(1:nwords_max),stat=istat)
!
        nwords_mine=(ite-its+1)*km
        kount=0
!
!------------------------------------------------------------
!***  Each task inserts its local values into the full circle
!***  and also bundles those values.
!------------------------------------------------------------
!
        do l=1,km
        do i=its,ite
          h_southpole(i,l)=hn(i,3,l)
!
          kount=kount+1
          h_send(kount)=h_southpole(i,l)
        enddo
        enddo
!
!------------------------------------------------------
!***  Each task sends its word count and polar winds to
!***  the other pole tasks.
!------------------------------------------------------
!
        do ipe=0,inpes-1    !<-- senders
          do jpe=0,inpes-1  !<-- receivers
!
            if(jpe/=ipe)then
!
              if(mype==ipe)then
!
                i_index(1)=its
                i_index(2)=ite
!
                call mpi_isend(i_index,2,mpi_integer,jpe,ipe &              !<-- send my word count to other tasks
                              ,mpi_comm_comp,ihandle(1),isend)
                call mpi_wait(ihandle(1),istatw,ierr)
!
                call mpi_isend(h_send,nwords_mine,mpi_real,jpe,ipe &        !<-- send my mass points to other tasks
                              ,mpi_comm_comp,ihandle(3),isend)
                call mpi_wait(ihandle(3),istatw,ierr)
!
!------------------------------------------------------------
!***  Each task receives the word count and mass point values
!***  from the other pole tasks.
!------------------------------------------------------------
!
              elseif(mype==jpe)then
                call mpi_irecv(i_index,2,mpi_integer,ipe,ipe &              !<-- receive word count from other tasks
                              ,mpi_comm_comp,ihandle(2),irecv)
                call mpi_wait(ihandle(2),istatw,ierr)
!
                nwords_remote=(i_index(2)-i_index(1)+1)*km                  !<-- number of words sent by remote task
                call mpi_irecv(h_recv,nwords_remote,mpi_real,ipe,ipe &      !<-- receive mass points from other tasks
                              ,mpi_comm_comp,ihandle(4),irecv)
                call mpi_wait(ihandle(4),istatw,ierr)
!
!------------------------------------------------------------------
!***  Each task inserts the mass values on the full latitude circle
!***  from the other pole tasks.
!------------------------------------------------------------------
!
                kount=0
                do l=1,km
                do i=i_index(1),i_index(2)
                  kount=kount+1
                  h_southpole(i,l)=h_recv(kount)
                enddo
                enddo
!
              endif
!
            endif
!
          enddo
        enddo
!
!-----------------------------------------------------------
!***  Update mass points with values on the opposite side of
!***  of the latitude circle.
!-----------------------------------------------------------
!
        iadd=(ide-3)/2
        do l=1,km
        do i=its,ite
          ind=i+iadd
          if(ind>ide)ind=ind-ide+3
          hn(i,1,l)=h_southpole(ind,l)
        enddo
        enddo
!
        deallocate(h_send,h_recv)
!
!-----------------------------------------------------------------
!
      endif south_tasks
!
!-----------------------------------------------------------------
!***  Carry out the same procedure for the tasks at the north pole
!***  using three latitude circles.
!-----------------------------------------------------------------
!
      north_tasks: if(n_bdy)then
!
!----------------------------------------------------------------
!***  Allocate the arrays that will hold all local mass points in
!***  the latitude circle of interest.
!----------------------------------------------------------------
!
        nwords_max=(ide-ids+1)*km
        allocate(h_recv(1:nwords_max),stat=istat)
        allocate(h_send(1:nwords_max),stat=istat)
!
        nwords_mine=(ite-its+1)*km
        kount=0
!
!--------------------------------------------------------------
!***  Each task inserts its local values into the full circle
!***  and bundles those values.
!--------------------------------------------------------------
!
        do l=1,km
        do i=its,ite
          h_northpole(i,l)=hn(i,jde-2,l)
!
          kount=kount+1
          h_send(kount)=h_northpole(i,l)
        enddo
        enddo
!
!---------------------------------------------------------
!***  Each task sends its word count and polar mass points
!***  to the other pole tasks.
!---------------------------------------------------------
!
        num_pes=inpes*jnpes
        npe_start=num_pes-inpes
        npe_end=num_pes-1
!
        do ipe=npe_start,npe_end       !<-- senders
          do jpe=npe_start,npe_end     !<-- receivers
!
            if(jpe/=ipe)then
!
              if(mype==ipe)then
!
                i_index(1)=its
                i_index(2)=ite
!
                call mpi_isend(i_index,2,mpi_integer,jpe,ipe &              !<-- send my word count to other tasks
                              ,mpi_comm_comp,ihandle(1),isend)
                call mpi_wait(ihandle(1),istatw,ierr)
!
                call mpi_isend(h_send,nwords_mine,mpi_real,jpe,ipe &        !<-- send my mass points to other tasks
                              ,mpi_comm_comp,ihandle(3),isend)
                call mpi_wait(ihandle(3),istatw,ierr)
!
!------------------------------------------------------
!***  Each task receives the word count and mass points
!***  from the other pole tasks.
!------------------------------------------------------
!
              elseif(mype==jpe)then
                call mpi_irecv(i_index,2,mpi_integer,ipe,ipe &              !<-- receive word count from other tasks
                              ,mpi_comm_comp,ihandle(2),irecv)
                call mpi_wait(ihandle(2),istatw,ierr)
!
                nwords_remote=(i_index(2)-i_index(1)+1)*km                  !<-- number of words sent by remote task
                call mpi_irecv(h_recv,nwords_remote,mpi_real,ipe,ipe &      !<-- receive mass points from other tasks
                              ,mpi_comm_comp,ihandle(4),irecv)
                call mpi_wait(ihandle(4),istatw,ierr)
!
!------------------------------------------------------------------
!***  Each task inserts the mass values on the full latitude circle
!***  from the other pole tasks.
!------------------------------------------------------------------
!
                kount=0
                do l=1,km
                do i=i_index(1),i_index(2)
                  kount=kount+1
                  h_northpole(i,l)=h_recv(kount)
                enddo
                enddo
!
              endif
!
            endif
!
          enddo
        enddo
!
!--------------------------------------------------------
!***  Update mass points with values on the opposite side
!***  of the latitude circle.
!--------------------------------------------------------
!
        iadd=(ide-3)/2
        do l=1,km
        do i=its,ite
          ind=i+iadd
          if(ind>ide)ind=ind-ide+3
          hn(i,jte,l)=h_northpole(ind,l)
        enddo
        enddo
!
        deallocate(h_send,h_recv)
!
!-----------------------------------------------------------------------
!
      endif north_tasks
!
!-----------------------------------------------------------------------
!
                        end subroutine polehn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine swapwn &
(wn,i_start,i_end,j_start,j_end,km,inpes)
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 i_start &
,i_end &
,j_start &
,j_end &
,km &
,inpes
!
real(kind=kfpt),dimension(i_start:i_end,j_start:j_end,km),intent(inout):: &
 wn
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
integer(kind=kint) :: &
 irecv &
,isend &
,j &
,l &
,length_e &
,length_w &
,ntask
!
integer(kind=kint),dimension(mpi_status_size) :: &
 jstat
!
real(kind=kfpt),dimension(jts_h1:jte_h1,km) :: &
 eastx 
!
real(kind=kfpt),dimension(2,jts_h1:jte_h1,km) :: &
 westx
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
      mype=mype_share
!
!-----------------------------------------------------------------------
!***  Each task along the eastern and western "boundaries" needs to
!***  exchange with its counterpart on the opposite side.
!
!***  What is being said here is that i=1 (west) and i=ite-2 (east)
!***  are coincident, i=2 (west) and i=ite-1 (east) are coincident,
!***  and i=3 (west) and i=ite (east) are coincident.
!-----------------------------------------------------------------------
!
      length_e=(jte_h1-jts_h1+1)*km
      length_w=2*length_e
!
      if(e_bdy)then
        do l=1,km
        do j=jts_h1,jte_h1
          eastx(j,l)=wn(ite-2,j,l)
        enddo
        enddo
!
        ntask=mype-inpes+1
        call mpi_send(eastx,length_e,mpi_real,ntask,mype &
                     ,mpi_comm_comp,isend)
!
        call mpi_recv(westx,length_w,mpi_real,ntask,ntask &
                     ,mpi_comm_comp,jstat,irecv)
!
        do l=1,km
        do j=jts_h1,jte_h1
          wn(ite-1,j,l)=westx(1,j,l)
          wn(ite,j,l)=westx(2,j,l)
        enddo
        enddo
!-----------------------------------------------------------------------
      elseif(w_bdy)then
!
        ntask=mype+inpes-1
        call mpi_recv(eastx,length_e,mpi_real,ntask,ntask &
                     ,mpi_comm_comp,jstat,irecv)
!
        do l=1,km
        do j=jts_h1,jte_h1
          wn(its,j,l)=eastx(j,l)
        enddo
        enddo
!
        do l=1,km
        do j=jts_h1,jte_h1
          westx(1,j,l)=wn(its+1,j,l)
          westx(2,j,l)=wn(its+2,j,l)
        enddo
        enddo
!
        call mpi_send(westx,length_w,mpi_real,ntask,mype &
                     ,mpi_comm_comp,isend)
!
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
!
                        endsubroutine swapwn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine polewn &
(u,v,i_start,i_end,j_start,j_end,km,inpes,jnpes)
!
!-----------------------------------------------------------------------
!***  Create polar wind arrays holding all values around the
!***  southernmost and northernmost latitude circles in the
!***  true integration regions.
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!include 'kind.inc'
!-----------------------------------------------------------------------
!
integer(kind=kint),intent(in) :: &
 i_end &
,i_start &
,inpes &
,j_end &
,j_start &
,jnpes &
,km
!
real(kind=kfpt),dimension(ims:ime,jms:jme,km),intent(inout):: &
 u &
,v
!
!-----------------------------------------------------------------------
!***  Local Variables
!-----------------------------------------------------------------------
!
integer(kind=kint) :: &
 i &
,iadd &
,ierr &
,ind &
,ipe &
,irecv &
,isend &
,istat &
,jpe &
,kount &
,l &
,npe_end &
,npe_start &
,num_pes &
,nwords_max &
,nwords_mine &
,nwords_remote
!
integer(kind=kint),dimension(2) :: &
 i_index
integer(kind=kint),dimension(4) :: &
 ihandle
!
real(kind=kfpt),dimension(ids-3:ide+3,km) :: &
 u_northpole &
,u_southpole &
,v_northpole &
,v_southpole
!
real(kind=kfpt),allocatable,dimension(:) :: &
 wind_recv &
,wind_send
!
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!***  First do the tasks at the south pole.
!-----------------------------------------------------------------------
!
      south_tasks: if(s_bdy)then
!
!--------------------------------------------------------------------
!***  Allocate the arrays that will hold all local wind components in
!***  the latitude circle of interest.
!--------------------------------------------------------------------
!
        nwords_max=(ide-ids+1)*km*2
        allocate(wind_recv(1:nwords_max),stat=istat)
        allocate(wind_send(1:nwords_max),stat=istat)
!
        nwords_mine=(ite-its+1)*km*2
        kount=-1
!
!--------------------------------------------------------------
!***  Each task inserts its local values into the full circle
!***  and also bundles those values.
!--------------------------------------------------------------
!
        do l=1,km
        do i=its,ite
          u_southpole(i,l)=u(i,2,l)
          v_southpole(i,l)=v(i,2,l)
!
          kount=kount+2
          wind_send(kount)  =u_southpole(i,l)
          wind_send(kount+1)=v_southpole(i,l)
        enddo
        enddo
!
!------------------------------------------------------
!***  Each task sends its word count and polar winds to
!***  the other pole tasks.
!------------------------------------------------------
!
        do ipe=0,inpes-1    !<-- senders
          do jpe=0,inpes-1  !<-- receivers
!
            if(jpe/=ipe)then
!
              if(mype==ipe)then
!
                i_index(1)=its
                i_index(2)=ite
!
                call mpi_isend(i_index,2,mpi_integer,jpe,ipe &              !<-- send my word count to other tasks
                              ,mpi_comm_comp,ihandle(1),isend)
                call mpi_wait(ihandle(1),istatw,ierr)
!
                call mpi_isend(wind_send,nwords_mine,mpi_real,jpe,ipe &     !<-- send my winds to other tasks
                              ,mpi_comm_comp,ihandle(3),isend)
                call mpi_wait(ihandle(3),istatw,ierr)
!
!-----------------------------------------------------
!***  Each task receives the word count and winds from
!***  the other pole tasks.
!-----------------------------------------------------
!
              elseif(mype==jpe)then
                call mpi_irecv(i_index,2,mpi_integer,ipe,ipe &              !<-- receive word count from other tasks
                              ,mpi_comm_comp,ihandle(2),irecv)
                call mpi_wait(ihandle(2),istatw,ierr)
!
                nwords_remote=(i_index(2)-i_index(1)+1)*km*2                !<-- number of words sent by remote task
                call mpi_irecv(wind_recv,nwords_remote,mpi_real,ipe,ipe &   !<-- receive winds from other tasks
                              ,mpi_comm_comp,ihandle(4),irecv)
                call mpi_wait(ihandle(4),istatw,ierr)
!
!------------------------------------------------------------------
!***  Each task inserts the wind values on the full latitude circle
!***  from the other pole tasks.
!------------------------------------------------------------------
!
                kount=-1
                do l=1,km
                do i=i_index(1),i_index(2)
                  kount=kount+2
                  u_southpole(i,l)=wind_recv(kount)
                  v_southpole(i,l)=wind_recv(kount+1)
                enddo
                enddo
!
              endif
!
            endif
!
          enddo
        enddo
!
!-----------------------------------------------------
!***  Update winds with values on the opposite side of
!***  the latitude circle.
!-----------------------------------------------------
!
        iadd=(ide-3)/2
        do l=1,km
        do i=its,ite
          ind=i+iadd
          if(ind>ide)ind=ind-ide+3
          u(i,1,l)=-u_southpole(ind,l)
          v(i,1,l)=-v_southpole(ind,l)
        enddo
        enddo
!
        deallocate(wind_send,wind_recv)
!
!-----------------------------------------------------------------
!
      endif south_tasks
!
!-----------------------------------------------------------------
!***  Carry out the same procedure for the tasks at the north pole
!***  using three latitude circles.
!-----------------------------------------------------------------
!
      north_tasks: if(n_bdy)then
!
!--------------------------------------------------------------------
!***  Allocate the arrays that will hold all local wind components in
!***  the latitude circle of interest.
!--------------------------------------------------------------------
!
        nwords_max=(ide-ids+1)*km*2
        allocate(wind_recv(1:nwords_max),stat=istat)
        allocate(wind_send(1:nwords_max),stat=istat)
!
        nwords_mine=(ite-its+1)*km*2
        kount=-1
!
!--------------------------------------------------------------
!***  Each task inserts its local values into the full circle
!***  and bundles those values.
!--------------------------------------------------------------
!
        do l=1,km
        do i=its,ite
          u_northpole(i,l)=u(i,jde-2,l)
          v_northpole(i,l)=v(i,jde-2,l)
!
          kount=kount+2
          wind_send(kount)  =u_northpole(i,l)
          wind_send(kount+1)=v_northpole(i,l)
        enddo
        enddo
!
!------------------------------------------------------
!***  Each task sends its word count and polar winds to
!***  the other pole tasks.
!------------------------------------------------------
!
        num_pes=inpes*jnpes
        npe_start=num_pes-inpes
        npe_end=num_pes-1
!
        do ipe=npe_start,npe_end       !<-- senders
          do jpe=npe_start,npe_end     !<-- receivers
!
            if(jpe/=ipe)then
!
              if(mype==ipe)then
!
                i_index(1)=its
                i_index(2)=ite
!
                call mpi_isend(i_index,2,mpi_integer,jpe,ipe &              !<-- send my word count to other tasks
                              ,mpi_comm_comp,ihandle(1),isend)
                call mpi_wait(ihandle(1),istatw,ierr)
!
                call mpi_isend(wind_send,nwords_mine,mpi_real,jpe,ipe &     !<-- send my winds to other tasks
                              ,mpi_comm_comp,ihandle(3),isend)
                call mpi_wait(ihandle(3),istatw,ierr)
!
!-----------------------------------------------------
!***  Each task receives the word count and winds from
!***  the other pole tasks.
!-----------------------------------------------------
!
              elseif(mype==jpe)then
                call mpi_irecv(i_index,2,mpi_integer,ipe,ipe &              !<-- receive word count from other tasks
                              ,mpi_comm_comp,ihandle(2),irecv)
                call mpi_wait(ihandle(2),istatw,ierr)
!
                nwords_remote=(i_index(2)-i_index(1)+1)*km*2                !<-- number of words sent by remote task
                call mpi_irecv(wind_recv,nwords_remote,mpi_real,ipe,ipe &   !<-- receive winds from other tasks
                              ,mpi_comm_comp,ihandle(4),irecv)
                call mpi_wait(ihandle(4),istatw,ierr)
!
!------------------------------------------------------------------
!***  Each task inserts the wind values on the full latitude circle
!***  from the other pole tasks.
!------------------------------------------------------------------
!
                kount=-1
                do l=1,km
                do i=i_index(1),i_index(2)
                  kount=kount+2
                  u_northpole(i,l)=wind_recv(kount)
                  v_northpole(i,l)=wind_recv(kount+1)
                enddo
                enddo
!
              endif
!
            endif
!
          enddo
        enddo
!
!-----------------------------------------------------
!***  Update winds with values on the opposite side of
!***  the latitude circle.
!-----------------------------------------------------
!
        iadd=(ide-3)/2
        do l=1,km
        do i=its,ite
          ind=i+iadd
          if(ind>ide)ind=ind-ide+3
          u(i,jte-1,l)=-u_northpole(ind,l)
          u(i,jte  ,l)=-u_northpole(ind,l)
          v(i,jte-1,l)=-v_northpole(ind,l)
          v(i,jte  ,l)=-v_northpole(ind,l)
        enddo
        enddo
!
        deallocate(wind_send,wind_recv)
!
!-----------------------------------------------------------------------
!
      endif north_tasks
!
!-----------------------------------------------------------------------
!
                        end subroutine polewn
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        subroutine read_bc &
(lm,lnsh,lnsv,ntsd,dt &
,runbc,idatbc,ihrstbc,tboco &
,pdbs,pdbn,pdbw,pdbe &
,tbs,tbn,tbw,tbe &
,qbs,qbn,qbw,qbe &
,wbs,wbn,wbw,wbe &
,ubs,ubn,ubw,ube &
,vbs,vbn,vbw,vbe)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
implicit none
!
!-----------------------------------------------------------------------
!
!include 'kind.inc'
!
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm  &                       ! total # of levels
,lnsh &                      ! # of boundary h lines for bc in reg. setup
,lnsv &                      ! # of boundary v lines for bc in reg. setup
,ntsd                        ! current timestep

integer(kind=kint),intent(out):: &
 ihrstbc                     ! boundary conditions starting time

integer(kind=kint),dimension(1:3),intent(out):: &
 idatbc                      ! date of boundary data, day, month, year

real(kind=kfpt),intent(in):: &
 dt                          ! dynamics time step

real(kind=kfpt),intent(out):: &
 tboco                       ! boundary conditions interval, hours

real(kind=kfpt),dimension(ims:ime,1:lnsh,1:2),intent(out):: &
 pdbn &                      ! pressure difference at northern boundary
,pdbs                        ! pressure difference at southern boundary

real(kind=kfpt),dimension(1:lnsh,jms:jme,1:2),intent(out):: &
 pdbe &                      ! pressure difference at eastern boundary
,pdbw                        ! pressure difference at western boundary

real(kind=kfpt),dimension(ims:ime,1:lnsh,1:lm,1:2),intent(out):: &
 tbn &                       ! temperature at northern boundary
,tbs &                       ! temperature at southern boundary
,qbn &                       ! specific humidity at northern boundary
,qbs &                       ! specific humidity at southern boundary
,wbn &                       ! condensate at northern boundary
,wbs                         ! condensate at southern boundary

real(kind=kfpt),dimension(1:lnsh,jms:jme,1:lm,1:2),intent(out):: &
 tbe &                       ! temperature at eastern boundary
,tbw &                       ! temperature at western boundary
,qbe &                       ! specific humidity at eastern boundary
,qbw &                       ! specific humidity at western boundary
,wbe &                       ! condensate at eastern boundary
,wbw                         ! condensate at western boundary

real(kind=kfpt),dimension(ims:ime,1:lnsv,1:lm,1:2),intent(out):: &
 ubn &                       ! u wind component at northern boundary
,ubs &                       ! u wind component at southern boundary
,vbn &                       ! v wind component at northern boundary
,vbs                         ! v wind component at southern boundary

real(kind=kfpt),dimension(1:lnsv,jms:jme,1:lm,1:2),intent(out):: &
 ube &                       ! u wind component at eastern boundary
,ubw &                       ! u wind component at western boundary
,vbe &                       ! v wind component at eastern boundary
,vbw                         ! v wind component at western boundary

logical(kind=klog) :: runbc
!-----------------------------------------------------------------------
!---local variables-----------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,i_hi &                      ! max i loop limit (cannot be > ide)
,i_lo &                      ! min i loop limit (cannot be < ids)
,ihr &                       ! current hour
,ihrbc &                     ! current hour for bc's
,istat &                     ! return status
,iunit &                     ! file unit number
,j_hi &                      ! max j loop limit (cannot be > jde)
,j_lo &                      ! min j loop limit (cannot be < jds)
,l &
,n1 &                        ! dimension 1 of working arrays
,n2 &                        ! dimension 2 of working arrays
,n3 &                        ! dimension 3 of working arrays
,n4                          ! dimension 4 of working arrays

real(kind=kfpt),allocatable,dimension(:,:,:) :: pdbe_g,pdbn_g,pdbs_g,pdbw_g
 
real(kind=kfpt),allocatable,dimension(:,:,:,:) :: qbe_g,qbn_g,qbs_g,qbw_g &
                                                 ,tbe_g,tbn_g,tbs_g,tbw_g &
                                                 ,ube_g,ubn_g,ubs_g,ubw_g &
                                                 ,vbe_g,vbn_g,vbs_g,vbw_g & 
                                                 ,wbe_g,wbn_g,wbs_g,wbw_g 
character(64) :: infile
logical(kind=klog) :: opened
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!***  Allocate temporary working arrays.
!-----------------------------------------------------------------------
!
      allocate(pdbe_g(1:lnsh,jds:jde,1:2),stat=i)
      allocate(pdbw_g(1:lnsh,jds:jde,1:2),stat=i)
      allocate(pdbn_g(ids:ide,1:lnsh,1:2),stat=i)
      allocate(pdbs_g(ids:ide,1:lnsh,1:2),stat=i)
!
      allocate(tbe_g(1:lnsh,jds:jde,1:lm,1:2),stat=istat)
      allocate(tbw_g(1:lnsh,jds:jde,1:lm,1:2),stat=istat)
      allocate(tbn_g(ids:ide,1:lnsh,1:lm,1:2),stat=istat)
      allocate(tbs_g(ids:ide,1:lnsh,1:lm,1:2),stat=istat)
      allocate(qbe_g(1:lnsh,jds:jde,1:lm,1:2),stat=istat)
      allocate(qbw_g(1:lnsh,jds:jde,1:lm,1:2),stat=istat)
      allocate(qbn_g(ids:ide,1:lnsh,1:lm,1:2),stat=istat)
      allocate(qbs_g(ids:ide,1:lnsh,1:lm,1:2),stat=istat)
      allocate(wbe_g(1:lnsh,jds:jde,1:lm,1:2),stat=istat)
      allocate(wbw_g(1:lnsh,jds:jde,1:lm,1:2),stat=istat)
      allocate(wbn_g(ids:ide,1:lnsh,1:lm,1:2),stat=istat)
      allocate(wbs_g(ids:ide,1:lnsh,1:lm,1:2),stat=istat)
!
      allocate(ube_g(1:lnsv,jds:jde,1:lm,1:2),stat=istat)
      allocate(ubw_g(1:lnsv,jds:jde,1:lm,1:2),stat=istat)
      allocate(ubn_g(ids:ide,1:lnsv,1:lm,1:2),stat=istat)
      allocate(ubs_g(ids:ide,1:lnsv,1:lm,1:2),stat=istat)
      allocate(vbe_g(1:lnsv,jds:jde,1:lm,1:2),stat=istat)
      allocate(vbw_g(1:lnsv,jds:jde,1:lm,1:2),stat=istat)
      allocate(vbn_g(ids:ide,1:lnsv,1:lm,1:2),stat=istat)
      allocate(vbs_g(ids:ide,1:lnsv,1:lm,1:2),stat=istat)
!
!-----------------------------------------------------------------------
!***  Because subdomains that lie along the global domain boundary
!***  may have haloes that extend beyond the global limits, create
!***  limits here that keep loops from reaching beyond those
!***  global limits.
!-----------------------------------------------------------------------
!
      i_lo=max(ims,ids)
      i_hi=min(ime,ide)
      j_lo=max(jms,jds)
      j_hi=min(jme,jde)
!
!-----------------------------------------------------------------------
!***  Read in boundary variables and arrays. 
!-----------------------------------------------------------------------
!
      ihr=nint(ntsd*dt/3600.)
      ihrbc=ihr
      write(infile,'(a,i3.3)')'boco.',ihrbc
!
      select_unit: do l=51,59
        inquire(l,opened=opened)
        if(.not.opened)then
          iunit=l
          exit select_unit
        endif
      enddo select_unit
!
      open(unit=iunit,file=infile,status='OLD' &
          ,form='UNFORMATTED',iostat=istat)
!
      read(iunit)runbc,idatbc,ihrstbc,tboco
!
      read(iunit)pdbs_g,pdbn_g,pdbw_g,pdbe_g
!
      do n3=1,2
        do n2=1,lnsh
!       do n1=ims,ime
        do n1=i_lo,i_hi
          pdbn(n1,n2,n3)=pdbn_g(n1,n2,n3)
          pdbs(n1,n2,n3)=pdbs_g(n1,n2,n3)
        enddo
        enddo
!
!       do n2=jms,jme
        do n2=j_lo,j_hi
        do n1=1,lnsh
          pdbe(n1,n2,n3)=pdbe_g(n1,n2,n3)
          pdbw(n1,n2,n3)=pdbw_g(n1,n2,n3)
        enddo
        enddo
      enddo
!
      deallocate(pdbe_g,pdbn_g,pdbs_g,pdbw_g)
!
      read(iunit)tbs_g,tbn_g,tbw_g,tbe_g
      read(iunit)qbs_g,qbn_g,qbw_g,qbe_g
      read(iunit)wbs_g,wbn_g,wbw_g,wbe_g
!
      do n4=1,2
      do n3=1,lm
        do n2=1,lnsh
!       do n1=ims,ime
        do n1=i_lo,i_hi
          tbn(n1,n2,n3,n4)=tbn_g(n1,n2,n3,n4)
          tbs(n1,n2,n3,n4)=tbs_g(n1,n2,n3,n4)
          qbn(n1,n2,n3,n4)=qbn_g(n1,n2,n3,n4)
          qbs(n1,n2,n3,n4)=qbs_g(n1,n2,n3,n4)
          wbn(n1,n2,n3,n4)=wbn_g(n1,n2,n3,n4)
          wbs(n1,n2,n3,n4)=wbs_g(n1,n2,n3,n4)
        enddo
        enddo
!
!       do n2=jms,jme
        do n2=j_lo,j_hi
        do n1=1,lnsh
          tbe(n1,n2,n3,n4)=tbe_g(n1,n2,n3,n4)
          tbw(n1,n2,n3,n4)=tbw_g(n1,n2,n3,n4)
          qbe(n1,n2,n3,n4)=qbe_g(n1,n2,n3,n4)
          qbw(n1,n2,n3,n4)=qbw_g(n1,n2,n3,n4)
          wbe(n1,n2,n3,n4)=wbe_g(n1,n2,n3,n4)
          wbw(n1,n2,n3,n4)=wbw_g(n1,n2,n3,n4)
        enddo
        enddo
      enddo
      enddo
!
      deallocate(qbe_g)
      deallocate(qbn_g)
      deallocate(qbs_g)
      deallocate(qbw_g)
      deallocate(tbe_g)
      deallocate(tbn_g)
      deallocate(tbs_g)
      deallocate(tbw_g)
      deallocate(wbe_g)
      deallocate(wbn_g)
      deallocate(wbs_g)
      deallocate(wbw_g)
!
      read(iunit)ubs_g,ubn_g,ubw_g,ube_g
      read(iunit)vbs_g,vbn_g,vbw_g,vbe_g
!
      do n4=1,2
      do n3=1,lm
        do n2=1,lnsv
!       do n1=ims,ime
        do n1=i_lo,i_hi
          ubn(n1,n2,n3,n4)=ubn_g(n1,n2,n3,n4)       
          ubs(n1,n2,n3,n4)=ubs_g(n1,n2,n3,n4)       
          vbn(n1,n2,n3,n4)=vbn_g(n1,n2,n3,n4)       
          vbs(n1,n2,n3,n4)=vbs_g(n1,n2,n3,n4)       
        enddo
        enddo
!
!       do n2=jms,jme
        do n2=j_lo,j_hi
        do n1=1,lnsv
          ube(n1,n2,n3,n4)=ube_g(n1,n2,n3,n4)
          ubw(n1,n2,n3,n4)=ubw_g(n1,n2,n3,n4)
          vbe(n1,n2,n3,n4)=vbe_g(n1,n2,n3,n4)
          vbw(n1,n2,n3,n4)=vbw_g(n1,n2,n3,n4)
        enddo
        enddo
      enddo
      enddo
!
      deallocate(ube_g)
      deallocate(ubn_g)
      deallocate(ubs_g)
      deallocate(ubw_g)
      deallocate(vbe_g)
      deallocate(vbn_g)
      deallocate(vbs_g)
      deallocate(vbw_g)
!
      close(iunit)
!
!-----------------------------------------------------------------------
                        end subroutine read_bc
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine bocoh &
(lm,lnsh &
,dt,pt &
,dsg2,pdsg1 &
,pd &
,pdbe,pdbn,pdbs,pdbw &
,tbe,tbn,tbs,tbw,qbe,qbn,qbs,qbw,wbe,wbn,wbs,wbw &
,t,q,cw &
,pint)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------

!include 'kind.inc'
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 wa=0.5 &                    ! weighting factor
,w1=wa*0.25 &                ! weighting factor
,w2=1.-wa                    ! weighting factor
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,lnsh                        ! blending area width, h points

real(kind=kfpt),intent(in):: &
 dt &                        ! dynamics time step
,pt                          ! pressure at the top

real(kind=kfpt),dimension(1:lm),intent(in):: &
 dsg2 &                      ! delta sigmas
,pdsg1                       ! delta pressures

real(kind=kfpt),dimension(ims:ime,1:lnsh,1:2),intent(inout):: &
 pdbn &                      ! pressure difference at northern boundary
,pdbs                        ! pressure difference at southern boundary

real(kind=kfpt),dimension(1:lnsh,jms:jme,1:2),intent(inout):: &
 pdbe &                      ! pressure difference at eastern boundary
,pdbw                        ! pressure difference at western boundary

real(kind=kfpt),dimension(ims:ime,1:lnsh,1:lm,1:2),intent(inout):: &
 tbn &                       ! temperature at northern boundary
,tbs &                       ! temperature at southern boundary
,qbn &                       ! specific humidity at northern boundary
,qbs &                       ! specific humidity at southern boundary
,wbn &                       ! condensate at northern boundary
,wbs                         ! condensate at southern boundary

real(kind=kfpt),dimension(1:lnsh,jms:jme,1:lm,1:2),intent(inout):: &
 tbe &                       ! temperature at eastern boundary
,tbw &                       ! temperature at western boundary
,qbe &                       ! specific humidity at eastern boundary
,qbw &                       ! specific humidity at western boundary
,wbe &                       ! condensate at eastern boundary
,wbw                         ! condensate at western boundary

real(kind=kfpt),dimension(ims:ime,jms:jme),intent(inout):: &
 pd                          ! sigma range pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 cw &                        ! condensate
,q &                         ! specific humidity
,t                           ! temperature

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm+1),intent(inout):: &
 pint                        ! pressure at interfaces
!-----------------------------------------------------------------------
!---local variables-----------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,ib &                        ! index in x direction, boundary zone
,ihe &                       ! ending index in x direction, boundaries
,ihs &                       ! starting index in x direction, boundaries
,j &                         ! index in y direction
,jb &                        ! index in y direction, boundary zone
,jhe &                       ! ending index in x direction, boundaries
,jhs &                       ! starting index in x direction, boundaries
,k &                         ! boundary line counter
,ks &                        ! smoothing counter
,l &                         ! index in p direction
,lines &                     ! boundary smoothing area
,nsmud                       ! number of smoothing passes

real(kind=kfpt),dimension(1:lnsh):: &
 wh(lnsh) &                  ! blending weighting function, temperature
,wq(lnsh)                    ! blending weighting function, moisture

real(kind=kfpt),dimension(ims:ime,jms:jme):: &
 pdr                         ! pressure difference

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 pr &                        ! pressure
,qr &                        ! specific humidity
,tr                          ! temperature
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      lines=lnsh
      nsmud=lines-1
!-----------------------------------------------------------------------
      wh(1)=1.
      wq(1)=1.
      do k=2,lnsh
        wh(k)=exp(-(real(k)-1.)/(real(lnsh)-1.)*7.)
        wq(k)=exp(-(real(k)-1.)/(real(lnsh)-1.)*7.)
      enddo
!-----------------------------------------------------------------------
!-------------update of boundary values at h points---------------------
!-----------------------------------------------------------------------
!-------------southern and northern boundary----------------------------
!-----------------------------------------------------------------------
      if(s_bdy)then
        do jb=1,lnsh
          do ib=its_h2,ite_h2
            pdbs(ib,jb,1)=pdbs(ib,jb,1)+pdbs(ib,jb,2)*dt
          enddo
        enddo
!
        do l=1,lm
          do jb=1,lnsh
            do ib=its_h2,ite_h2
              tbs(ib,jb,l,1)=tbs(ib,jb,l,1)+tbs(ib,jb,l,2)*dt
              qbs(ib,jb,l,1)=qbs(ib,jb,l,1)+qbs(ib,jb,l,2)*dt
              wbs(ib,jb,l,1)=wbs(ib,jb,l,1)+wbs(ib,jb,l,2)*dt
            enddo
          enddo
        enddo
      endif
!
      if(n_bdy)then
        do jb=1,lnsh
          do ib=its_h2,ite_h2
            pdbn(ib,jb,1)=pdbn(ib,jb,1)+pdbn(ib,jb,2)*dt
          enddo
        enddo
!
        do l=1,lm
          do jb=1,lnsh
            do ib=its_h2,ite_h2
              tbn(ib,jb,l,1)=tbn(ib,jb,l,1)+tbn(ib,jb,l,2)*dt
              qbn(ib,jb,l,1)=qbn(ib,jb,l,1)+qbn(ib,jb,l,2)*dt
              wbn(ib,jb,l,1)=wbn(ib,jb,l,1)+wbn(ib,jb,l,2)*dt
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------western and eastern boundary------------------------------
!-----------------------------------------------------------------------
      if(w_bdy)then
        do jb=jts_h2,jte_h2
          do ib=1,lnsh
            pdbw(ib,jb,1)=pdbw(ib,jb,1)+pdbw(ib,jb,2)*dt
          enddo
        enddo
!
        do l=1,lm
          do jb=jts_h2,jte_h2
            do ib=1,lnsh
              tbw(ib,jb,l,1)=tbw(ib,jb,l,1)+tbw(ib,jb,l,2)*dt
              qbw(ib,jb,l,1)=qbw(ib,jb,l,1)+qbw(ib,jb,l,2)*dt
              wbw(ib,jb,l,1)=wbw(ib,jb,l,1)+wbw(ib,jb,l,2)*dt
            enddo
          enddo
        enddo
      endif
!
      if(e_bdy)then
        do jb=jts_h2,jte_h2
          do ib=1,lnsh
            pdbe(ib,jb,1)=pdbe(ib,jb,1)+pdbe(ib,jb,2)*dt
          enddo
        enddo
!
        do l=1,lm
          do jb=jts_h2,jte_h2
            do ib=1,lnsh
              tbe(ib,jb,l,1)=tbe(ib,jb,l,1)+tbe(ib,jb,l,2)*dt
              qbe(ib,jb,l,1)=qbe(ib,jb,l,1)+qbe(ib,jb,l,2)*dt
              wbe(ib,jb,l,1)=wbe(ib,jb,l,1)+wbe(ib,jb,l,2)*dt
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------southern boundary-----------------------------------------
!-----------------------------------------------------------------------
      if(s_bdy)then
        do j=1,lnsh
          jb=j
          ihs=jb
          ihe=ide+1-jb
          do i=max(its_h2,ihs),min(ite_h2,ihe)
            ib=i
            pd(i,j)=pdbs(ib,jb,1)*wh(jb)+pd(i,j)*(1.-wh(jb))
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=1,lnsh
            jb=j
            ihs=jb
            ihe=ide+1-jb
            do i=max(its_h2,ihs),min(ite_h2,ihe)
              ib=i
              t(i,j,l)=tbs(ib,jb,l,1)*wh(jb)+t(i,j,l)*(1.-wh(jb))
              q(i,j,l)=qbs(ib,jb,l,1)*wq(jb)+q(i,j,l)*(1.-wq(jb))
              cw(i,j,l)=wbs(ib,jb,l,1)*wq(jb)+cw(i,j,l)*(1.-wq(jb))
              pint(i,j,l+1)=(pint(i,j,l) &
                            +(dsg2(l)*pd(i,j)+pdsg1(l)))*wh(jb) &
                           +pint(i,j,l+1)*(1.-wh(jb))
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------northern boundary-----------------------------------------
!-----------------------------------------------------------------------
      if(n_bdy)then
        do j=jde-lnsh+1,jde
          jb=j-jde+lnsh
          ihs=1-jb+lnsh
          ihe=ide+jb-lnsh
          do i=max(its_h2,ihs),min(ite_h2,ihe)
            ib=i
            pd(i,j)=pdbn(ib,jb,1)*wh(lnsh+1-jb)+pd(i,j)*(1.-wh(lnsh+1-jb))
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jde-lnsh+1,jde
            jb=j-jde+lnsh
            ihs=1-jb+lnsh
            ihe=ide+jb-lnsh
            do i=max(its_h2,ihs),min(ite_h2,ihe)
              ib=i
              t(i,j,l)=tbn(ib,jb,l,1)*wh(lnsh+1-jb) &
                      +t(i,j,l)*(1.-wh(lnsh+1-jb))
              q(i,j,l)=qbn(ib,jb,l,1)*wq(lnsh+1-jb) &
                      +q(i,j,l)*(1.-wq(lnsh+1-jb))
              cw(i,j,l)=wbn(ib,jb,l,1)*wq(lnsh+1-jb) &
                       +cw(i,j,l)*(1.-wq(lnsh+1-jb))
              pint(i,j,l+1)=(pint(i,j,l) &
                            +(dsg2(l)*pd(i,j)+pdsg1(l)))*wh(lnsh+1-jb) &                    
                           +pint(i,j,l+1)*(1.-wh(lnsh+1-jb))
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------western boundary------------------------------------------
!-----------------------------------------------------------------------
      if(w_bdy)then
        do i=1,lnsh
          ib=i
          jhs=1+ib
          jhe=jde-ib
          do j=max(jts_h2,jhs),min(jte_h2,jhe)
            jb=j
            pd(i,j)=pdbw(ib,jb,1)*wh(ib)+pd(i,j)*(1.-wh(ib))
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do i=1,lnsh
            ib=i
            jhs=1+ib
            jhe=jde-ib
            do j=max(jts_h2,jhs),min(jte_h2,jhe)
              jb=j
              t(i,j,l)=tbw(ib,jb,l,1)*wh(ib)+t(i,j,l)*(1.-wh(ib))
              q(i,j,l)=qbw(ib,jb,l,1)*wq(ib)+q(i,j,l)*(1.-wq(ib))
              cw(i,j,l)=wbw(ib,jb,l,1)*wq(ib)+cw(i,j,l)*(1.-wq(ib))
              pint(i,j,l+1)=(pint(i,j,l) &
                            +(dsg2(l)*pd(i,j)+pdsg1(l)))*wh(ib) &
                           +pint(i,j,l+1)*(1.-wh(ib))
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------eastern boundary------------------------------------------
!-----------------------------------------------------------------------
      if(e_bdy)then
        do i=ide+1-lnsh,ide
          ib=i-ide+lnsh
          jhs=2+lnsh-ib
          jhe=jde-lnsh-1+ib
          do j=max(jts_h2,jhs),min(jte_h2,jhe)
            jb=j
            pd(i,j)=pdbe(ib,jb,1)*wh(lnsh+1-ib)+pd(i,j)*(1.-wh(lnsh+1-ib))
          enddo
        enddo
!
        do l=1,lm
          do i=ide+1-lnsh,ide
            ib=i-ide+lnsh
            jhs=2+lnsh-ib
            jhe=jde-lnsh-1+ib
            do j=max(jts_h2,jhs),min(jte_h2,jhe)
              jb=j
              t(i,j,l)=tbe(ib,jb,l,1)*wh(lnsh+1-ib) &
                      +t(i,j,l)*(1.-wh(lnsh+1-ib))
              q(i,j,l)=qbe(ib,jb,l,1)*wq(lnsh+1-ib) &
                      +q(i,j,l)*(1.-wq(lnsh+1-ib))
              cw(i,j,l)=wbe(ib,jb,l,1)*wq(lnsh+1-ib) &
                       +cw(i,j,l)*(1.-wq(lnsh+1-ib))
              pint(i,j,l+1)=(pint(i,j,l) &
                            +(dsg2(l)*pd(i,j)+pdsg1(l)))*wh(lnsh+1-ib) &
                           +pint(i,j,l+1)*(1.-wh(lnsh+1-ib))
            enddo
          enddo
        enddo
      endif
!
!-----------------------------------------------------------------------
!***  The four corner points.
!-----------------------------------------------------------------------
!
      if(s_bdy.and.w_bdy)then
        pd(its+1,jts+1)=(pd(its+1,jts)+pd(its,jts+1) &
                        +pd(its+2,jts+1)+pd(its+1,jts+2))*0.25
        do l=1,lm
          t(its+1,jts+1,l)=(t(its+1,jts,l)+t(its,jts+1,l) &
                           +t(its+2,jts+1,l)+t(its+1,jts+2,l))*0.25
          q(its+1,jts+1,l)=(q(its+1,jts,l)+q(its,jts+1,l) &
                           +q(its+2,jts+1,l)+q(its+1,jts+2,l))*0.25
          cw(its+1,jts+1,l)=(cw(its+1,jts,l)+cw(its,jts+1,l) &
                           + cw(its+2,jts+1,l)+cw(its+1,jts+2,l))*0.25
        enddo
      endif
!
      if(s_bdy.and.e_bdy)then
        pd(ite-1,jts+1)=(pd(ite-1,jts)+pd(ite-2,jts+1) &
                        +pd(ite,jts+1)+pd(ite-1,jts+2))*0.25
        do l=1,lm
          t(ite-1,jts+1,l)=(t(ite-1,jts,l)+t(ite-2,jts+1,l) &
                           +t(ite,jts+1,l)+t(ite-1,jts+2,l))*0.25
          q(ite-1,jts+1,l)=(q(ite-1,jts,l)+q(ite-2,jts+1,l) &
                           +q(ite,jts+1,l)+q(ite-1,jts+2,l))*0.25
          cw(ite-1,jts+1,l)=(cw(ite-1,jts,l)+cw(ite-2,jts+1,l) &
                            +cw(ite,jts+1,l)+cw(ite-1,jts+2,l))*0.25
        enddo
      endif
!
      if(n_bdy.and.w_bdy)then
        pd(its+1,jte-1)=(pd(its+1,jte-2)+pd(its,jte-1) &
                        +pd(its+2,jte-1)+pd(its+1,jte))*0.25
        do l=1,lm
          t(its+1,jte-1,l)=(t(its+1,jte-2,l)+t(its,jte-1,l) &
                           +t(its+2,jte-1,l)+t(its+1,jte,l))*0.25
          q(its+1,jte-1,l)=(q(its+1,jte-2,l)+q(its,jte-1,l) &
                           +q(its+2,jte-1,l)+q(its+1,jte,l))*0.25
          cw(its+1,jte-1,l)=(cw(its+1,jte-2,l)+cw(its,jte-1,l) &
                            +cw(its+2,jte-1,l)+cw(its+1,jte,l))*0.25
        enddo
      endif
!
      if(n_bdy.and.e_bdy)then
        pd(ite-1,jte-1)=(pd(ite-1,jte-2)+pd(ite-2,jte-1) &
                        +pd(ite,jte-1)+pd(ite-1,jte))*0.25
        do l=1,lm
          t(ite-1,jte-1,l)=(t(ite-1,jte-2,l)+t(ite-2,jte-1,l) &
                           +t(ite,jte-1,l)+t(ite-1,jte,l))*0.25
          q(ite-1,jte-1,l)=(q(ite-1,jte-2,l)+q(ite-2,jte-1,l) &
                           +q(ite,jte-1,l)+q(ite-1,jte,l))*0.25
          cw(ite-1,jte-1,l)=(cw(ite-1,jte-2,l)+cw(ite-2,jte-1,l) &
                            +cw(ite,jte-1,l)+cw(ite-1,jte,l))*0.25
        enddo
      endif
!
!-----------------------------------------------------------------------
!
      if(s_bdy)then
        do j=jts,jts+1
          do i=its_h2,ite_h2
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jts,jts+1
            do i=its_h2,ite_h2
              pint(i,j,l+1)=pint(i,j,l)+(dsg2(l)*pd(i,j)+pdsg1(l))
            enddo
          enddo
        enddo
      endif
!
      if(n_bdy)then
        do j=jte-1,jte
          do i=its_h2,ite_h2
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jte-1,jte
            do i=its_h2,ite_h2
              pint(i,j,l+1)=pint(i,j,l)+(dsg2(l)*pd(i,j)+pdsg1(l))
            enddo
          enddo
        enddo
      endif
!
      if(w_bdy)then
        do j=jts_h2,jte_h2
          do i=its,its+1
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jts_h2,jte_h2
            do i=its,its+1
              pint(i,j,l+1)=pint(i,j,l)+(dsg2(l)*pd(i,j)+pdsg1(l))
            enddo
          enddo
        enddo
      endif
!
      if(e_bdy)then
        do j=jts_h2,jte_h2
          do i=ite-1,ite
            pint(i,j,1)=pt
          enddo
        enddo
!
        do l=1,lm
          do j=jts_h2,jte_h2
            do i=ite-1,ite
              pint(i,j,l+1)=pint(i,j,l)+(dsg2(l)*pd(i,j)+pdsg1(l))
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!
      if(nsmud>=1)then
!
        smooth: do ks=1,nsmud
!-----------------------------------------------------------------------
          if(s_bdy)then
            do j=jts+1,jts-1+lines
              do i=max(its_h2,ids+1),min(ite_h2,ide-1)
                pdr(i,j)=(pd(i,j-1)+pd(i-1,j) &
                         +pd(i+1,j)+pd(i,j+1))*w1+w2*pd(i,j)       
              enddo
            enddo
!
            do l=1,lm
              do j=jts+1,jts-1+lines
                do i=max(its_h2,ids+1),min(ite_h2,ide-1)
                  tr(i,j,l)=(t(i,j-1,l)+t(i-1,j,l) &
                            +t(i+1,j,l)+t(i,j+1,l))*w1 &
                            +w2*t(i,j,l)       
                  qr(i,j,l)=(q(i,j-1,l)+q(i-1,j,l) &
                            +q(i+1,j,l)+q(i,j+1,l))*w1 &
                            +w2*q(i,j,l)
                  pr(i,j,l)=(pint(i,j-1,l+1)+pint(i-1,j,l+1) &
                            +pint(i+1,j,l+1)+pint(i,j+1,l+1))*w1 &
                            +w2*pint(i,j,l+1)       
                enddo
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do j=jte-lines+1,jte-1
              do i=max(its_h2,ids+1),min(ite_h2,ide-1)
                pdr(i,j)=(pd(i,j-1)+pd(i-1,j) &
                         +pd(i+1,j)+pd(i,j+1))*w1+w2*pd(i,j)       
              enddo
            enddo
!
            do l=1,lm
              do j=jte-lines+1,jte-1
                do i=max(its_h2,ids+1),min(ite_h2,ide-1)
                  tr(i,j,l)=(t(i,j-1,l)+t(i-1,j,l) &
                            +t(i+1,j,l)+t(i,j+1,l))*w1 &
                            +w2*t(i,j,l)       
                  qr(i,j,l)=(q(i,j-1,l)+q(i-1,j,l) &
                            +q(i+1,j,l)+q(i,j+1,l))*w1 &
                            +w2*q(i,j,l)
                  pr(i,j,l)=(pint(i,j-1,l+1)+pint(i-1,j,l+1) &
                            +pint(i+1,j,l+1)+pint(i,j+1,l+1))*w1 &
                            +w2*pint(i,j,l+1)       
                enddo
              enddo
            enddo
          endif
!
          if(w_bdy)then
            do j=max(jts_h2,jds+lines),min(jte_h2,jde-lines)
              do i=its+1,its-1+lines
                pdr(i,j)=(pd(i,j-1)+pd(i-1,j) &
                         +pd(i+1,j)+pd(i,j+1))*w1+w2*pd(i,j)       
              enddo
            enddo
!
            do l=1,lm
              do j=max(jts,jds+lines),min(jte,jde-lines)
                do i=its+1,its-1+lines
                  tr(i,j,l)=(t(i,j-1,l)+t(i-1,j,l) &
                            +t(i+1,j,l)+t(i,j+1,l))*w1 &
                            +w2*t(i,j,l)       
                  qr(i,j,l)=(q(i,j-1,l)+q(i-1,j,l) &
                            +q(i+1,j,l)+q(i,j+1,l))*w1 &
                            +w2*q(i,j,l)
                  pr(i,j,l)=(pint(i,j-1,l+1)+pint(i-1,j,l+1) &
                            +pint(i+1,j,l+1)+pint(i,j+1,l+1))*w1 &
                            +w2*pint(i,j,l+1)       
                enddo
              enddo
            enddo
          endif
!
          if(e_bdy)then
            do j=max(jts_h2,jds+lines),min(jte_h2,jde-lines)
              do i=ite-lines+1,ite-1
                pdr(i,j)=(pd(i,j-1)+pd(i-1,j) &
                         +pd(i+1,j)+pd(i,j+1))*w1+w2*pd(i,j)       
              enddo
            enddo
!
            do l=1,lm
              do j=max(jts_h2,jds+lines),min(jte_h2,jde-lines)
                do i=ite-lines+1,ite-1
                  tr(i,j,l)=(t(i,j-1,l)+t(i-1,j,l) &
                            +t(i+1,j,l)+t(i,j+1,l))*w1 &
                            +w2*t(i,j,l)       
                  qr(i,j,l)=(q(i,j-1,l)+q(i-1,j,l) &
                            +q(i+1,j,l)+q(i,j+1,l))*w1 &
                            +w2*q(i,j,l)
                  pr(i,j,l)=(pint(i,j-1,l+1)+pint(i-1,j,l+1) &
                            +pint(i+1,j,l+1)+pint(i,j+1,l+1))*w1 &
                            +w2*pint(i,j,l+1)       
                enddo
              enddo
            enddo
          endif
!-----------------------------------------------------------------------
          if(s_bdy)then
            do j=jts+1,jts-1+lines
              do i=max(its_h2,ids+1),min(ite_h2,ide-1)
                pd(i,j)=pdr(i,j)
              enddo
            enddo
!
            do l=1,lm
              do j=jts+1,jts-1+lines
                do i=max(its_h2,ids+1),min(ite_h2,ide-1)
                  t(i,j,l)=tr(i,j,l)
                  q(i,j,l)=qr(i,j,l)
                  pint(i,j,l+1)=pr(i,j,l)
                enddo
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do j=jte-lines+1,jte-1
              do i=max(its_h2,ids+1),min(ite_h2,ide-1)
                pd(i,j)=pdr(i,j)
              enddo
            enddo
!
            do l=1,lm
              do j=jte-lines+1,jte-1
                do i=max(its_h2,ids+1),min(ite_h2,ide-1)
                  t(i,j,l)=tr(i,j,l)
                  q(i,j,l)=qr(i,j,l)
                  pint(i,j,l+1)=pr(i,j,l)
                enddo
              enddo
            enddo
          endif
!
          if(w_bdy)then
            do j=max(jts_h2,jds+lines),min(jte_h2,jde-lines)
              do i=its+1,its-1+lines
                pd(i,j)=pdr(i,j)
              enddo
            enddo
!
            do l=1,lm
              do j=max(jts_h2,jds+lines),min(jte_h2,jde-lines)
                do i=its+1,its-1+lines
                  t(i,j,l)=tr(i,j,l)
                  q(i,j,l)=qr(i,j,l)
                  pint(i,j,l+1)=pr(i,j,l)
                enddo
              enddo
            enddo
          endif
!
          if(e_bdy)then
            do j=max(jts_h2,jds+lines),min(jte_h2,jde-lines)
              do i=ite-lines+1,ite-1
                pd(i,j)=pdr(i,j)
              enddo
            enddo
!
            do l=1,lm
              do j=max(jts_h2,jds+lines),min(jte_h2,jde-lines)
                do i=ite-lines+1,ite-1
                  t(i,j,l)=tr(i,j,l)
                  q(i,j,l)=qr(i,j,l)
                  pint(i,j,l+1)=pr(i,j,l)
                enddo
              enddo
            enddo
          endif
!-----------------------------------------------------------------------
        enddo smooth
!
      endif
!-----------------------------------------------------------------------
!
                        endsubroutine bocoh
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        subroutine bocov &
(lm,lnsv &
,dt &
,ube,ubn,ubs,ubw,vbe,vbn,vbs,vbw &
,u,v)
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
implicit none
!-----------------------------------------------------------------------

!include 'kind.inc'
!-----------------------------------------------------------------------
real(kind=kfpt),parameter:: &
 wa=0.5 &                    ! weighting factor
,w1=wa*0.25 &                ! weighting factor
,w2=1.-wa                    ! weighting factor
!-----------------------------------------------------------------------
integer(kind=kint),intent(in):: &
 lm &                        ! total # of levels
,lnsv                        ! blending area width, v points

real(kind=kfpt),intent(in):: &
 dt                          ! dynamics time step

real(kind=kfpt),dimension(ims:ime,1:lnsv,1:lm,1:2),intent(inout):: &
 ubn &                       ! u wind component at northern boundary
,ubs &                       ! u wind component at southern boundary
,vbn &                       ! v wind component at northern boundary
,vbs                         ! v wind component at southern boundary

real(kind=kfpt),dimension(1:lnsv,jms:jme,1:lm,1:2),intent(inout):: &
 ube &                       ! u wind component at eastern boundary
,ubw &                       ! u wind component at western boundary
,vbe &                       ! v wind component at eastern boundary
,vbw                         ! v wind component at western boundary

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm),intent(inout):: &
 u &                         ! u wind component
,v                           ! v wind component
!-----------------------------------------------------------------------
!---local variables-----------------------------------------------------
!-----------------------------------------------------------------------
integer(kind=kint):: &
 i &                         ! index in x direction
,ib &                        ! index in x direction, boundary zone
,ive &                       ! ending index in x direction, boundaries
,ivs &                       ! starting index in x direction, boundaries
,j &                         ! index in y direction
,jb &                        ! index in y direction, boundary zone
,jve &                       ! ending index in x direction, boundaries
,jvs &                       ! starting index in x direction, boundaries
,k &                         ! boundary line counter
,ks &                        ! smoothing counter
,l &                         ! index in p direction
,lines &                     ! boundary smoothing area
,nsmud                       ! number of smoothing passes

real(kind=kfpt),dimension(1:lnsv):: &
 wv(lnsv)                    ! blending weighting function, wind

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm):: &
 ur &                        ! u wind component
,vr                          ! v wind component
!-----------------------------------------------------------------------
!***********************************************************************
!-----------------------------------------------------------------------
      lines=lnsv
      nsmud=lines-1
!-----------------------------------------------------------------------
      wv(1)=1.
      do k=2,lnsv
        wv(k)=exp(-(real(k)-1.)/(real(lnsv)-1.)*7.)
      enddo
!-----------------------------------------------------------------------
!-------------update boundary values at v points------------------------
!-----------------------------------------------------------------------
!-------------southern and northern boundaries--------------------------
!-----------------------------------------------------------------------
      if(s_bdy)then
        do l=1,lm
          do jb=1,lnsv
            do ib=its_h2,min(ite_h2,ide-1)
              ubs(ib,jb,l,1)=ubs(ib,jb,l,1)+ubs(ib,jb,l,2)*dt
              vbs(ib,jb,l,1)=vbs(ib,jb,l,1)+vbs(ib,jb,l,2)*dt
            enddo
          enddo
        enddo
      endif
!
      if(n_bdy)then
        do l=1,lm
          do jb=1,lnsv
            do ib=its_h2,min(ite_h2,ide-1)
              ubn(ib,jb,l,1)=ubn(ib,jb,l,1)+ubn(ib,jb,l,2)*dt
              vbn(ib,jb,l,1)=vbn(ib,jb,l,1)+vbn(ib,jb,l,2)*dt
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------western and eastern boundaries----------------------------
!-----------------------------------------------------------------------
      if(w_bdy)then
        do l=1,lm
          do jb=jts_h2,min(jte_h2,jde-1)
            do ib=1,lnsv
              ubw(ib,jb,l,1)=ubw(ib,jb,l,1)+ubw(ib,jb,l,2)*dt
              vbw(ib,jb,l,1)=vbw(ib,jb,l,1)+vbw(ib,jb,l,2)*dt
            enddo
          enddo
        enddo
      endif
!
      if(e_bdy)then
        do l=1,lm
          do jb=jts_h1,min(jte_h2,jde-1)
            do ib=1,lnsv
              ube(ib,jb,l,1)=ube(ib,jb,l,1)+ube(ib,jb,l,2)*dt
              vbe(ib,jb,l,1)=vbe(ib,jb,l,1)+vbe(ib,jb,l,2)*dt
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------southern boundary-----------------------------------------
!-----------------------------------------------------------------------
      if(s_bdy)then
        do l=1,lm
          do j=jts,jts-1+lnsv
            jb=j
            ivs=max(its_h1,jb)
            ive=min(ite_h1,ide-jb)
            do i=ivs,ive
              ib=i
              u(i,j,l)=ubs(ib,jb,l,1)*wv(jb)+u(i,j,l)*(1.-wv(jb))
              v(i,j,l)=vbs(ib,jb,l,1)*wv(jb)+v(i,j,l)*(1.-wv(jb))
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------northern boundary-----------------------------------------
!-----------------------------------------------------------------------
      if(n_bdy)then
        do l=1,lm
          do j=jte-lnsv,jte-1
            jb=j-jte+lnsv+1
            ivs=max(its_h1,lnsv-jb+1)
            ive=min(ite_h1,ide+jb-lnsv-1)
            do i=ivs,ive
              ib=i
              u(i,j,l)=ubn(ib,jb,l,1)*wv(lnsv+1-jb) &
                      +u(i,j,l)*(1.-wv(lnsv+1-jb))
              v(i,j,l)=vbn(ib,jb,l,1)*wv(lnsv+1-jb) &
                      +v(i,j,l)*(1.-wv(lnsv+1-jb))
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------western boundary------------------------------------------
!-----------------------------------------------------------------------
      if(w_bdy)then
        do l=1,lm
          do i=its,its-1+lnsv
            ib=i
            jvs=max(jts_h1,1+ib)
            jve=min(jte_h1,jde-1-ib)
            do j=jvs,jve
              jb=j
              u(i,j,l)=ubw(ib,jb,l,1)*wv(ib)+u(i,j,l)*(1.-wv(ib))
              v(i,j,l)=vbw(ib,jb,l,1)*wv(ib)+v(i,j,l)*(1.-wv(ib))
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
!-------------eastern boundary------------------------------------------
!-----------------------------------------------------------------------
      if(e_bdy)then
        do i=ite-lnsv,ite-1
          ib=i-ide+lnsv+1
          jvs=max(jts_h1,lnsv-ib+2)
          jve=min(jte_h1,jde-lnsv+ib-2)
          do j=jvs,jve
            jb=j
            do l=1,lm
              u(i,j,l)=ube(ib,jb,l,1)*wv(lnsv+1-ib) &
                      +u(i,j,l)*(1.-wv(lnsv+1-ib))
              v(i,j,l)=vbe(ib,jb,l,1)*wv(lnsv+1-ib) &
                      +v(i,j,l)*(1.-wv(lnsv+1-ib))
            enddo
          enddo
        enddo
      endif
!-----------------------------------------------------------------------
      if(nsmud>=1)then
!
        smooth: do ks=1,nsmud
!-----------------------------------------------------------------------
          if(s_bdy)then
            do l=1,lm
              do j=jts+1,jts-1+lines
                do i=max(its_h1,ids+1),min(ite_h1,ide-2)
                  ur(i,j,l)=(u(i,j-1,l)+u(i,j-1,l) &
                            +u(i+1,j,l)+u(i,j+1,l))*w1+w2*u(i,j,l)       
                  vr(i,j,l)=(v(i,j-1,l)+v(i,j-1,l) &
                            +v(i+1,j,l)+v(i,j+1,l))*w1+w2*v(i,j,l)       
                enddo
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do l=1,lm
              do j=jte-lines,jte-2
                do i=max(its_h1,ids+1),min(ite_h1,ide-2)
                  ur(i,j,l)=(u(i,j-1,l)+u(i,j-1,l) &
                            +u(i+1,j,l)+u(i,j+1,l))*w1+w2*u(i,j,l)       
                  vr(i,j,l)=(v(i,j-1,l)+v(i,j-1,l) &
                            +v(i+1,j,l)+v(i,j+1,l))*w1+w2*v(i,j,l)       
                enddo
              enddo
            enddo
          endif
!
          if(w_bdy)then
            do l=1,lm
              do j=max(jts_h1,jds+lines),min(jte_h1,jde-lines-1)
                do i=its+1,its-1+lines
                  ur(i,j,l)=(u(i,j-1,l)+u(i,j-1,l) &
                            +u(i+1,j,l)+u(i,j+1,l))*w1+w2*u(i,j,l)       
                  vr(i,j,l)=(v(i,j-1,l)+v(i,j-1,l) &
                            +v(i+1,j,l)+v(i,j+1,l))*w1+w2*v(i,j,l)       
                enddo
              enddo
            enddo
          endif
!
          if(e_bdy)then
            do l=1,lm
              do j=max(jts_h1,jds+lines),min(jte_h1,jde-lines-1)
                do i=ite-lines,ite-2
                  ur(i,j,l)=(u(i,j-1,l)+u(i,j-1,l) &
                            +u(i+1,j,l)+u(i,j+1,l))*w1+w2*u(i,j,l)       
                  vr(i,j,l)=(v(i,j-1,l)+v(i,j-1,l) &
                            +v(i+1,j,l)+v(i,j+1,l))*w1+w2*v(i,j,l)       
                enddo
              enddo
            enddo
          endif
!-----------------------------------------------------------------------
          if(s_bdy)then
            do l=1,lm
              do j=jts+1,jts-1+lines
                do i=max(its_h1,ids+1),min(ite_h1,ide-2)
                  u(i,j,l)=ur(i,j,l)       
                  v(i,j,l)=vr(i,j,l)       
                enddo
              enddo
            enddo
          endif
!
          if(n_bdy)then
            do l=1,lm
              do j=jte-lines,jte-2
                do i=max(its_h1,ids+1),min(ite_h1,ide-2)
                  u(i,j,l)=ur(i,j,l)       
                  v(i,j,l)=vr(i,j,l)       
                enddo
              enddo
            enddo
          endif
!
          if(w_bdy)then
            do l=1,lm
              do j=max(jts_h1,jds+lines),min(jte_h1,jde-lines-1)
                do i=its+1,its-1+lines
                  u(i,j,l)=ur(i,j,l)       
                  v(i,j,l)=vr(i,j,l)       
                enddo
              enddo
            enddo
          endif
!
          if(e_bdy)then
            do l=1,lm
              do j=max(jts_h1,jds+lines),min(jte_h1,jde-lines-1)
                do i=ite-lines,ite-2
                  u(i,j,l)=ur(i,j,l)       
                  v(i,j,l)=vr(i,j,l)       
                enddo
              enddo
            enddo
          endif
        enddo smooth
!-----------------------------------------------------------------------
      endif
!-----------------------------------------------------------------------
!
                        endsubroutine bocov
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
                        endmodule module_fltbnds
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

