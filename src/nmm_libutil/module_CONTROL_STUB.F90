
                        module module_control
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
!      public NMMB_Finalize
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!---look-up tables------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine consts &
      (global,idtad,idtadt,secdif &
      ,smag2,smag4,codamp,wcor &
      ,pt &
      ,tph0d,tlm0d &
      ,dphd,dlmd &
      ,dxh,rdxh &
      ,dxv,rdxv &
      ,dyh,rdyh &
      ,dyv,rdyv &
      ,ddv,rddv &
      ,ddmpu,ddmpv &
      ,ef4t,wpdar &
      ,fcp,fdiv &
      ,curv,f &
      ,fad,fah &
      ,dare,rare &
      ,glat,glon &
      ,vlat,vlon &
      ,hdacx,hdacy &
      ,hdacvx,hdacvy &
      ,lnsh,lnsad &
      ,nboco,tboco)
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
integer(kind=kint),intent(in) :: &
 lnsh

integer(kind=kint),intent(out) :: &
 idtad &
,idtadt &
,lnsad &
,nboco

real(kind=kfpt),intent(in) :: &
 codamp &    ! divergence damping coefficient
,pt &        ! Pressure at top of domain (Pa)
,smag2 &     ! Smagorinsky coefficient for 2nd order diffusion 
,smag4 &     ! Smagorinsky coefficient for 4th order diffusion
,tlm0d &
,tph0d &
,wcor

logical(kind=klog),intent(in) :: &
 global
 
real(kind=kfpt),intent(out) :: &
 ddmpv &
,dlmd &
,dphd &
,dyh &
,dyv &
,ef4t &
,rdyh &
,rdyv &
,tboco
 
real(kind=kfpt),dimension(jds:jde),intent(out) :: &
 curv &
,dare &
,ddmpu &
,ddv &
,dxh &
,dxv &
,fad &
,fah &
,fcp &
,fdiv &
,rare &
,rddv &
,rdxh &
,rdxv &
,wpdar
 
real(kind=kfpt),dimension(ims:ime,jms:jme),intent(out) :: &
 f &
,glat &
,glon &
,vlat &
,vlon &
,hdacx &
,hdacy &
,hdacvx &
,hdacvy
!
logical(kind=klog),intent(in) :: &
 secdif
!
!-----------------------------------------------------------------------
!
                        end subroutine consts       
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine init &
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
      ,t,q,cw,psgdt &
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
       )
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
,indx_q2

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

real(kind=kfpt),dimension(ims:ime,jms:jme,1:lm-1),intent(out) :: &
 psgdt

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
!
                        end subroutine init
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine exptbl
!     ******************************************************************
!     *                                                                *
!     *               exponential function table                       *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
                        end subroutine exptbl
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        function zjexp(x)
!     ******************************************************************
!     *                                                                *
!     *               exponential function table                       *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
real(kind=kfpt):: &
 zjexp

!-----------------------------------------------------------------------
!--local variables------------------------------------------------------
!-----------------------------------------------------------------------

real(kind=kfpt):: &
,x                           ! argument
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        end function zjexp
!
!-----------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
                        subroutine tablep
!     ******************************************************************
!     *                                                                *
!     *    generates the table for finding pressure from               *
!     *    saturation specific humidity and potential temperature      *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
                        end subroutine tablep
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
                        subroutine tablet
!     ******************************************************************
!     *                                                                *
!     *    generates the table for finding temperature from            *
!     *    pressure and equivalent potential temperature               *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!-----------------------------------------------------------------------
!
                        end subroutine tablet
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine spline(jtb,nold,xold,yold,y2,nnew,xnew,ynew,p,q)           
!     ******************************************************************
!     *                                                                *
!     *  this is a one-dimensional cubic spline fitting routine        *
!     *  programed for a small scalar machine.                         *
!     *                                                                *
!     *  programer: z. janjic, yugoslav fed. hydromet. inst., beograd  *
!     *                                                                *
!     *                                                                *
!     *                                                                *
!     *  nold - number of given values of the function.  must be ge 3. *
!     *  xold - locations of the points at which the values of the     *
!     *         function are given.  must be in ascending order.       *
!     *  yold - the given values of the function at the points xold.   *
!     *  y2   - the second derivatives at the points xold.  if natural *
!     *         spline is fitted y2(1)=0. and y2(nold)=0. must be      *
!     *         specified.                                             *
!     *  nnew - number of values of the function to be calculated.     *
!     *  xnew - locations of the points at which the values of the     *
!     *         function are calculated.  xnew(k) must be ge xold(1)   *
!     *         and le xold(nold).                                     *
!     *  ynew - the values of the function to be calculated.           *
!     *  p, q - auxiliary vectors of the length nold-2.                *
!     *                                                                *
!     ******************************************************************
!-----------------------------------------------------------------------
!
      integer,intent(in) :: jtb,nnew,nold
!
      real,dimension(jtb),intent(in) :: xnew,xold,yold
      real,dimension(jtb),intent(inout) :: p,q,y2
      real,dimension(jtb),intent(out) :: ynew
!
!-----------------------------------------------------------------------
!
                        endsubroutine spline    
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine psitbl
!     ******************************************************************
!     *                                                                *
!     *               surface layer integral functions                 *
!     *               responsible person: z.janjic                     *
!     *                                                                *
!     ******************************************************************
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-----------------------------------------------------------------------
!
      implicit none
!
!
                        endsubroutine psitbl
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      FUNCTION TIMEF()
!
!-----------------------------------------------------------------------
!
      REAL*8 TIMEF
!
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!
      END FUNCTION TIMEF
!
!-----------------------------------------------------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-----------------------------------------------------------------------
!
      subroutine NMMB_Finalize
      integer irc
      stop 911
!
      end subroutine NMMB_Finalize
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
!
                       end module module_control
!
!-----------------------------------------------------------------------
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!-----------------------------------------------------------------------
