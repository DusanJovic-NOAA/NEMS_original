      subroutine rayleigh_damp(im,ix,iy,km,a,b,c,u1,v1,dt,cp,
     &                         levr,pgr,prsl,prslrd0)
!
!   ********************************************************************
! ----->  i m p l e m e n t a t i o n    v e r s i o n   <----------
!
!          --- rayleigh friction with total energy conservation ---
!              ----------------     -----------------------
!
!------ friction coefficient is based on deldif ----
!----------------------------------------------------------------------c
!    use
!        routine is called from gbphys  (after call to gwdps)
!
!    purpose
!        using the gwd parameterizations of ps-glas and ph-
!        gfdl technique.  the time tendencies of u v are 
!        altered to include/mimic the effect of non-stationary 
!        gravity wave drag from convection, frontgenosis,
!        wind shear etc.  loss of kinetic energy form gwd drag
!        is converted into internal energy.   
!
!  input
!        a(iy,km)  non-lin tendency for v wind component
!        b(iy,km)  non-lin tendency for u wind component
!        c(iy,km)  non-lin tendency for temperature
!        u1(ix,km) zonal wind m/sec  at t0-dt
!        v1(ix,km) meridional wind m/sec at t0-dt
!        t1(ix,km) temperature deg k at t0-dt
!
!        dt           time step    secs
!        pgr(ix)      surface pressure (pa)               
!        prsl(ix,km)  pressure at middle of layer (pa)
!
!  output
!        a, b, c as augmented by tendency due to rayleigh friction
!   ********************************************************************
      use machine , only : kind_phys
      implicit none
!
      integer,intent(in)  :: im, ix, iy, km,levr
      real(kind=kind_phys),intent(in)    :: dt, cp, prslrd0
      real(kind=kind_phys),intent(in)    :: pgr(ix),prsl(ix,km)
      real(kind=kind_phys),intent(in)    :: u1(ix,km), v1(ix,km)
      real(kind=kind_phys),intent(inout) :: a(iy,km), b(iy,km), c(iy,km)

!--- local variables
      real(kind=kind_phys) rtrd(ix,km),slrd0
      real(kind=kind_phys), parameter :: cons1=1.0, cons2=2.0, half=0.5
      real(kind=kind_phys) dtaux, dtauy, wrk1, rtrd1, rfactrd
      real(kind=kind_phys) eng0, eng1, tem1, tem2, dti, hfbcpdt
      integer i, k
!
! change prslrd0(2mb) to pascal
!      prslrd0  = 200.
!-----initialize some arrays
      slrd0=prslrd0/100000.0
!
!     rtrd1 = 1./(5*86400) ! reciprocal of time scale per scale height
      rtrd1 = 1./(10*86400) ! reciprocal of time scale per scale height
                            ! above beginning sigma level for rayleigh damping
!
! pressure in pascal
      do i = 1,im
       do k=1,km
        if(prsl(i,k)/pgr(i) < slrd0) then
          wrk1 = log(slrd0*pgr(i)/prsl(i,k))
          if (k > levr) then
            rtrd(i,k) = rtrd1 * wrk1 * wrk1
          else
            rtrd(i,k) = rtrd1 * wrk1
          endif
        else
          rtrd(i,k) = 0
        endif
       enddo
      enddo

      dti = cons1 / dt
      hfbcpdt = half / (cp*dt)

      do k = 1,km
        do i = 1,im
          rfactrd = cons1 / (cons1+dt*rtrd(i,k)) - cons1
          dtaux  = u1(i,k) * rfactrd
          dtauy  = v1(i,k) * rfactrd
          eng0   = u1(i,k)*u1(i,k) + v1(i,k)*v1(i,k)
          tem1   = u1(i,k) + dtaux
          tem2   = v1(i,k) + dtauy
          eng1   = tem1*tem1 + tem2*tem2
          a(i,k) = a(i,k) + dtauy * dti
          b(i,k) = b(i,k) + dtaux * dti
          c(i,k) = c(i,k) + (eng0-eng1) * hfbcpdt
        enddo
      enddo


      return
      end
