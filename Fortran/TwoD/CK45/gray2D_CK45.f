c     Solves the Gray-Scott problem from the text
c     Adaptive time step using Cash-Karp scheme with Integrating Factor

      program reaction_diffusion
      implicit none

* PROGRAM LOCAL PARAMETER & VARIABLE 
c     n are the number of modes (a power of 2)
      integer n, mxneq 
      parameter (n=128, mxneq=2048)
      integer len, lenh, lenhmax,istep
      parameter(len=2*n, lenh=len/2, lenhmax=len/2,istep=6)
      integer nprint, k, kn, mm, khalfp, i, m, ncount,ij,it,count
      double precision length, dy, radius2, dt,time
      double precision yf(lenh*n), xf(lenh*n)
* SERVICE VARIABLES
      double precision ee2(n,2,n), u(n,2,n), 
     +     u1(n,2,n), u2(n,2,n), u3(n,2,n), u4(n,2,n), 
     +     u5(n,2,n), u6(n,2,n), ustar(n,2,n),
     +     a(n,2,n), b(n,2,n), c(n,2,n), d(n,2,n), e(n,2,n), f(n,2,n),
     +     fftfftu(n,2,n), fftfftu2(n,2,n), fftfftustar(n,2,n),
     +     fftfftu3(n,2,n), fftfftu4(n,2,n),uold(n,2,n),
     +     fftfftu5(n,2,n),fftfftu6(n,2,n), 
     +     fftffta(n,2,n), fftfftb(n,2,n), 
     +     fftfftc(n,2,n), fftfftd(n,2,n),
     +     fftffte(n,2,n), fftfftf(n,2,n)
      double precision ffta(n,2,n), fftat(n,2,n), fftfftat(n,2,n)
      double precision crk(istep),crkstar(istep),brk(istep,istep),dtnew
     + ,ark(istep),etol,dtmax,dtmin,uerror,verror,tend,tfinal,tprint
      integer ifail
* COMMON BLOCKS
      double precision trigy(2,mxneq)
      integer nfay, ifay(20)
      common /yfac/ trigy, nfay, ifay
      double precision biga, bigb, epsilon, delta,aa,bb
      common /parms/ biga, bigb, epsilon
      double precision wavey1(mxneq), wavey2(mxneq)
      common /wvz/ wavey1, wavey2
      double precision pi
      common /pival/ pi
      double precision asn60, root
      common /roots/ asn60, root
      double precision rhs1,rhs2

c     useful variables
      asn60=0.5d0*dsqrt(3d0)
      root=1d0/dsqrt(2d0)
      pi=4d0*datan(1d0)
      call rkweight(crk,crkstar,brk,ark)

c     tracks the maximal error and timestep
      open(40,file='error.dat',status='unknown')
c     tracks the failed timesteps
      open(41,file='failures.dat',status='unknown')
c     stores the data from 5th order scheme
      open(51,file='gray_CK45.dat',status='unknown')
c     stores the data from 4th order scheme
      open(52,file='graystar_CK45.dat',status='unknown')

c     gridsize, pde evaluated on -length<x,y<length
c     dt is initial time step, ncount is number of timesteps before data output
c     nprint is the number of the times data is to be output
c     tfinal is final time and this sets tprint,
c     etol is error tolerance, dtmin and dtmax
c     are min and max allowable timesteps.
      tfinal=500d0
c      write(*,*) 'tfinal?'
c      read*,tfinal
      length=40d0
      dt=0.1d0
      nprint=10
      ncount=1000000
      tprint=tfinal/dble(nprint)
      ifail=0
      etol=1d-4
      dtmin=1d-10
      dtmax=0.9d0
      
c     ++++++++++++++++++ problem dependent parameters +++++++++++++++++
      epsilon=0.01d0
      delta=epsilon
      aa=9d0
      bb=0.4d0
      biga=aa*epsilon
      bigb=bb*(epsilon**(1d0/3d0))

      dy=2.0d0*pi/dble(n)
c     +++++++++++++++++++ initial condition ++++++++++++++++++++++++++
      do i=1,lenh
         xf(i)=((i-1)*dy)
         xf(i)=length*(xf(i)-pi)/pi
         do k=1,n
            yf(k)=((k-1)*dy)
            yf(k)=length*(yf(k)-pi)/pi
            radius2=yf(k)**2d0+xf(i)**2d0
         u1(i,1,k)=1d0-0.5d0*dexp(-0.05d0*radius2)
         u1(i,2,k)=0.25d0*dexp(-0.05d0*radius2)
         enddo
      enddo

c++++++++++++++++fft derivative initialization +++++++++++++++++
      call prefft(n,nfay,ifay,trigy)
      kn=n
      khalfp=n/2+1
      do i=1,n
         mm=i/khalfp
         m=mm*n+1
         wavey1(i)=(i-m)*pi/length
         wavey2(i)=wavey1(i)*wavey1(i)
      enddo
      wavey1(khalfp)=0.0

c+++++++++++++++++++the time stepping routine ++++++++++++++++++

      time=0d0
      uerror=0d0
      verror=0d0
c     Cash-Karp adaptive Runge-Kutta
      do 90 it=1,nprint

c     +++++++++++++++++++ a loop taken until print out ++++++++++
         do count=1,ncount

          tend=dble(it)*tprint
          if(dt.gt.tend-time) dt=tend-time

      do k=1,n
         do ij=1,lenh
            uold(ij,1,k)=u1(ij,1,k)
            uold(ij,2,k)=u1(ij,2,k)
         enddo
      enddo


 5        do k=1,n
               do ij=1,lenh
            ee2(ij,1,k)=dexp(-(wavey2(k)+wavey2(ij))*dt)
            ee2(ij,2,k)=dexp(-delta*(wavey2(k)+wavey2(ij))*dt)
                  a(ij,1,k)=dt*rhs1(u1(ij,1,k),u1(ij,2,k))
                  a(ij,2,k)=dt*rhs2(u1(ij,1,k),u1(ij,2,k))
               enddo
            enddo

            call doubletransform(a,fftffta,n,len,ffta,fftat,fftfftat)
            call doubletransform(u1,fftfftu,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
      fftfftu2(ij,1,k)=fftfftu(ij,1,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     +*ark(2))+brk(2,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*ark(2))
     +*fftffta(ij,1,k)
      fftfftu2(ij,2,k)=fftfftu(ij,2,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     + *ark(2)*delta)+brk(2,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(2)))*fftffta(ij,2,k)
               enddo
            enddo

            call inversetransform(fftfftu2,u2,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
                  b(ij,1,k)=dt*rhs1(u2(ij,1,k),u2(ij,2,k))
                  b(ij,2,k)=dt*rhs2(u2(ij,1,k),u2(ij,2,k))
               enddo
            enddo

            call doubletransform(b,fftfftb,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
      fftfftu3(ij,1,k)=fftfftu(ij,1,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     +*ark(3))+brk(3,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*ark(3))
     +*fftffta(ij,1,k)
     ++brk(3,2)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(3)-ark(2)))
     +*fftfftb(ij,1,k)
      fftfftu3(ij,2,k)=fftfftu(ij,2,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     + *ark(3)*delta)+brk(3,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*ark(3))*fftffta(ij,2,k)
     + +brk(3,2)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(3)-ark(2)))*fftfftb(ij,2,k)
               enddo
            enddo

            call inversetransform(fftfftu3,u3,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
                  c(ij,1,k)=dt*rhs1(u3(ij,1,k),u3(ij,2,k))
                  c(ij,2,k)=dt*rhs2(u3(ij,1,k),u3(ij,2,k))
               enddo
            enddo

            call doubletransform(c,fftfftc,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
      fftfftu4(ij,1,k)=fftfftu(ij,1,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     +*ark(4))+brk(4,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*ark(4))
     +*fftffta(ij,1,k)
     + +brk(4,2)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(4)-ark(2)))
     +*fftfftb(ij,1,k)
     ++brk(4,3)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(4)-ark(3)))
     +*fftfftc(ij,1,k)
      fftfftu4(ij,2,k)=fftfftu(ij,2,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     + *ark(4)*delta)+brk(4,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(4)))*fftffta(ij,2,k)
     + +brk(4,2)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(4)-ark(2)))*fftfftb(ij,2,k)
     + +brk(4,3)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(4)-ark(3)))*fftfftc(ij,2,k)
               enddo
            enddo

            call inversetransform(fftfftu4,u4,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
                  d(ij,1,k)=dt*rhs1(u4(ij,1,k),u4(ij,2,k))
                  d(ij,2,k)=dt*rhs2(u4(ij,1,k),u4(ij,2,k))
               enddo
            enddo

            call doubletransform(d,fftfftd,n,len,ffta,fftat,fftfftat)


            do k=1,n
               do ij=1,lenh
      fftfftu5(ij,1,k)=fftfftu(ij,1,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     +*ark(5))+brk(5,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*ark(5))
     +*fftffta(ij,1,k)
     + +brk(5,2)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(5)-ark(2)))
     +*fftfftb(ij,1,k)
     ++brk(5,3)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(5)-ark(3)))
     +*fftfftc(ij,1,k)
     ++brk(5,4)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(5)-ark(4)))
     +*fftfftd(ij,1,k)
      fftfftu5(ij,2,k)=fftfftu(ij,2,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     + *ark(5)*delta)+brk(5,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*ark(5))*fftffta(ij,2,k)
     + +brk(5,2)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(5)-ark(2)))*fftfftb(ij,2,k)
     + +brk(5,3)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(5)-ark(3)))*fftfftc(ij,2,k)
     + +brk(5,4)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(5)-ark(4)))*fftfftd(ij,2,k)
               enddo
            enddo

            call inversetransform(fftfftu5,u5,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
                  e(ij,1,k)=dt*rhs1(u5(ij,1,k),u5(ij,2,k))
                  e(ij,2,k)=dt*rhs2(u5(ij,1,k),u5(ij,2,k))
               enddo
            enddo

            call doubletransform(e,fftffte,n,len,ffta,fftat,fftfftat)


            do k=1,n
               do ij=1,lenh
      fftfftu6(ij,1,k)=fftfftu(ij,1,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     +*ark(6))+brk(6,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*ark(6))
     +*fftffta(ij,1,k)
     + +brk(6,2)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(6)-ark(2)))
     +*fftfftb(ij,1,k)
     ++brk(6,3)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(6)-ark(3)))
     +*fftfftc(ij,1,k)
     ++brk(6,4)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(6)-ark(4)))
     +*fftfftd(ij,1,k)
     ++brk(6,5)*dexp(-(wavey2(k)+wavey2(ij))*dt*(ark(6)-ark(5)))
     +*fftffte(ij,1,k)
      fftfftu6(ij,2,k)=fftfftu(ij,2,k)*dexp(-(wavey2(k)+wavey2(ij))*dt
     + *ark(6)*delta)+brk(6,1)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*ark(6))*fftffta(ij,2,k)
     + +brk(6,2)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(6)-ark(2)))*fftfftb(ij,2,k)
     + +brk(6,3)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(6)-ark(3)))*fftfftc(ij,2,k)
     + +brk(6,4)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(6)-ark(4)))*fftfftd(ij,2,k)
     + +brk(6,5)*dexp(-(wavey2(k)+wavey2(ij))*dt*
     +delta*(ark(6)-ark(5)))*fftffte(ij,2,k)
               enddo
            enddo

            call inversetransform(fftfftu6,u6,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
                  f(ij,1,k)=dt*rhs1(u6(ij,1,k),u6(ij,2,k))
                  f(ij,2,k)=dt*rhs2(u6(ij,1,k),u6(ij,2,k))
               enddo
            enddo

            call doubletransform(f,fftfftf,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
        fftfftustar(ij,1,k)=ee2(ij,1,k)*(fftfftu(ij,1,k)+
     +  crkstar(1)*fftffta(ij,1,k)
     ++crkstar(3)*fftfftc(ij,1,k)*dexp((wavey2(k)+wavey2(ij))*dt*ark(3))
     ++crkstar(4)*fftfftd(ij,1,k)*dexp((wavey2(k)+wavey2(ij))*dt*ark(4))
     ++crkstar(5)*fftffte(ij,1,k)*dexp((wavey2(k)+wavey2(ij))*dt*ark(5))
     ++crkstar(6)*fftfftf(ij,1,k)*dexp((wavey2(k)+wavey2(ij))*dt*ark(6))
     +)
        fftfftustar(ij,2,k)=ee2(ij,2,k)*(fftfftu(ij,2,k)+
     + crkstar(1)*fftffta(ij,2,k)
     ++crkstar(3)*fftfftc(ij,2,k)*
     +dexp(delta*(wavey2(k)+wavey2(ij))*dt*ark(3))
     ++crkstar(4)*fftfftd(ij,2,k)*
     +dexp(delta*(wavey2(k)+wavey2(ij))*dt*ark(4))
     ++crkstar(5)*fftffte(ij,2,k)*
     +dexp(delta*(wavey2(k)+wavey2(ij))*dt*ark(5))
     ++crkstar(6)*fftfftf(ij,2,k)*
     +dexp(delta*(wavey2(k)+wavey2(ij))*dt*ark(6)))

        fftfftu(ij,1,k)=ee2(ij,1,k)*(fftfftu(ij,1,k)+
     +  crk(1)*fftffta(ij,1,k)
     +  +crk(3)*fftfftc(ij,1,k)*dexp((wavey2(k)+wavey2(ij))*dt*ark(3))
     +  +crk(4)*fftfftd(ij,1,k)*dexp((wavey2(k)+wavey2(ij))*dt*ark(4))
     +  +crk(6)*fftfftf(ij,1,k)*dexp((wavey2(k)+wavey2(ij))*dt*ark(6)))
       fftfftu(ij,2,k)=ee2(ij,2,k)*(fftfftu(ij,2,k)+
     + crk(1)*fftffta(ij,2,k)+crk(3)*fftfftc(ij,2,k)*
     +dexp(delta*(wavey2(k)+wavey2(ij))*dt*ark(3))
     ++crk(4)*fftfftd(ij,2,k)*
     +dexp(delta*(wavey2(k)+wavey2(ij))*dt*ark(4))
     ++crk(6)*fftfftf(ij,2,k)*
     +dexp(delta*(wavey2(k)+wavey2(ij))*dt*ark(6)))
               enddo
            enddo

      call inversetransform(fftfftu,u,n,len,ffta,fftat,fftfftat)
      call inversetransform(fftfftustar,ustar,n,len,ffta,fftat,fftfftat)

         do k=1,n
         do ij=1,lenh
         if(dabs(u(ij,1,k)-ustar(ij,1,k)).gt.uerror) then
             uerror=dabs(u(ij,1,k)-ustar(ij,1,k))
         endif
         if(dabs(u(ij,2,k)-ustar(ij,2,k)).gt.verror) then
             verror=dabs(u(ij,2,k)-ustar(ij,2,k))
         endif
         enddo
         enddo

c failed time steps
          if(dabs(max(uerror,verror)).gt.etol) then
          dtnew=0.95d0*dt*(dabs(etol/max(uerror,verror)))**0.25d0
          uerror=0d0
          verror=0d0
c          write(*,*) 'Failed time step ', time,dt,dtnew
          ifail=ifail+1
          write(41,*) time,dt,dtnew,uerror,verror,ifail
          if((dtnew.lt.dtmin).and.(tend-time.gt.dtmin)) then
          write(*,*) 'Too small...'
            stop
          endif
         dt=min(dtnew,tend-time)
        do k=1,n
          do ij=1,lenh
             u1(ij,1,k)=uold(ij,1,k)
             u1(ij,2,k)=uold(ij,2,k)
       enddo
        enddo

         goto 5
         endif

c successful time step
        do k=1,n
          do ij=1,lenh
             u1(ij,1,k)=u(ij,1,k)
             u1(ij,2,k)=u(ij,2,k)
        enddo
        enddo

       time=time+dt
       dtnew=0.95d0*dt*(dabs(etol/max(uerror,verror)))**0.2d0

       write(40,*) time,dt,max(uerror,verror),dtnew
c      write(*,*) time,max(uerror,verror),dtnew,ifail

      if(dabs(time-tend).lt.1e-14) then
         time=tend
      write(*,*) 'Time=',time,' Failed steps=',ifail

c     printing out data
         do k=1,n
            do ij=1,lenh
               write(51,*) u1(ij,1,k),u1(ij,2,k)
               write(52,*) ustar(ij,1,k),ustar(ij,2,k)
            enddo
         enddo

        ifail=0
        goto 90
        endif

c        write(*,*) dtnew,dtmin,tend,time,tend-time
          if((dtnew.lt.dtmin).and.(tend-time.gt.dtmin)) then
             write(*,*) 'Required time step too small - terminating ..'
             stop
          endif
          dt=min(dtnew,dtmax,tend-time)

         uerror=0d0
         verror=0d0

         enddo
c     +++++++++++++++++ end of time loop  ++++++++++++++++++++++++

 90   continue


      end


c     +++++++++++++++++ the nonlinear terms in the pde +++++++++++++
      double precision function rhs1(u,v)
      implicit none
* PROCEDURE FORMAL ARGUMENTS
      double precision u, v
* COMMON BLOCKS
      double precision biga,bigb, epsilon
      common /parms/ biga, bigb, epsilon
      rhs1=-u*v*v+biga*(1d0-u)
      return
      end

      double precision function rhs2(u,v)
      implicit none
* PROCEDURE FORMAL ARGUMENTS
      double precision u, v
* COMMON BLOCKS
      double precision biga,bigb, epsilon
      common /parms/ biga, bigb, epsilon
      rhs2=u*v*v-bigb*v
      return
      end

c     ++++++++++++++++ useful slave routines +++++++++++++++++++++

      subroutine transpose(u,v,n)
      implicit none
* SUBROUTINE FORMAL ARGUMENTS
      integer n
      double precision u(n,2,n), v(n,2,n)
* SUBROUTINE LOCAL VARIABLE      
      integer k,ij
      do k=1,n
         do ij=1,n
            v(ij,1,k)=u(k,1,ij)
            v(ij,2,k)=u(k,2,ij)
         enddo
      enddo
      return
      end

      subroutine doubletransform(a,fftffta,n,len,ffta,fftat,fftfftat)
      implicit none
* SUBROUTINE FORMAL ARGUMENTS
      integer n, len
      double precision a(n,2,n), fftffta(n,2,n)
      double precision ffta(n,2,n), fftat(n,2,n), fftfftat(n,2,n)
* COMMON BLOCK
      integer mxneq
      parameter(mxneq=2048)
      double precision trigy(2,mxneq)
      integer nfay, ifay(20)
      common /yfac/ trigy, nfay, ifay
      
      external transpose
      call fft1(a,ffta,n,nfay,ifay,-1,trigy,len/2)
      call transpose(ffta,fftat,n)
      call fft1(fftat,fftfftat,n,nfay,ifay,-1,trigy,len/2)
      call transpose(fftfftat,fftffta,n)

      return
      end

      subroutine inversetransform(fftffta,a,n,len,ffta,fftat,fftfftat)
      implicit none
* SUBROUTINE FORMAL ARGUMENTS
      integer n, len
      double precision a(n,2,n), fftffta(n,2,n)
      double precision ffta(n,2,n), fftat(n,2,n), fftfftat(n,2,n)
* COMMON BLOCK
      integer mxneq
      parameter(mxneq=2048)
      double precision trigy(2,mxneq)
      integer nfay, ifay(20)
      common /yfac/ trigy, nfay, ifay

      external transpose
      call transpose(fftffta,fftfftat,n)
      call fft1(fftfftat,fftat,n,nfay,ifay,+1,trigy,len/2)
      call transpose(fftat,ffta,n)
      call fft1(ffta,a,n,nfay,ifay,+1,trigy,len/2)

      return
      end


c     +++++++++++++++++++FFT routines ++++++++++++++++++++++++++++++
c     These are primarily from Canuto et al 1988 Spectral Methods in
c     Fluid Mechanics, Springer-Verlag

      subroutine fft1(a,c,n,nfax,ifax,isign,trig,len)
      implicit none
      integer mxneq,n,i,ifac,la,nfax,isign,len,ij
      parameter(mxneq=2048)
      double precision a(n,2,0:n-1),c(n,2,0:n-1)
      double precision trig(2,0:mxneq-1),pi,xni
      integer ifax(*)
      logical odd
      common/pival/pi

      la=1
      odd=.true.
      do 10 i=1,nfax
         ifac=ifax(i)
         if(odd)then
            call pass1(a,c,n,isign,ifac,la,trig,len)
         else
            call pass1(c,a,n,isign,ifac,la,trig,len)
         endif
         odd=.not. odd
         la=la*ifac
 10   continue
      if(odd)then
         do 30 i=0,n-1
            do 20 ij=1,len
               c(ij,1,i)=a(ij,1,i)
               c(ij,2,i)=a(ij,2,i)
 20         continue
 30      continue
      endif
      if(isign.eq.-1) then
         xni=1./n
         do 50 i=0,n-1
            do 40 ij=1,len
               c(ij,1,i)=xni*c(ij,1,i)
               c(ij,2,i)=xni*c(ij,2,i)
 40         continue
 50      continue
      endif
      return
      end

      subroutine pass1(a,c,n,isign,ifac,la,trig,len)
      implicit none 
      integer mxneq,n
      parameter(mxneq=2048)
      double precision a(n,2,0:n-1),c(n,2,0:n-1),cc,ss,sn60,root,asn60
     +     ,trig(2,0:mxneq-1),t1,t2,s1,s2,c1,c2,ta1,ta2,ap1,ap2,am1,am2
      integer ind(0:20),jnd(0:20),j2,i2,ij,j0,j1,i0,i1,l,jump,i,k,j,m
     +     ,ifac,la,isign,len
      common/roots/asn60,root
      sn60=dble(isign)*asn60
      m=n/ifac

      do 10 k=0,ifac-1
         ind(k)=k*m
         jnd(k)=k*la
 10   continue

      i=0
      j=0
      jump=(ifac-1)*la
      do 130 k=0,m-la,la
         do 120 l=1,la
            if(ifac.eq.2)then
               i0=ind(0)+i
               i1=ind(1)+i
               j0=jnd(0)+j
               j1=jnd(1)+j
               cc=trig(1,k)
               ss=isign*trig(2,k)
               if(k.eq.0) then
                  do 20 ij=1,len
                     c(ij,1,j0)=a(ij,1,i0)+a(ij,1,i1)
                     c(ij,2,j0)=a(ij,2,i0)+a(ij,2,i1)
                     c(ij,1,j1)=a(ij,1,i0)-a(ij,1,i1)
                     c(ij,2,j1)=a(ij,2,i0)-a(ij,2,i1)
 20               continue
               else
                  do 50 ij=1,len
                     c(ij,1,j0)=a(ij,1,i0)+a(ij,1,i1)
                     c(ij,2,j0)=a(ij,2,i0)+a(ij,2,i1)
                     am1=a(ij,1,i0)-a(ij,1,i1)
                     am2=a(ij,2,i0)-a(ij,2,i1)
                     c(ij,1,j1)=cc*am1-ss*am2
                     c(ij,2,j1)=ss*am1+cc*am2
 50               continue
               endif
            elseif(ifac.eq.3) then
               i0=ind(0)+i
               i1=ind(1)+i
               i2=ind(2)+i
               j0=jnd(0)+j
               j1=jnd(1)+j
               j2=jnd(2)+j
               if(k.eq.0)then
                  do 60 ij=1,len
                     ap1=a(ij,1,i1)+a(ij,1,i2)
                     ap2=a(ij,2,i1)+a(ij,2,i2)
                     c(ij,1,j0)=a(ij,1,i0)+ap1
                     c(ij,2,j0)=a(ij,2,i0)+ap2
                     ta1=a(ij,1,i0)-0.5*ap1
                     ta2=a(ij,2,i0)-0.5*ap2
                     am1=sn60*(a(ij,1,i1)-a(ij,1,i2))
                     am2=sn60*(a(ij,2,i1)-a(ij,2,i2))
                     c(ij,1,j1)=ta1-am2
                     c(ij,2,j1)=ta2+am1
                     c(ij,1,j2)=ta1+am2
                     c(ij,2,j2)=ta2-am1
 60               continue
               else
                  c1=trig(1,k)
                  c2=trig(1,2*k)
                  s1=isign*trig(2,k)
                  s2=isign*trig(2,2*k)
                  do 70 ij=1,len
                     ap1=a(ij,1,i1)+a(ij,1,i2)
                     ap2=a(ij,2,i1)+a(ij,2,i2)
                     c(ij,1,j0)=a(ij,1,i0)+ap1
                     c(ij,2,j0)=a(ij,2,i0)+ap2
                     ta1=a(ij,1,i0)-0.5*ap1
                     ta2=a(ij,2,i0)-0.5*ap2
                     am1=sn60*(a(ij,1,i1)-a(ij,1,i2))
                     am2=sn60*(a(ij,2,i1)-a(ij,2,i2))
                     t1=ta1-am2
                     t2=ta2+am1
                     c(ij,1,j1)=c1*t1-s1*t2
                     c(ij,2,j1)=s1*t1+c1*t2
                     t1=ta1+am2
                     t2=ta2-am1
                     c(ij,1,j2)=c2*t1-s2*t2
                     c(ij,2,j2)=s2*t1+c2*t2
 70               continue
               endif
            endif
            i=i+1
            j=j+1
 120     continue
         j=j+jump
 130  continue
      return
      end


      subroutine prefft(n,nfax,ifax,trig)
      implicit none
      integer mxneq,k,n,nfax
      parameter(mxneq=2048)
      double precision trig(2,0:mxneq-1),arg,pi
      integer ifax(*)
      common/pival/pi
      call factor(n,nfax,ifax)
      do k=0,n-1
         arg=2d0*pi*dble(k)/dble(n)
         trig(1,k)=dcos(arg)
         trig(2,k)=dsin(arg)
      enddo
      return
      end

      subroutine factor(n,nfax,ifax)
      implicit none
      integer ifax(*),ii,n,nfax,nn
      nfax=0
      nn=n

      do 10 ii=1,20
         if(nn.eq.3*(nn/3)) then
            nfax=nfax+1
            ifax(nfax)=3
            nn=nn/3
         else
            goto 20
         endif
 10   continue
 20   continue

      do 30 ii=nfax+1,20
         if(nn.eq.2*(nn/2)) then
            nfax=nfax+1
            ifax(nfax)=2
            nn=nn/2
         else
            goto 40
         endif
 30   continue
 40   continue
      if(nn.ne.1)then
         stop
      endif
      return
      end


      subroutine rkweight(crk,crkstar,brk,ark)
      implicit double precision (a-h,o-z)
      double precision crk(6),crkstar(6),brk(6,6),ark(6)

      ark(2)=1d0/5d0
      ark(3)=3d0/10d0
      ark(4)=3d0/5d0
      ark(5)=1d0
      ark(6)=7d0/8d0
      crk(1)=37d0/378d0
      crk(2)=0d0
      crk(3)=250d0/621d0
      crk(4)=125d0/594d0
      crk(5)=0d0
      crk(6)=512d0/1771d0
      crkstar(1)=2825d0/27648d0
      crkstar(2)=0d0
      crkstar(3)=18575d0/48384d0
      crkstar(4)=13525d0/55296d0
      crkstar(5)=277d0/14336d0
      crkstar(6)=1d0/4d0
      brk(2,1)=1d0/5d0
      brk(3,1)=3d0/40d0
      brk(3,2)=9d0/40d0
      brk(4,1)=3d0/10d0
      brk(4,2)=-9d0/10d0
      brk(4,3)=6d0/5d0
      brk(5,1)=-11d0/54d0
      brk(5,2)=5d0/2d0
      brk(5,3)=-70d0/27d0
      brk(5,4)=35d0/27d0
      brk(6,1)=1631d0/55296d0
      brk(6,2)=175d0/512d0
      brk(6,3)=575d0/13824d0
      brk(6,4)=44275d0/110592d0
      brk(6,5)=253d0/4096d0
      return
      end
