c     Solves the labyrinthine pattern problem from the text
c     Constant time step: 4th order Runge-Kutta with Integrating factor

      program reaction_diffusion
      implicit none

* PROGRAM LOCAL PARAMETER & VARIABLE 
c     n are the number of modes (a power of 2)
      integer n, mxneq 
      parameter (n=128, mxneq=2048)
      integer len, lenh, lenhmax
      parameter(len=2*n, lenh=len/2, lenhmax=len/2)
      integer nprint, k, kn, mm, khalfp, i, m, ncount,ij,it,count
      double precision length, uminus, dy, radius, dt,time
      double precision u(lenhmax,2,n), yf(lenh*n), xf(lenh*n)
* SERVICE VARIABLES
      double precision ee(n,2,n), ee2(n,2,n), 
     +     u1(n,2,n), u2(n,2,n), u3(n,2,n), u4(n,2,n), 
     +     a(n,2,n), b(n,2,n), c(n,2,n), d(n,2,n), 
     +     fftfftu(n,2,n), fftfftu2(n,2,n), 
     +     fftfftu3(n,2,n), fftfftu4(n,2,n), 
     +     fftffta(n,2,n), fftfftb(n,2,n), 
     +     fftfftc(n,2,n), fftfftd(n,2,n)
      double precision ffta(n,2,n), fftat(n,2,n), fftfftat(n,2,n)
* COMMON BLOCKS
      double precision trigy(2,mxneq)
      integer nfay, ifay(20)
      common /yfac/ trigy, nfay, ifay
      double precision a0, a1, epsilon, delta, vminus
      common /parms/ a0, a1, epsilon, delta, vminus
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

      open(51,file='laby_RK4.dat',status='unknown')

c     gridsize, pde evaluated on -length<x,y<length
c     dt is time step, ncount is number of timesteps before data output
c     nprint is the number of the times data is to be output
      length=100d0
      dt=0.2d0
      nprint=10
      ncount=100

c     ++++++++++++++++++ problem dependent parameters +++++++++++++++++
      a0=-0.1d0
      a1=2d0
      epsilon=0.05d0
      delta=4d0
      uminus=-dsqrt(-(1d0-a1)/a1)

      do i=1,1000
         uminus=uminus-((a1*uminus**3)+uminus*(1d0-a1)-a0)/
     +        (3d0*a1*uminus**2+1d0-a1)
      enddo
      vminus=(uminus-a0)/a1

      dy=2.0d0*pi/dble(n)
c     +++++++++++++++++++ initial condition ++++++++++++++++++++++++++
      do i=1,lenh
         xf(i)=((i-1)*dy)
         xf(i)=length*(xf(i)-pi)/pi
         do k=1,n
            yf(k)=((k-1)*dy)
            yf(k)=length*(yf(k)-pi)/pi
            radius=dsqrt(0.01d0*yf(k)**2+xf(i)**2)
            u(i,1,k)=a1*vminus+a0-
     +           2d0*a1*((2d0*vminus*dexp(-0.1*radius**2)))
            u(i,2,k)=vminus-(2d0*vminus*dexp(-0.1d0*radius**2))
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

      do k=1,n
         do ij=1,lenh
            u1(ij,1,k)=u(ij,1,k)
            u1(ij,2,k)=u(ij,2,k)
            ee(ij,1,k)=dexp(-(wavey2(k)+wavey2(ij))*dt/2d0)
            ee2(ij,1,k)=dexp(-(wavey2(k)+wavey2(ij))*dt)
            ee(ij,2,k)=dexp(-delta*(wavey2(k)+wavey2(ij))*dt/2d0)
            ee2(ij,2,k)=dexp(-delta*(wavey2(k)+wavey2(ij))*dt)
         enddo
      enddo

c     standard 4th order Runge-Kutta
      do 90 it=1,nprint

c     +++++++++++++++++++ a loop taken until print out ++++++++++
         do count=1,ncount 

            time=dble(it*count)*dt

            do k=1,n
               do ij=1,lenh
                  a(ij,1,k)=dt*rhs1(u1(ij,1,k),u1(ij,2,k))
                  a(ij,2,k)=dt*rhs2(u1(ij,1,k),u1(ij,2,k))
               enddo
            enddo

            call doubletransform(a,fftffta,n,len,ffta,fftat,fftfftat)
            call doubletransform(u1,fftfftu,n,len,ffta,fftat,fftfftat)

            do k=1,n
               do ij=1,lenh
                  fftfftu2(ij,1,k)=ee(ij,1,k)*(fftfftu(ij,1,k)+
     +                 fftffta(ij,1,k)/2d0)
                  fftfftu2(ij,2,k)=ee(ij,2,k)*(fftfftu(ij,2,k)+
     +                 fftffta(ij,2,k)/2d0)
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
                  fftfftu3(ij,1,k)=ee(ij,1,k)*fftfftu(ij,1,k)+
     +                 fftfftb(ij,1,k)/2d0
                  fftfftu3(ij,2,k)=ee(ij,2,k)*fftfftu(ij,2,k)+
     +                 fftfftb(ij,2,k)/2d0
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
                  fftfftu4(ij,1,k)=ee2(ij,1,k)*fftfftu(ij,1,k)+
     +                 ee(ij,1,k)*fftfftc(ij,1,k)
                  fftfftu4(ij,2,k)=ee2(ij,2,k)*fftfftu(ij,2,k)+
     +                 ee(ij,2,k)*fftfftc(ij,2,k)
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
                  fftfftu(ij,1,k)=ee2(ij,1,k)*fftfftu(ij,1,k)+
     +                 (ee2(ij,1,k)*fftffta(ij,1,k)
     .                 +2d0*ee(ij,1,k)*(fftfftb(ij,1,k)+
     +                 fftfftc(ij,1,k))+fftfftd(ij,1,k))/6d0
                  fftfftu(ij,2,k)=ee2(ij,2,k)*fftfftu(ij,2,k)+
     +                 (ee2(ij,2,k)*fftffta(ij,2,k)
     .                 +2d0*ee(ij,2,k)*(fftfftb(ij,2,k)+
     +                 fftfftc(ij,2,k))+fftfftd(ij,2,k))/6d0
               enddo
            enddo

            call inversetransform(fftfftu,u1,n,len,ffta,fftat,fftfftat)
         enddo
c     +++++++++++++++++ end of time loop  ++++++++++++++++++++++++

c     printing out data
         do k=1,n
            do ij=1,lenh
               write(51,*) u1(ij,1,k),u1(ij,2,k)
            enddo
         enddo
         write(*,*) 'Time= ',time

 90   continue


      end


c     +++++++++++++++++ the nonlinear terms in the pde +++++++++++++
      double precision function rhs1(u,v)
      implicit none
* PROCEDURE FORMAL ARGUMENTS
      double precision u, v
* COMMON BLOCKS
      double precision a0, a1, epsilon, delta, vminus
      common /parms/ a0, a1, epsilon, delta, vminus
      rhs1=u-u**3-v
      return
      end

      double precision function rhs2(u,v)
      implicit none
* PROCEDURE FORMAL ARGUMENTS
      double precision u, v
* COMMON BLOCKS
      double precision a0, a1, epsilon, delta, vminus
      common /parms/ a0, a1, epsilon, delta, vminus
      rhs2=epsilon*(u-a1*v-a0)
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

