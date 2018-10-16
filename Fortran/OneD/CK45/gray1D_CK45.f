c Solves the Gray-Scott 1D problem from the text
c Cash-Karp adaptive scheme
      program reaction_diffusion
      implicit none

* PROGRAM LOCAL PARAMETER & VARIABLE 
c     n are the number of modes (a power of 2)
      integer n,mxneq,istep
      parameter(mxneq=2048,n=256,istep=6)
      integer nprint,it,count,m,khalfp,mm,kn,i,k,ncount
      double precision length,time,dt,dy
* SERVICE VARIABLES
      double precision u(2,n),yf(mxneq),ee(n),ee2(n),u1(2,n),u2(2,n)
     .,u3(2,n),u4(2,n),u5(2,n),u6(2,n),fftu(2,n),fftu2(2,n),fftu3(2,n),
     .fftu4(2,n),fftu6(2,n), fftustar(2,n), ustar(2,n),
     .fftu5(2,n), a(2,n),b(2,n),c(2,n),d(2,n),e(2,n),f(2,n),ffta(2,n),
     .fftb(2,n),fftc(2,n),uold(2,n),fftd(2,n),ffte(2,n),fftf(2,n)
     .,eev(n),eev2(n)
      double precision crk(istep),crkstar(istep),brk(istep,istep),dtnew
     + ,ark(istep),etol,dtmax,dtmin,uerror,verror,tend,tfinal,tprint
      integer ifail
* COMMON BLOCKS
      double precision trigy(2,mxneq)
      integer nfay, ifay(20)
      common/yfac/trigy,nfay,ifay
      double precision wavey1(mxneq), wavey2(mxneq)
      common/wvz/wavey1,wavey2
      common/space/yf
      double precision asn60,root
      common/roots/asn60,root
      double precision pi
      common/pival/pi
      double precision biga,bigb,epsilon,delta,delta2
      common/parms/biga,bigb,epsilon


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
      open(51,file='gray1D_CK45.dat',status='unknown')
c     stores the data from 4th order scheme
      open(52,file='gray1Dstar_CK45.dat',status='unknown')

      tfinal=1000d0
       dy=2d0*pi/dble(n)
       dt=0.1d0
       length=50d0
       epsilon=0.01d0
       biga=9d0
       bigb=0.4d0
      delta2=epsilon
      delta=dsqrt(delta2)
      biga=biga*delta2
      bigb=bigb*(delta**(2d0/3d0))
        ncount=1000000
       nprint=10
      tprint=tfinal/dble(nprint)
      ifail=0
      etol=1d-5
      dtmin=1d-10
      dtmax=0.9d0

c +++++++++++++++++++initial condition++++++++++++++++++++++++++
      do k=1,n
            yf(k)=(dble(k-1)*dy)
            yf(k)=length*(yf(k)-pi)/pi
       u1(1,k)=1d0-0.5d0*(dsin(pi*(yf(k)-length)/(2d0*length)))**100d0
       u1(2,k)=0.25d0*(dsin(pi*(yf(k)-length)/(2d0*length)))**100d0
      enddo

c++++++++++++++++fft derivative initialization +++++++++++++++++
      call prefft(n,nfay,ifay,trigy)
      kn=n
      khalfp=n/2+1
      do i=1,n
         mm=i/khalfp
         m=mm*n+1
         wavey1(i)=dble(i-m)*pi/length
         wavey2(i)=wavey1(i)*wavey1(i)
      enddo
         wavey1(khalfp)=0.0d0

c+++++++++++++++++++the time stepping routine ++++++++++++++++++

      time=0d0
      uerror=0d0
      verror=0d0

       do 90 it=1,nprint
c     +++++++++++++++++++ a loop taken until print out ++++++++++
         do count=1,ncount 

          tend=dble(it)*tprint
          if(dt.gt.tend-time) dt=tend-time

      do k=1,n
            uold(1,k)=u1(1,k)
            uold(2,k)=u1(2,k)
      enddo

 5    do k=1,n
         ee(k)=dexp(-wavey2(k)*dt/2d0)
         ee2(k)=dexp(-wavey2(k)*dt)
         eev(k)=dexp(-epsilon*wavey2(k)*dt/2d0)
         eev2(k)=dexp(-epsilon*wavey2(k)*dt)
           call righthandside(a(1,k),a(2,k),u1(1,k),u1(2,k),dt)
      enddo

       call fft1(a,ffta,n,nfay,ifay,-1,trigy)
       call fft1(u1,fftu,n,nfay,ifay,-1,trigy)

       do k=1,n
      fftu2(1,k)=fftu(1,k)*dexp(-wavey2(k)*dt
     +*ark(2))+brk(2,1)*dexp(-(wavey2(k))*dt*ark(2))
     +*ffta(1,k)
      fftu2(2,k)=fftu(2,k)*dexp(-(wavey2(k))*dt
     + *ark(2)*epsilon)+brk(2,1)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(2)))*ffta(2,k)
       enddo
       call fft1(fftu2,u2,n,nfay,ifay,+1,trigy)

       do k=1,n
          call righthandside(b(1,k),b(2,k),u2(1,k),u2(2,k),dt)
       enddo
       call fft1(b,fftb,n,nfay,ifay,-1,trigy)

       do k=1,n
      fftu3(1,k)=fftu(1,k)*dexp(-(wavey2(k))*dt
     +*ark(3))+brk(3,1)*dexp(-(wavey2(k))*dt*ark(3))
     +*ffta(1,k)
     ++brk(3,2)*dexp(-(wavey2(k))*dt*(ark(3)-ark(2)))
     +*fftb(1,k)          
      fftu3(2,k)=fftu(2,k)*dexp(-(wavey2(k))*dt
     + *ark(3)*epsilon)+brk(3,1)*dexp(-(wavey2(k))*dt*
     +epsilon*ark(3))*ffta(2,k)+brk(3,2)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(3)-ark(2)))*fftb(2,k)
       enddo
       call fft1(fftu3,u3,n,nfay,ifay,+1,trigy)

       do k=1,n
          call righthandside(c(1,k),c(2,k),u3(1,k),u3(2,k),dt)
       enddo
       call fft1(c,fftc,n,nfay,ifay,-1,trigy)

       do k=1,n
      fftu4(1,k)=fftu(1,k)*dexp(-(wavey2(k))*dt
     +*ark(4))+brk(4,1)*dexp(-(wavey2(k))*dt*ark(4))
     +*ffta(1,k)
     + +brk(4,2)*dexp(-(wavey2(k))*dt*(ark(4)-ark(2)))
     +*fftb(1,k)
     ++brk(4,3)*dexp(-(wavey2(k))*dt*(ark(4)-ark(3)))
     +*fftc(1,k)
      fftu4(2,k)=fftu(2,k)*dexp(-(wavey2(k))*dt
     + *ark(4)*epsilon)+brk(4,1)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(4)))*ffta(2,k)
     + +brk(4,2)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(4)-ark(2)))*fftb(2,k)
     + +brk(4,3)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(4)-ark(3)))*fftc(2,k)
       enddo
       call fft1(fftu4,u4,n,nfay,ifay,+1,trigy)

          do k=1,n
          call righthandside(d(1,k),d(2,k),u4(1,k),u4(2,k),dt)
          enddo
       call fft1(d,fftd,n,nfay,ifay,-1,trigy)

        do k=1,n  
       fftu5(1,k)=fftu(1,k)*dexp(-(wavey2(k))*dt
     +*ark(5))+brk(5,1)*dexp(-(wavey2(k))*dt*ark(5))
     +*ffta(1,k)
     + +brk(5,2)*dexp(-(wavey2(k))*dt*(ark(5)-ark(2)))
     +*fftb(1,k)
     ++brk(5,3)*dexp(-(wavey2(k))*dt*(ark(5)-ark(3)))
     +*fftc(1,k)
     ++brk(5,4)*dexp(-(wavey2(k))*dt*(ark(5)-ark(4)))
     +*fftd(1,k)
       fftu5(2,k)=fftu(2,k)*dexp(-(wavey2(k))*dt
     + *ark(5)*epsilon)+brk(5,1)*dexp(-(wavey2(k))*dt*
     +epsilon*ark(5))*ffta(2,k)
     + +brk(5,2)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(5)-ark(2)))*fftb(2,k)
     + +brk(5,3)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(5)-ark(3)))*fftc(2,k)
     + +brk(5,4)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(5)-ark(4)))*fftd(2,k)
          enddo

       call fft1(fftu5,u5,n,nfay,ifay,+1,trigy)

          do k=1,n
          call righthandside(e(1,k),e(2,k),u5(1,k),u5(2,k),dt)
          enddo
       call fft1(e,ffte,n,nfay,ifay,-1,trigy)

       do k=1,n
      fftu6(1,k)=fftu(1,k)*dexp(-(wavey2(k))*dt
     +*ark(6))+brk(6,1)*dexp(-(wavey2(k))*dt*ark(6))*ffta(1,k)
     ++brk(6,2)*dexp(-(wavey2(k))*dt*(ark(6)-ark(2)))*fftb(1,k)
     ++brk(6,3)*dexp(-(wavey2(k))*dt*(ark(6)-ark(3)))*fftc(1,k)
     ++brk(6,4)*dexp(-(wavey2(k))*dt*(ark(6)-ark(4)))*fftd(1,k)
     ++brk(6,5)*dexp(-(wavey2(k))*dt*(ark(6)-ark(5)))*ffte(1,k)
      fftu6(2,k)=fftu(2,k)*dexp(-(wavey2(k))*dt
     + *ark(6)*epsilon)+brk(6,1)*dexp(-(wavey2(k))*dt*
     +epsilon*ark(6))*ffta(2,k)
     + +brk(6,2)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(6)-ark(2)))*fftb(2,k)
     + +brk(6,3)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(6)-ark(3)))*fftc(2,k)
     + +brk(6,4)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(6)-ark(4)))*fftd(2,k)
     + +brk(6,5)*dexp(-(wavey2(k))*dt*
     +epsilon*(ark(6)-ark(5)))*ffte(2,k)
      enddo

       call fft1(fftu6,u6,n,nfay,ifay,+1,trigy)

          do k=1,n
          call righthandside(f(1,k),f(2,k),u6(1,k),u6(2,k),dt)
          enddo
       call fft1(f,fftf,n,nfay,ifay,-1,trigy)

            do k=1,n
        fftustar(1,k)=ee2(k)*(fftu(1,k)+
     +  crkstar(1)*ffta(1,k)
     ++crkstar(3)*fftc(1,k)*dexp((wavey2(k))*dt*ark(3))
     ++crkstar(4)*fftd(1,k)*dexp((wavey2(k))*dt*ark(4))
     ++crkstar(5)*ffte(1,k)*dexp((wavey2(k))*dt*ark(5))
     ++crkstar(6)*fftf(1,k)*dexp((wavey2(k))*dt*ark(6)))
        fftustar(2,k)=eev2(k)*(fftu(2,k)+
     + crkstar(1)*ffta(2,k)
     ++crkstar(3)*fftc(2,k)*
     +dexp(epsilon*(wavey2(k))*dt*ark(3))
     ++crkstar(4)*fftd(2,k)*
     +dexp(epsilon*(wavey2(k))*dt*ark(4))
     ++crkstar(5)*ffte(2,k)*
     +dexp(epsilon*(wavey2(k))*dt*ark(5))
     ++crkstar(6)*fftf(2,k)*
     +dexp(epsilon*(wavey2(k))*dt*ark(6)))
        fftu(1,k)=ee2(k)*(fftu(1,k)+crk(1)*ffta(1,k)
     +  +crk(3)*fftc(1,k)*dexp((wavey2(k))*dt*ark(3))
     +  +crk(4)*fftd(1,k)*dexp((wavey2(k))*dt*ark(4))
     +  +crk(6)*fftf(1,k)*dexp((wavey2(k))*dt*ark(6)))
       fftu(2,k)=eev2(k)*(fftu(2,k)+
     + crk(1)*ffta(2,k)+crk(3)*fftc(2,k)*
     +dexp(epsilon*(wavey2(k))*dt*ark(3))
     ++crk(4)*fftd(2,k)*
     +dexp(epsilon*(wavey2(k))*dt*ark(4))
     ++crk(6)*fftf(2,k)*
     +dexp(epsilon*(wavey2(k))*dt*ark(6)))
        enddo

       call fft1(fftu,u,n,nfay,ifay,+1,trigy)
       call fft1(fftustar,ustar,n,nfay,ifay,+1,trigy)

        do k=1,n
         if(dabs(u(1,k)-ustar(1,k)).gt.uerror) then
             uerror=dabs(u(1,k)-ustar(1,k))
         endif
         if(dabs(u(2,k)-ustar(2,k)).gt.verror) then
             verror=dabs(u(2,k)-ustar(2,k))
         endif
         enddo

c failed time steps
          if(dabs(max(uerror,verror)).gt.etol) then
          dtnew=0.95d0*dt*(dabs(etol/max(uerror,verror)))**0.25d0
          ifail=ifail+1
          write(41,*) time,dt,dtnew,uerror,verror,ifail
          if((dtnew.lt.dtmin).and.(tend-time.gt.dtmin)) then
          write(*,*) 'Too small...'
            stop
          endif
          uerror=0d0
          verror=0d0
         dt=min(dtnew,tend-time)
        do k=1,n
              u1(1,k)=uold(1,k)
              u1(2,k)=uold(2,k)
          enddo

         goto 5
         endif

c successful time step
        do k=1,n
             u1(1,k)=u(1,k)
             u1(2,k)=u(2,k)
        enddo

       time=time+dt
       dtnew=0.95d0*dt*(dabs(etol/max(uerror,verror)))**0.2d0
       write(40,*) time,dt,max(uerror,verror),dtnew

      if(dabs(time-tend).lt.1e-14) then
         time=tend
      write(*,*) 'Time=',time,' Failed steps=',ifail

c     printing out data
         do k=1,n
               write(51,*) sngl(yf(k)),u1(1,k),u(2,k)
               write(52,*) yf(k),ustar(1,k),ustar(2,k)
         enddo

        ifail=0
        goto 90
        endif
          if((dtnew.lt.dtmin).and.(tend-time.gt.dtmin)) then
             write(*,*) 'Required time step too small - terminating ..'
             stop
          endif
          dt=min(dtnew,dtmax,tend-time)

         uerror=0d0
         verror=0d0
         enddo
c     +++++++++++++++++ end of time loop  ++++++++++++++++++++++++

 90      continue

      end

c     +++++++++++++++++ the nonlinear terms in the pde +++++++++++++
      subroutine righthandside(rhs1,rhs2,u,v,dt)
      implicit none
      double precision rhs1,rhs2,u,v,tmp,dt
* COMMON BLOCKS
      double precision biga,bigb,epsilon
      common/parms/biga,bigb,epsilon
      tmp=u*v*v
      rhs1=dt*(-tmp+biga*(1d0-u))
      rhs2=dt*(tmp-bigb*v)
      return
      end

c     +++++++++++++++++++FFT routines ++++++++++++++++++++++++++++++
c     These are primarily from Canuto et al 1988 Spectral Methods in
c     Fluid Mechanics, Springer-Verlag

      subroutine fft1(a,c,n,nfax,ifax,isign,trig)
      implicit none
      integer mxneq,n,i,ifac,la,nfax,isign,ij
      parameter(mxneq=2048)
      double precision a(2,0:n-1),c(2,0:n-1)
      double precision trig(2,0:mxneq-1),pi,xni
      integer ifax(*)
      logical odd
      common/pival/pi

      la=1
      odd=.true.
      do 10 i=1,nfax
         ifac=ifax(i)
         if(odd)then
            call pass1(a,c,n,isign,ifac,la,trig,1)
         else
            call pass1(c,a,n,isign,ifac,la,trig,1)
         endif
         odd=.not. odd
         la=la*ifac
 10   continue
      if(odd)then
         do 30 i=0,n-1
            do 20 ij=1,2
               c(ij,i)=a(ij,i)
 20         continue
 30      continue
      endif
      if(isign.eq.-1) then
         xni=1./n
         do 50 i=0,n-1
            do 40 ij=1,2
               c(ij,i)=xni*c(ij,i)
 40         continue
 50      continue
      endif
      return
      end

      subroutine pass1(a,c,n,isign,ifac,la,trig,len)
      implicit none 
      integer mxneq,n
      parameter(mxneq=2048)
      double precision a(1,2,0:n-1),c(1,2,0:n-1),cc,ss,sn60,root,asn60
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
