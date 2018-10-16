function fisher1D_RK4(N,Nfinal,dt,ckeep,L,delta)
if nargin<6
    disp('Using default parameters');
    N=256;
    Nfinal=200;
    dt=0.1;
    ckeep=10;
    L=50;
    delta=0.25;
end

x=(2*L/N)*(-N/2:N/2-1)';
u=initial(x,delta);
uhat=fft(u);
ukeep=zeros(N,1+Nfinal/ckeep);
ukeep(:,1)=u;
tkeep=dt*[0:ckeep:Nfinal];

ksq=((pi/L)*[0:N/2 -N/2+1:-1]').^2;
%-----------------Runge-Kutta----------------------------------
E=exp(-dt*ksq/2); E2=E.^2;
for n = 1:Nfinal
    k1=dt*fft(rhside(u));
    u2=real(ifft(E.*(uhat+k1/2)));
    k2=dt*fft(rhside(u2));
    u3=real(ifft(E.*uhat+k2/2));
    k3=dt*fft(rhside(u3));
    u4=real(ifft(E2.*uhat+E.*k3));
    k4=dt*fft(rhside(u4));
    uhat=E2.*uhat+(E2.*k1+2*E.*(k2+k3)+k4)/6;
    u=real(ifft(uhat));
    if mod(n,ckeep)==0,
        ukeep(:,1+n/ckeep)=u;
    end
end
save('fisher1D_RK4.mat','tkeep','ukeep','N','L','x','dt')
%----------------------Figures---------------------------------
figure(1)
plot(x,ukeep,'r')
title('Fishers equation using integrating factor')
xlabel('x');ylabel('u');axis([-50 50 0 1])
figure(2)
mesh(tkeep,x,ukeep)
xlabel('t');ylabel('x');zlabel('z')
%--------------Initial Condition ---------------------------------
function u=initial(x,delta)
u=1./(2*cosh(delta*x));
%---------------Right Hand Side-----------------------------------
function rhs2=rhside(u)
rhs2=u.*(1-u);
