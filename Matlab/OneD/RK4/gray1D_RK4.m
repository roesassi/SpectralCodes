function gray1D_RK4(N,Nfinal,dt,ckeep,L,epsilon,a,b)
if nargin<8;
    disp('Using default parameters');
    N=512;
    Nfinal=10000;
    dt=0.2;
    ckeep=10;
    L=50;
    epsilon=0.01;
    a=9*epsilon;
    b=0.4*epsilon^(1/3);
end

x=(2*L/N)*(-N/2:N/2-1)';
u=initial(x,L);
uhat=fft(u);
ukeep=zeros(N,2,1+Nfinal/ckeep);
ukeep(:,:,1)=u;
tkeep=dt*[0:ckeep:Nfinal];

ksq=((pi/L)*[0:N/2 -N/2+1:-1]').^2;
%-----------------Runge-Kutta----------------------------------
E=[exp(-dt*ksq/2) exp(-epsilon*dt*ksq/2)]; E2=E.^2;
for n=1:Nfinal
    k1=dt*fft(rhside(u,a,b));
    u2=real(ifft(E.*(uhat+k1/2)));
    k2=dt*fft(rhside(u2,a,b));
    u3=real(ifft(E.*uhat+k2/2));
    k3=dt*fft(rhside(u3,a,b));
    u4=real(ifft(E2.*uhat+E.*k3));
    k4=dt*fft(rhside(u4,a,b));
    uhat=E2.*uhat+(E2.*k1+2*E.*(k2+k3)+k4)/6;
    u=real(ifft(uhat));
    if mod(n,ckeep)==0,
        ukeep(:,:,1+n/ckeep)=u;
    end
end
save('gray1D_RK4.mat','tkeep','ukeep','N','L','x')
%----------------------Figures---------------------------------
figure(1)
plot(x,squeeze(ukeep(:,1,1:100:end)),'g',x,squeeze(ukeep(:,2,1:100:end)),'r')
title('Solution using integrating factor')
xlabel('x');ylabel('u,v');axis('tight')
figure(2)
mesh(tkeep,x,squeeze(ukeep(:,1,:)));view([60,75]);
xlabel('t');ylabel('x');zlabel('z');title('Surface plot of u')
figure(3)
mesh(tkeep,x,squeeze(ukeep(:,2,:)));view([60,75]);
xlabel('t');ylabel('x');zlabel('z');title('Surface plot of v')
%--------------Initial Condition ---------------------------------
function u=initial(x,L)
u=[1-0.5*(sin(pi*(x-L)/(2*L)).^100) 0.25*(sin(pi*(x-L)/(2*L)).^100)];
%---------------Right Hand Side-----------------------------------
function rhs2=rhside(u,a,b)
t1=u(:,1).*u(:,2).*u(:,2);
rhs2=[-t1+a*(1-u(:,1)) t1-b*u(:,2)];
