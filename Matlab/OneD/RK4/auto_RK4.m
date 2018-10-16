function auto_RK4(N,Nfinal,dt,ckeep,L,epsilon)if nargin<6    disp('Using default parameters')    N = 512;    Nfinal=25000;    dt=0.02;    ckeep=500;    L=120;    epsilon=0.1;endx=(2*L/N)*(-N/2:N/2-1)';u=initial(x);uhat=fft(u);ukeep=zeros(N,2,1+Nfinal/ckeep);ukeep(:,:,1)=u;tkeep=dt*[0:ckeep:Nfinal];ksq=((pi/L)*[0:N/2 -N/2+1:-1]').^2;%-----------------Runge-Kutta----------------------------------------E=[exp(-dt*ksq/2) exp(-epsilon*dt*ksq/2)]; E2=E.^2;for n = 1:Nfinal    k1=dt.*fft(rhside(u));    u2=real(ifft(E.*(uhat+k1/2)));    k2=dt.*fft(rhside(u2));    u3=real(ifft(E.*uhat+k2/2));    k3=dt.*fft(rhside(u3));    u4=real(ifft(E2.*uhat+E.*k3));    k4=dt.*fft(rhside(u4));    uhat=E2.*uhat+(E2.*k1+2*E.*(k2+k3)+k4)/6;    u=real(ifft(uhat));    if mod(n,ckeep)==0,        ukeep(:,:,1+n/ckeep)=u;    endendsave('auto_RK4.mat','tkeep','ukeep','N','L','x')%----------------------Figures---------------------------------------figure(1)plot(x,squeeze(ukeep(:,1,1:10:end)),'g',x,squeeze(ukeep(:,2,1:10:end)),'r')title('Solution using integrating factor')xlabel('x');ylabel('u,v');axis('tight')figure(2)mesh(tkeep,x,squeeze(ukeep(:,1,:)));view([60,75]);xlabel('t');ylabel('x');zlabel('z')title('Surface plot of u')figure(3)mesh(tkeep,x,squeeze(ukeep(:,2,:)));view([60,75]);xlabel('t');ylabel('x');zlabel('z')title('Surface plot of v')%--------------Initial Condition ------------------------------------function u0=initial(x);u0=[0.5*(1+tanh((10-abs(x))/0.1)) 1-0.25*(1+tanh((10-abs(x))/0.1))];%---------------Right Hand Side--------------------------------------function rhs2=rhside(u)fu=max(u(:,1),0).^11; t1=u(:,2).*fu;rhs2=[t1 -t1];