function fisher1D_ETDRK4_B(N,Nfinal,dt,ckeep,L,delta)
if nargin<6
    disp('Using default parameters');
    N=256;
    Nfinal=200;
    dt=0.1;
    ckeep=10;
    L=50;
    delta=0.25;
end

rhside_par={}; 
x=(2*L/N)*(-N/2:N/2-1)';
u=initial(x,delta);
uhat=fft(u);
ukeep=zeros(N,1+Nfinal/ckeep);
ukeep(:,1)=u;
tkeep=dt*[0:ckeep:Nfinal];

ksq=((pi/L)*[0:N/2 -N/2+1:-1]').^2;
%-----------------Coefficients ----------------------------------
E2=exp(-dt*ksq/2); E=exp(-dt*ksq);
M = 16; % half the no. of points used in complex contour integration
r = exp(1i*pi*((1:M)-0.5)/M); % roots of 1i: (-1).^(1/(2*M))*exp(1i*[0:pi/M:pi-pi/M])
LR = -dt*ksq(:,ones(M,1)) + r(ones(N,1),:);
B21 = sum(real((exp(LR/2)-1)./LR),2)/M;
B31 = sum(real(((LR-4).*exp(LR/2)+LR+4)./(LR.^2)),2)/M;
B32 = sum(real((4*exp(LR/2)-2*LR-4)./(LR.^2)),2)/M;
B41 = sum(real(((LR-2).*exp(LR)+LR+2)./(LR.^2)),2)/M;
B43 = sum(real((2*exp(LR)-2*LR-2)./(LR.^2)),2)/M;
C1 = sum(real((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3),2)/M;
C23 = 2*sum(real((2+LR+exp(LR).*(-2+LR))./LR.^3),2)/M;
C4 = sum(real((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3),2)/M;
% The Ai equivalent would have been A1=0, A2=1/2, A3=1/2, A4=1
disp('Coefficients preparation completed')
%-----------------ETD Runge-Kutta----------------------------------
for n=1:Nfinal,
    k1=dt*fft(rhside(u,rhside_par{:}));
    u2=real(ifft(E2.*uhat + B21.*k1));
    k2=dt*fft(rhside(u2,rhside_par{:}));
    u3=real(ifft(E2.*uhat + B31.*k1 + B32.*k2));
    k3=dt*fft(rhside(u3,rhside_par{:}));
    u4=real(ifft(E.*uhat + B41.*k1 + B43.*k3));
    k4=dt*fft(rhside(u4,rhside_par{:}));
    uhat = E.*uhat + k1.*C1 + (k2+k3).*C23 + k4.*C4;
    u=real(ifft(uhat)); 
    if mod(n,ckeep)==0,
        ukeep(:,1+n/ckeep)=u;
    end
    if mod(n,10)==0,
        pause(0)
        disp(int2str(n))
    end
end
save('fisher1D_ETDRK4_B.mat','tkeep','ukeep','N','L','x','dt')
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
