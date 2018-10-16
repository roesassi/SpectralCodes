function gray1D_ETDRK4(N,Nfinal,dt,ckeep,L,epsilon,a,b)
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

rhside_par={a,b};
x=(2*L/N)*(-N/2:N/2-1)';
uv=initial(x,L);
uvhat=fft(uv);
uvkeep=zeros(N,2,1+Nfinal/ckeep);
uvkeep(:,:,1)=uv;
tkeep=dt*[0:ckeep:Nfinal];

ksq=[((pi/L)*[0:N/2 -N/2+1:-1]').^2, epsilon*((pi/L)*[0:N/2 -N/2+1:-1]').^2];
%-----------------Coefficients ----------------------------------
E2=exp(-dt*ksq/2); E=exp(-dt*ksq);
M = 16; % half the no. of points used in complex contour integration
r = exp(1i*pi*((1:M)-0.5)/M); % roots of 1i: (-1).^(1/(2*M))*exp(1i*[0:pi/M:pi-pi/M])
for n=1:2
    LR = -dt*ksq(:,n*ones(M,1)) + r(ones(N,1),:);
    B21(:,n) = sum(real((exp(LR/2)-1)./LR ),2)/M;
    C1(:,n) = sum(real((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3),2)/M;
    C23(:,n) = 2*sum(real((2+LR+exp(LR).*(-2+LR))./LR.^3),2)/M;
    C4(:,n) = sum(real((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3),2)/M;
end,
% The Ai equivalent would have been A1=0, A2=1/2, A3=1/2, A4=1
disp('Coefficients preparation completed')
%-----------------ETD Runge-Kutta----------------------------------
for n=1:Nfinal,
    k1=dt*fft(rhside(uv,rhside_par{:}));
    uvhat2=E2.*uvhat + B21.*k1;
    uv2=real(ifft(uvhat2));
    k2=dt*fft(rhside(uv2,rhside_par{:}));
    uv3=real(ifft(E2.*uvhat + B21.*k2));
    k3=dt*fft(rhside(uv3,rhside_par{:}));
    uv4=real(ifft(E2.*uvhat2 + B21.*(2*k3-k1)));
    k4=dt*fft(rhside(uv4,rhside_par{:}));
    uvhat = E.*uvhat + k1.*C1 + (k2+k3).*C23 + k4.*C4;
    uv=real(ifft(uvhat));
    if mod(n,ckeep)==0,
        uvkeep(:,:,1+n/ckeep)=uv;
    end
    if mod(n,10)==0,
        pause(0)
        disp(int2str(n))
    end
end
save('gray1D_ETDRK4.mat','tkeep','uvkeep','N','L','x')
%----------------------Figures---------------------------------
figure(1)
plot(x,squeeze(uvkeep(:,1,1:100:end)),'g',x,squeeze(uvkeep(:,2,1:100:end)),'r')
title('Solution using integrating factor')
xlabel('x');ylabel('u,v');axis('tight')
figure(2)
mesh(tkeep,x,squeeze(uvkeep(:,1,:)));view([60,75]);
xlabel('t');ylabel('x');zlabel('z');title('Surface plot of u')
figure(3)
mesh(tkeep,x,squeeze(uvkeep(:,2,:)));view([60,75]);
xlabel('t');ylabel('x');zlabel('z');title('Surface plot of v')
%--------------Initial Condition ---------------------------------
function u=initial(x,L)
u=[1-0.5*(sin(pi*(x-L)/(2*L)).^100) 0.25*(sin(pi*(x-L)/(2*L)).^100)];
%---------------Right Hand Side-----------------------------------
function rhs2=rhside(u,a,b)
t1=u(:,1).*u(:,2).*u(:,2);
rhs2=[-t1+a*(1-u(:,1)) t1-b*u(:,2)];
