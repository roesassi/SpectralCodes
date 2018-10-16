function fisher2D_ETDRK4_B(N,Nfinal,dt,ckeep,L,delta)
if nargin<6
    disp('Using default parameters')
    Nfinal=200;
    dt=0.1;
    ckeep=10;
    N=128;
    L=20;
    delta=0.25;
end

rhside_par={delta};
x=(2*L/N)*(-N/2:N/2-1)'; y=x;
u=initial(x,y,rhside_par{:});
uhat=fft2(u);
ukeep=zeros(N,N,1+Nfinal/ckeep); ukeep(:,:,1)=u;
tkeep=dt*[0:ckeep:Nfinal];

kx=(pi/L)*[0:N/2 -N/2+1:-1]'; ky=kx;
[kxx,kyy]=meshgrid(kx,ky);
ksq=(kxx.^2+kyy.^2);
%-----------------Coefficients ----------------------------------
E2=exp(-dt*ksq/2); E=exp(-dt*ksq); %Eu2.^2;
M = 16; % half the no. of points used in complex contour integration
r = exp(1i*pi*((1:M)-.5)/M); % roots of 1i: (-1).^(1/(2*M))*exp(1i*[0:pi/M:pi-pi/M])
for n=1:N/2+1,
    LR = -dt*ksq(n:N/2+1,n*ones(M,1)) + r(ones(N/2+1-n+1,1),:);
    B21(n:N/2+1,n) = sum(real( (exp(LR/2)-1)./LR ),2)/M;
    B31(n:N/2+1,n) = sum(real(((LR-4).*exp(LR/2)+LR+4)./(LR.^2)),2)/M;
    B32(n:N/2+1,n) = sum(real((4*exp(LR/2)-2*LR-4)./(LR.^2)),2)/M;
    B41(n:N/2+1,n) = sum(real(((LR-2).*exp(LR)+LR+2)./(LR.^2)),2)/M;
    B43(n:N/2+1,n) = sum(real((2*exp(LR)-2*LR-2)./(LR.^2)),2)/M;
    C1(n:N/2+1,n) = sum(real((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3),2)/M;
    C23(n:N/2+1,n) = 2*sum(real((2+LR+exp(LR).*(-2+LR))./LR.^3),2)/M;
    C4(n:N/2+1,n) = sum(real((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3),2)/M;
end,
% The Ai equivalent would have been A1=0, A2=1/2, A3=1/2, A4=1
B21=rsym(B21); B31=rsym(B31); B32=rsym(B32); B41=rsym(B41); B43=rsym(B43);
C1=rsym(C1); C23=rsym(C23); C4=rsym(C4);
disp('Coefficients preparation completed')
%-----------------ETD Runge-Kutta----------------------------------
for n=1:Nfinal,
    k1=dt*fft2(rhside(u,rhside_par{:}));
    u2=real(ifft2(E2.*uhat + B21.*k1));
    k2=dt*fft2(rhside(u2,rhside_par{:}));
    u3=real(ifft2(E2.*uhat + B31.*k1 + B32.*k2));
    k3=dt*fft2(rhside(u3,rhside_par{:}));
    u4=real(ifft2(E.*uhat + B41.*k1 + B43.*k3));
    k4=dt*fft2(rhside(u4,rhside_par{:}));
    uhat = E.*uhat + k1.*C1 + (k2+k3).*C23 + k4.*C4;
    u=real(ifft2(uhat)); 
    if mod(n,ckeep)==0,
        ukeep(:,:,1+n/ckeep)=u;
    end
    if mod(n,10)==0,
        pause(0)
        disp(int2str(n))
    end
end
save('fisher2D_ETDRK4_B.mat','tkeep','ukeep','N','L','x','y','Nfinal','ckeep')
%--------------Initial Condition ---------------------------------
function u=initial(x,y,delta)
[xx,yy]=meshgrid(x,y);
u=exp(-delta*(xx.^2+yy.^2))/2;
%---------------Right Hand Side-----------------------------------
function rhsu=rhside(u,delta)
rhsu=u.*(1-u);
%---------------Service routine-----------------------------------
function Qout=rsym(Qin);
Qout=tril(Qin,-1)+tril(Qin,-1)'+diag(diag(Qin));
Qout=[Qout, fliplr(Qout(:,2:end-1))];
Qout=[Qout; flipud(Qout(2:end-1,:))];
