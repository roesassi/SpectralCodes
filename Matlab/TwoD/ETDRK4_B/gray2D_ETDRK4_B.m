function gray2D_ETDRK4_B(N,Nfinal,dt,ckeep,L,epsilon,a,b)
if nargin<8;
    disp('Using default parameters');
    N=256;
    Nfinal=5000;
    dt=0.1;
    ckeep=500;
    L=40;
    epsilon=0.01;
    a=9*epsilon;
    b=0.4*epsilon^(1/3);
end

rhside_par={a,b};
x=(2*L/N)*(-N/2:N/2-1)'; y=x;
[u,v]=initial(x,y,rhside_par{:});
uhat=fft2(u);
vhat=fft2(v);
ukeep=zeros(N,N,1+Nfinal/ckeep); ukeep(:,:,1)=u;
vkeep=zeros(N,N,1+Nfinal/ckeep); vkeep(:,:,1)=v;
tkeep=dt*[0:ckeep:Nfinal];

kx=(pi/L)*[0:N/2 -N/2+1:-1]'; ky=kx;
[kxx,kyy]=meshgrid(kx,ky);
ksq=(kxx.^2+kyy.^2);
%-----------------Coefficients ----------------------------------
Eu2=exp(-dt*ksq/2); Eu=exp(-dt*ksq); %Eu2.^2;
Ev2=exp(-epsilon*dt*ksq/2); Ev=exp(-epsilon*dt*ksq);
M = 16; % half the no. of points used in complex contour integration
r = exp(1i*pi*((1:M)-.5)/M); % roots of 1i: (-1).^(1/(2*M))*exp(1i*[0:pi/M:pi-pi/M])
for n=1:N/2+1,
    LR = -dt*ksq(n:N/2+1,n*ones(M,1)) + r(ones(N/2+1-n+1,1),:);
    B21u(n:N/2+1,n) = sum(real( (exp(LR/2)-1)./LR ),2)/M;
    B31u(n:N/2+1,n) = sum(real(((LR-4).*exp(LR/2)+LR+4)./(LR.^2)),2)/M;
    B32u(n:N/2+1,n) = sum(real((4*exp(LR/2)-2*LR-4)./(LR.^2)),2)/M;
    B41u(n:N/2+1,n) = sum(real(((LR-2).*exp(LR)+LR+2)./(LR.^2)),2)/M;
    B43u(n:N/2+1,n) = sum(real((2*exp(LR)-2*LR-2)./(LR.^2)),2)/M;
    C1u(n:N/2+1,n) = sum(real((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3),2)/M;
    C23u(n:N/2+1,n) = 2*sum(real((2+LR+exp(LR).*(-2+LR))./LR.^3),2)/M;
    C4u(n:N/2+1,n) = sum(real((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3),2)/M;
    LR = -dt*(epsilon*ksq(n:N/2+1,n*ones(M,1))) + r(ones(N/2+1-n+1,1),:);
    B21v(n:N/2+1,n) = sum(real((exp(LR/2)-1)./LR),2)/M;
    B31v(n:N/2+1,n) = sum(real(((LR-4).*exp(LR/2)+LR+4)./(LR.^2)),2)/M;
    B32v(n:N/2+1,n) = sum(real((4*exp(LR/2)-2*LR-4)./(LR.^2)),2)/M;
    B41v(n:N/2+1,n) = sum(real(((LR-2).*exp(LR)+LR+2)./(LR.^2)),2)/M;
    B43v(n:N/2+1,n) = sum(real((2*exp(LR)-2*LR-2)./(LR.^2)),2)/M;
    C1v(n:N/2+1,n) = sum(real((-4-LR+exp(LR).*(4-3*LR+LR.^2))./LR.^3),2)/M;
    C23v(n:N/2+1,n) = 2*sum(real((2+LR+exp(LR).*(-2+LR))./LR.^3),2)/M;
    C4v(n:N/2+1,n) = sum(real((-4-3*LR-LR.^2+exp(LR).*(4-LR))./LR.^3),2)/M;
end,
% The Ai equivalent would have been A1=0, A2=1/2, A3=1/2, A4=1
B21u=rsym(B21u); B31u=rsym(B31u); B32u=rsym(B32u); B41u=rsym(B41u); B43u=rsym(B43u);
C1u=rsym(C1u); C23u=rsym(C23u); C4u=rsym(C4u);
B21v=rsym(B21v); B31v=rsym(B31v); B32v=rsym(B32v); B41v=rsym(B41v); B43v=rsym(B43v);
C1v=rsym(C1v); C23v=rsym(C23v); C4v=rsym(C4v);
disp('Coefficients preparation completed')
%-----------------ETD Runge-Kutta----------------------------------
for n=1:Nfinal,
    [dui,dvi]=rhside(u,v,rhside_par{:});
    k1u=dt*fft2(dui); k1v=dt*fft2(dvi);
    u2=real(ifft2(Eu2.*uhat + B21u.*k1u));
    v2=real(ifft2(Ev2.*vhat + B21v.*k1v));
    [dui,dvi]=rhside(u2,v2,rhside_par{:});
    k2u=dt*fft2(dui); k2v=dt*fft2(dvi);
    u3=real(ifft2(Eu2.*uhat + B31u.*k1u + B32u.*k2u));
    v3=real(ifft2(Ev2.*vhat + B31v.*k1v + B32v.*k2v));
    [dui,dvi]=rhside(u3,v3,rhside_par{:});
    k3u=dt*fft2(dui); k3v=dt*fft2(dvi);
    u4=real(ifft2(Eu.*uhat + B41u.*k1u + B43u.*k3u));
    v4=real(ifft2(Ev.*vhat + B41v.*k1v + B43v.*k3v));
    [dui,dvi]=rhside(u4,v4,rhside_par{:});
    k4u=dt*fft2(dui); k4v=dt*fft2(dvi);
    uhat = Eu.*uhat + k1u.*C1u + (k2u+k3u).*C23u + k4u.*C4u;
    vhat = Ev.*vhat + k1v.*C1v + (k2v+k3v).*C23v + k4v.*C4v;
    u=real(ifft2(uhat)); v=real(ifft2(vhat));
    if mod(n,ckeep)==0,
        ukeep(:,:,1+n/ckeep)=u;
        vkeep(:,:,1+n/ckeep)=v;
    end
    if mod(n,10)==0,
        pause(0)
        disp(int2str(n))
    end
end
save('gray2D_ETDRK4_B.mat','tkeep','ukeep','vkeep','N','L','x','y','Nfinal','ckeep')
%--------------Initial Condition ---------------------------------
function [u,v]=initial(x,y,a,b);
[xx,yy]=meshgrid(x,y);
common_exp=exp(-0.05*(xx.^2+yy.^2));
u=1-0.5*common_exp;
v=0.25*common_exp;
%---------------Right Hand Side-----------------------------------
function [rhsu,rhsv]=rhside(u,v,a,b);
uv2=u.*(v.^2);
rhsu=-uv2+a*(1-u);
rhsv=uv2-b*v;
%---------------Service routine-----------------------------------
function Qout=rsym(Qin);
Qout=tril(Qin,-1)+tril(Qin,-1)'+diag(diag(Qin));
Qout=[Qout, fliplr(Qout(:,2:end-1))];
Qout=[Qout; flipud(Qout(2:end-1,:))];
