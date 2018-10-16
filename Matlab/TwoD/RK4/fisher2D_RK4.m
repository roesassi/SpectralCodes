function fisher2D_RK4(N,Nfinal,dt,ckeep,L,delta)
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
%-----------------Runge-Kutta----------------------------------
Eu=exp(-dt*ksq/2); Eu2=Eu.^2;
for n=1:Nfinal,
    k1u=dt*fft2(rhside(u,rhside_par{:}));
    u2=real(ifft2(Eu.*(uhat+k1u/2)));
    k2u=dt*fft2(rhside(u2,rhside_par{:}));
    u3=real(ifft2(Eu.*uhat+k2u/2));
    k3u=dt*fft2(rhside(u3,rhside_par{:}));
    u4=real(ifft2(Eu2.*uhat+Eu.*k3u));
    k4u=dt*fft2(rhside(u4,rhside_par{:}));
    uhat=Eu2.*uhat+(Eu2.*k1u+2*Eu.*(k2u+k3u)+k4u)/6;
    u=real(ifft2(uhat));
    if mod(n,ckeep)==0,
        ukeep(:,:,1+n/ckeep)=u;
    end
%    if mod(n,10)==0,
%        pause(0)
%        disp(int2str(n))
%    end
end
save('fisher2D_RK4.mat','tkeep','ukeep','N','L','x','y','Nfinal','ckeep')
%--------------Initial Condition ---------------------------------
function u=initial(x,y,delta)
[xx,yy]=meshgrid(x,y);
u=exp(-delta*(xx.^2+yy.^2))/2;
%---------------Right Hand Side-----------------------------------
function rhsu=rhside(u,delta)
rhsu=u.*(1-u);
