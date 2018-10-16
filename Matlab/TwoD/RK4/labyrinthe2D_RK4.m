function labyrinthe2D_RK4(N,Nfinal,dt,ckeep,L,a0,a1,varepsilon,delta)
if nargin<9
    disp('Using default parameters')
    Nfinal=10000; %Nfinal=20000;
    dt=0.1;
    ckeep=1000;
    N=128; %N=64;
    L=100;
    a0=-0.1;
    a1=2;
    varepsilon=0.05;
    delta=4;
end

epsilon=delta;
rhside_par={a0,a1,varepsilon};

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
%-----------------Runge-Kutta----------------------------------
Eu=exp(-dt*ksq/2); Eu2=Eu.^2;
Ev=exp(-epsilon*dt*ksq/2); Ev2=Ev.^2;
for n=1:Nfinal,
    [dui,dvi]=rhside(u,v,rhside_par{:}); k1u=dt*fft2(dui); k1v=dt*fft2(dvi);
    u2=real(ifft2(Eu.*(uhat+k1u/2))); v2=real(ifft2(Ev.*(vhat+k1v/2)));
    [dui,dvi]=rhside(u2,v2,rhside_par{:}); k2u=dt*fft2(dui); k2v=dt*fft2(dvi);
    u3=real(ifft2(Eu.*uhat+k2u/2)); v3=real(ifft2(Ev.*vhat+k2v/2));
    [dui,dvi]=rhside(u3,v3,rhside_par{:}); k3u=dt*fft2(dui); k3v=dt*fft2(dvi);
    u4=real(ifft2(Eu2.*uhat+Eu.*k3u)); v4=real(ifft2(Ev2.*vhat+Ev.*k3v));
    [dui,dvi]=rhside(u4,v4,rhside_par{:}); k4u=dt*fft2(dui); k4v=dt*fft2(dvi);
    uhat=Eu2.*uhat+(Eu2.*k1u+2*Eu.*(k2u+k3u)+k4u)/6;
    vhat=Ev2.*vhat+(Ev2.*k1v+2*Ev.*(k2v+k3v)+k4v)/6;
    u=real(ifft2(uhat)); v=real(ifft2(vhat));
    if mod(n,ckeep)==0,
        ukeep(:,:,1+n/ckeep)=u;
        vkeep(:,:,1+n/ckeep)=v;
    end
%    if mod(n,100)==0,
%        pause(0)
%        disp(int2str(n))
%    end
end
save('labyrinthe2D_RK4.mat','tkeep','ukeep','vkeep','N','L','x','y','Nfinal','ckeep')
%--------------Initial Condition ---------------------------------
function [u,v]=initial(x,y,a0,a1,varepsilon)
u_minus=roots([a1 0 (1-a1) -a0]);
u_minus=min(u_minus(find(isreal(u_minus))));
v_minus=(u_minus-a0)/a1;
[xx,yy]=meshgrid(x,y);
common_exp=exp(-0.1*(xx.^2+0.01*yy.^2));
u = a1*v_minus+a0-4*a1*v_minus*common_exp;
v = v_minus-2*v_minus*common_exp;
%---------------Right Hand Side-----------------------------------
function [rhsu,rhsv]=rhside(u,v,a0,a1,varepsilon)
rhsu=u-u.^3-v; rhsv=varepsilon*(u-a1*v-a0);
