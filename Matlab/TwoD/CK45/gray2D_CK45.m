function gray2D_CK45(N,tfinal,h,hkeep,L,epsilon,a,b)
if nargin<8;
    disp('Using default parameters');
    N=256;
    tfinal=500;
    h=0.1;
    hkeep=50;
    L=40;
    epsilon=0.01;
    a=9*epsilon;
    b=0.4*epsilon^(1/3);
end

rhside_par={a,b};
x=(2*L/N)*(-N/2:N/2-1)'; y=x;
[u,v]=initial(x,y,rhside_par{:});
uvhat=[fft2(u), fft2(v)];
ukeep(:,:,1)=u; vkeep(:,:,1)=v;
tkeep=0; tcheck=0;

% Cash-Karp coefficients
rk_a=[0 1/5 3/10 3/5 1 7/8];
rk_b=[0 0 0 0 0; 1/5 0 0 0 0; 3/40 9/40 0 0 0; 3/10 -9/10 6/5 0 0; ...
        -11/54 5/2 -70/27 35/27 0; 1631/55296 175/512 575/13824 44275/110592 253/4096];
rk_c=[37/378 0 250/621 125/594 0 512/1771];
rk_d=[2825/27648 0 18575/48384 13525/55296 277/14336 1/4];

kx=(pi/L)*[0:N/2 -N/2+1:-1]'; ky=kx;
[kxx,kyy]=meshgrid(kx,ky);
ksq=(kxx.^2+kyy.^2);
hmax=0.9; h_old=0; t=0;
MINTIMESTEP=1e-10; EPSTOL=1e-4;
rejected_steps=0; accepted_steps=0;
while h >= MINTIMESTEP,
    if abs(h-h_old)>MINTIMESTEP,
        E_ai(:,:,1)=ones(N,2*N);
        for i=2:6, E_ai(:,:,i)=[exp(-ksq*rk_a(i)*h), exp(-epsilon*ksq*rk_a(i)*h)]; end,
        E=[exp(-ksq*h), exp(-epsilon*ksq*h)];
    end,

    % Embedded RK45 step
    [fuv,guv]=rhside(u,v,rhside_par{:});
    munutilde(:,:,1)=h*[fft2(fuv), fft2(guv)];
    deltauvhat=(rk_c(1)-rk_d(1))*munutilde(:,:,1);
    for i=2:6,
        sumbmunutilde=rk_b(i,1)*munutilde(:,:,1);
        for j=2:i-1, sumbmunutilde=sumbmunutilde+rk_b(i,j)*munutilde(:,:,j); end,
        uvhat_ai=E_ai(:,:,i).*(uvhat+sumbmunutilde);
        [fuv,guv]=rhside(real(ifft2(uvhat_ai(:,1:N))),real(ifft2(uvhat_ai(:,N+1:end))),rhside_par{:});
        munutilde(:,:,i)=h*[fft2(fuv), fft2(guv)]./E_ai(:,:,i);
        deltauvhat=deltauvhat+(rk_c(i)-rk_d(i))*munutilde(:,:,i);
    end,
    deltauvhat=E.*deltauvhat;
    deltauv=[real(ifft2(deltauvhat(:,1:N))), real(ifft2(deltauvhat(:,N+1:end)))];
    errmax = max(max(abs(deltauv)));

    h_old=h;
    if errmax>EPSTOL, % rejected step: the time step is reduced
        disp(['                       h = ' num2str(h) ' WAS REJECTED (errmax=' num2str(errmax) ')'])
        pause(0),
        h=h*max(0.95*(EPSTOL/errmax)^(1/4),1/10);
        if h<MINTIMESTEP, error(['Stepsize underflow; h = ' num2str(h)]), end
        rejected_steps=rejected_steps+1;
    else, % completed step: the time step is increased
        t=t+h;
        disp(['Now at t=' num2str(t) ', using h = ' num2str(h) ' (errmax=' num2str(errmax) ')'])
        pause(0),
        for i=[1 3 4 6], uvhat=uvhat+rk_c(i)*munutilde(:,:,i); end,
        uvhat=E.*uvhat;
        u = real(ifft2(uvhat(:,1:N)));
        v = real(ifft2(uvhat(:,N+1:end)));
        if t >= min(tkeep(end)+hkeep,tfinal),
            ukeep = cat(3,ukeep,u);
            vkeep = cat(3,vkeep,v);
            tkeep = [tkeep, t];
        end,
        h=min([h*min(0.95*(EPSTOL/errmax)^(1/5),5), hmax, tkeep(end)+hkeep-t, tfinal-t]);
        accepted_steps=accepted_steps+1;
        tcheck=[tcheck; t];
    end,
end
disp(['rejected_steps=' int2str(rejected_steps), ...
        ', accepted_steps=' int2str(accepted_steps)]),
save('gray2D_CK45.mat','tkeep','ukeep','vkeep','tcheck','N','L','x','y','tfinal','hkeep')
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
