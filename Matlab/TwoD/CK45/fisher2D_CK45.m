function fisher2D_CK45(N,tfinal,h,hkeep,L,delta)
if nargin<6
    disp('Using default parameters')
    tfinal=20;
    h=0.1;
    hkeep=1;
    N=128;
    L=20;
    delta=0.25;
end

rhside_par={delta};
x=(2*L/N)*(-N/2:N/2-1)'; y=x;
u=initial(x,y,rhside_par{:});
uhat=fft2(u);
ukeep(:,:,1)=u;
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
        E_ai(:,:,1)=ones(N,N);
        for i=2:6, E_ai(:,:,i)=exp(-ksq*rk_a(i)*h); end,
        E=exp(-ksq*h);
    end,
    
    % Embedded RK45 step
    fu=rhside(u,rhside_par{:}); 
    mutilde(:,:,1)=h*fft2(fu);
    deltauhat=(rk_c(1)-rk_d(1))*mutilde(:,:,1);
    for i=2:6,
        sumbmutilde=rk_b(i,1)*mutilde(:,:,1);
        for j=2:i-1, sumbmutilde=sumbmutilde+rk_b(i,j)*mutilde(:,:,j); end,
        uhat_ai=E_ai(:,:,i).*(uhat+sumbmutilde);
        fu=rhside(real(ifft2(uhat_ai)),rhside_par{:}); 
        mutilde(:,:,i)=h*fft2(fu)./E_ai(:,:,i);
        deltauhat=deltauhat+(rk_c(i)-rk_d(i))*mutilde(:,:,i);
    end,
    deltauhat=E.*deltauhat;
    deltau=real(ifft2(deltauhat));
    errmax = max(max(abs(deltau)));

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
        for i=[1 3 4 6], uhat=uhat+rk_c(i)*mutilde(:,:,i); end,
        uhat=E.*uhat;
        u = real(ifft2(uhat));
        if t >= min(tkeep(end)+hkeep,tfinal),
            ukeep = cat(3,ukeep,u); 
            tkeep = [tkeep, t];
        end,
        h=min([h*min(0.95*(EPSTOL/errmax)^(1/5),5), hmax, tkeep(end)+hkeep-t, tfinal-t]);
        accepted_steps=accepted_steps+1;
        tcheck=[tcheck; t];
    end,
end
disp(['rejected_steps=' int2str(rejected_steps), ...
        ', accepted_steps=' int2str(accepted_steps)]),
save('fisher2D_CK45.mat','tkeep','ukeep','tcheck','N','L','x','y','tfinal','hkeep')
%--------------Initial Condition ---------------------------------
function u=initial(x,y,delta)
[xx,yy]=meshgrid(x,y);
u=exp(-delta*(xx.^2+yy.^2))/2;
%---------------Right Hand Side-----------------------------------
function rhsu=rhside(u,delta)
rhsu=u.*(1-u);
