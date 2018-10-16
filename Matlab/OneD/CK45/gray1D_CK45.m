function gray1D_CK45(N,tfinal,h,hkeep,L,epsilon,a,b)
if nargin<8;
    disp('Using default parameters');
    N=512;
    tfinal=2000;
    h=0.2;
    hkeep=2;
    L=50;
    epsilon=0.01;
    a=9*epsilon;
    b=0.4*epsilon^(1/3);
end

rhside_par={a,b};
x=(2*L/N)*(-N/2:N/2-1)';
uv=initial(x,L);
uvhat=fft(uv);
uvkeep(:,:,1)=uv;
tkeep=0; tcheck=0;

% Cash-Karp coefficients
rk_a=[0 1/5 3/10 3/5 1 7/8];
rk_b=[0 0 0 0 0; 1/5 0 0 0 0; 3/40 9/40 0 0 0; 3/10 -9/10 6/5 0 0; ...
        -11/54 5/2 -70/27 35/27 0; 1631/55296 175/512 575/13824 44275/110592 253/4096];
rk_c=[37/378 0 250/621 125/594 0 512/1771];
rk_d=[2825/27648 0 18575/48384 13525/55296 277/14336 1/4];

ksq=((pi/L)*[0:N/2 -N/2+1:-1]').^2;
hmax=0.9; h_old=0; t=0;
MINTIMESTEP=1e-10; EPSTOL=1e-4;
rejected_steps=0; accepted_steps=0;
while h >= MINTIMESTEP,
    if abs(h-h_old)>MINTIMESTEP,
        E_ai(:,:,1)=ones(N,2);
        for i=2:6, E_ai(:,:,i)=[exp(-ksq*rk_a(i)*h), exp(-epsilon*ksq*rk_a(i)*h)]; end,
        E=[exp(-ksq*h), exp(-epsilon*ksq*h)];
    end,

    % Embedded RK45 step
    fuv=rhside(uv,rhside_par{:});
    munutilde(:,:,1)=h*fft(fuv);
    deltauvhat=(rk_c(1)-rk_d(1))*munutilde(:,:,1);
    for i=2:6,
        sumbmunutilde=rk_b(i,1)*munutilde(:,:,1);
        for j=2:i-1, sumbmunutilde=sumbmunutilde+rk_b(i,j)*munutilde(:,:,j); end,
        uvhat_ai=E_ai(:,:,i).*(uvhat+sumbmunutilde);
        fuv=rhside(real(ifft(uvhat_ai)),rhside_par{:});
        munutilde(:,:,i)=h*fft(fuv)./E_ai(:,:,i);
        deltauvhat=deltauvhat+(rk_c(i)-rk_d(i))*munutilde(:,:,i);
    end,
    deltauvhat=E.*deltauvhat;
    deltauv=real(ifft(deltauvhat));
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
        uv = real(ifft(uvhat));
        if t >= min(tkeep(end)+hkeep,tfinal),
            uvkeep = cat(3,uvkeep,uv);
            tkeep = [tkeep, t];
        end,
        h=min([h*min(0.95*(EPSTOL/errmax)^(1/5),5), hmax, tkeep(end)+hkeep-t, tfinal-t]);
        accepted_steps=accepted_steps+1;
        tcheck=[tcheck; t];
    end,
end
disp(['rejected_steps=' int2str(rejected_steps), ...
        ', accepted_steps=' int2str(accepted_steps)]),
save('gray1D_CK45.mat','tkeep','uvkeep','tcheck','N','L','x','tfinal','hkeep')
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
