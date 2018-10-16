% ADI for Fisher's equation in 2D 

tic
L=50;% domain is [-L L], periodic
N=256;
h=2*pi/N; x=h*(1:N); x=L*(x-pi)/pi;
y=x;
[xx,yy]=meshgrid(x,y);
%Fourier Differentiation matrix (see Trefethen's Spectral Method book)
column=[-pi^2/(3*h^2)-1/6 ...
    -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2];
D2=(pi/L)^2*toeplitz(column); 
time=0;
% time step
dt=0.2;
matrix1=eye(N)+dt/2*D2;
matrix2=inv(eye(N)-dt/2*D2);
matrix3=eye(N)+dt/2*D2';
matrix4=inv(eye(N)-dt/2*D2');
%initial condition
u0 = [1*exp(-0.25*(xx.^2+yy.^2))];
count=0;
ukeep=u0;

for i=1:50
count=count+1;
time=time+dt;
rhs0=dt/2*u0.*(1-u0);
rhs=0;
u0=matrix2*u0*matrix3+matrix2*rhs0;

rhs0=dt/2*u0.*(1-u0);
rhs0=0;
u0=matrix1*u0*matrix4+rhs0*matrix4;

if count==5 
time
ukeep=[ukeep u0];
count=0;
end

end

tau=toc

figure(1)
temp=reshape(ukeep,256,256,11);
plot(x,squeeze(temp(N/2,:,:)))
ukeep2=ukeep;
save('ADI2.mat','ukeep2','x','y','temp','dt','tau');
