function fout=fourierupsample2D(fin,newNx,newNy);
% Given a periodic function fin in 2D computed at equidistant
% nodes Nx x Ny, then fout is its upsampled version on
% newNx x newNy nodes.

[Ny,Nx]=size(fin);
HiFx=(Nx-mod(Nx,2))/2+1; HiFy=(Ny-mod(Ny,2))/2+1;
fftfin=max((newNx/Nx)*(newNy/Ny),1)*fft2(fin);
fout=real(ifft2([fftfin(1:HiFy,1:HiFx), ...
	zeros(HiFy,newNx-Nx), fftfin(1:HiFy,HiFx+1:end); ...
	zeros(newNy-Ny,newNx); fftfin(HiFy+1:end,1:HiFx), ...
	zeros(Ny-HiFy,newNx-Nx), fftfin(HiFy+1:end,HiFx+1:end)]));
