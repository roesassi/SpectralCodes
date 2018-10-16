function fout=fourierupsample(fin,newN);
% Given a periodic function fin, computed at N equispaced nodes in
% the periodic domain [-L,L], fout is its upsampled version on newN
% nodes onto the same domain.

N=length(fin); HiF=(N-mod(N,2))/2+1;
fftfin=max((newN/N),1)*fft(fin);
fout=real(ifft([fftfin(1:HiF); zeros(newN-N,1); fftfin(HiF+1:N)]));
