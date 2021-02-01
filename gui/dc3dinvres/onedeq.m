function A=onedeq(N,Bg,z,proz)

% ONEDEQ - 1D equivalent models
% onedeq(Rho,z,proz)

if nargin<4, proz=6; end
R=fwd3d1d(N,Bg,z);
for k=1:length(Bg),
    maBg=Bg;miBg=Bg;
    while (rms(R,fwd3d1d(N,maBg,z))<proz), maBg(k)=maBg(k)*(1+proz/100); end
    fprintf('.');
    while (rms(R,fwd3d1d(N,miBg,z))<proz), 
        miBg(k)=miBg(k)/(1+proz/100); 
        if miBg(k)<1, break; end
    end
    fprintf('.');
    mi(k)=miBg(k);ma(k)=maBg(k);
end
A=round([mi(:) Bg(:) ma(:)]);
if nargout<1,
    zm=z+0.5;
    loglog(Bg,zm,'r-x',mi,zm,'b:',ma,zm,'b:');axis ij
    loglog(Bg,zm,'r-x');axis ij
    hold on
    for k=1:length(Bg), line([mi(k) ma(k)],[zm(k) zm(k)]); end
    hold off
    xlabel('resistivity [Ohm.m]');ylabel('depth [m]');grid on;
end