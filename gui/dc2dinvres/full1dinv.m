function [Bg,R]=full1dinv(N,z,lam,niter)

% FULL1DINV - Full 1D inversion of 3D data
% Bg=full1dinv(N,z)

if nargin<4, niter=10; end
if nargin<3, lam=10; end
if nargin<2, global Model;z=Model.z; end
if nargin<1, global N; end
S=sens1d(z,N);
STS=S'*S;
rho=median(N.r);
Bg=ones(length(z),1)*rho;
R=ones(length(N.a),1)*rho;
i=0;
fprintf('It %d: min/max=%g/%g RMS=%.2f%%\n',i,rndig(min(Bg)),rndig(max(Bg)),rms(N.r,R,1));
for i=1:niter,
    dR=log(N.r)-log(R);
    dBg=(STS+lam*eye(length(z)))\(S'*dR);
    nBg=Bg.*exp(dBg);
    oldR=R;
    R=fwd2d1d(N,nBg,z);
    tau=linesearch(N,oldR,R);
    Bg=Bg.*exp(tau*dBg);
    R=fwd2d1d(N,Bg,z);
    fprintf('It %d: tau=%.2f min/max=%g/%g RMS=%.2f%%\n',...
        i,tau,rndig(min(Bg)),rndig(max(Bg)),rms(N.r,R,1));
    ddr=(log(R)-log(oldR)-S*dBg(:))/(dBg(:)'*dBg(:));ddr=ddr(:);
    for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dBg(i); end
end
% fprintf('%d ',round(Bg)');fprintf('\n');
% onedeq(N,Bg,z,3);