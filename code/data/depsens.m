function [s,zmaxs,k]=depsens(A,B,M,N,z);

% DEPSENS - depth sensitivity function
% s         = depsens(A,B,M,N,z)
% [s,zmaxs] = ...

zq4=z.^2*4;
AMQ=sum((A-M).^2,2);
ANQ=sum((A-N).^2,2);
BMQ=sum((B-M).^2,2);
BNQ=sum((B-N).^2,2);
k=2*pi./(1./sqrt(AMQ)-1./sqrt(ANQ)-1./sqrt(BMQ)+1./sqrt(BNQ));
s=(1./sqrt(zq4+AMQ).^3-1./sqrt(zq4+ANQ).^3-1./sqrt(zq4+BMQ).^3+1./sqrt(zq4+BNQ).^3)*4*k.*z/2/pi;
zmaxs=mean(z(s==max(s)));
if nargout<1,
    plot(z,s);
end