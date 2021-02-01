function Lay=get1dmodel(N,z,dR,rho,lam)
% GET1DMODEL - Starting model from data by 1D-Inversion
% Lay=get1dmodel(N,z,dR,rho)

if nargin<5, lam=0.1; end
if nargin<4, rho=median(N.r); end
if nargin<3, dR=log(N.r)-log(rho); end

amq=(N.elec(N.a,1)-N.elec(N.m,1)).^2;
anq=ones(size(amq))*Inf;
bmq=anq;
bnq=anq;
fn=find(N.n>0);
anq(fn)=(N.elec(N.a(fn),1)-N.elec(N.n(fn),1)).^2;
fb=find(N.b>0);
bmq(fb)=(N.elec(N.b(fb),1)-N.elec(N.m(fb),1)).^2;
bnq(fb)=(N.elec(N.b(fb),1)-N.elec(N.n(fb),1)).^2;

nz=length(z)-1;
%eindsens=zeros(length(N.a),nz);
eindsens=[];
scal=1./sqrt(amq)-1./sqrt(anq)-1./sqrt(bmq)+1./sqrt(bnq);
oben=scal;
for l = 1:nz,
    zzq=ones(size(amq))*(4*z(l+1)^2);
    unten=1./sqrt(amq+zzq)-1./sqrt(anq+zzq)-1./sqrt(bmq+zzq)+1./sqrt(bnq+zzq);
    %eindsens(:,1)=(oben-unten)./scal;
    eindsens=[eindsens rho*(oben-unten)./scal./N.r];
    oben=unten;
end
Lay=exp((eindsens'*eindsens+lam*eye(nz))\(eindsens'*dR))*rho;
%Lay=exp(invshifted(eindsens,dR,2,0.1,0.5))*rho;