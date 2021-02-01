function [out,dM]=onedinv(N,Mod,R)

% ONEDINV - 1D inversion of dataset
% Bg = onedinv(N,Mod[,rho_0])
%
% [newModel,dM] = onedinv(N,Model)

%minKonf=min(abs(N.k)); 
%rq=median(N.r(find(abs(N.k)==minKonf)));

if nargin<3, R=Mod.Bg(1); end
if nargin<2, error('Two input arguments required'); end
rq=median(N.r);
if isstruct(Mod),
    z=Mod.z;
    if nargin>2,
        Bg=Mod.Bg;
        rq=R;
    else
        Bg=ones(length(z),1)*rq;
    end
else
    z=Mod;
    Bg=ones(length(z),1)*rq;
end
dR=log(N.r)-log(rq);
S=eindsens(z,N);
global INV
if isfield(INV,'auto'), inv=INV.auto; end
INV.auto=1;
D=100;
% C=1; % no weighting
nr=length(N.r);
if isfield(N,'err'), D=spdiags(1./log(1+N.err(:)),nr,nr); end
%dBg=invshifted(S,dR,1,0.01,0.8);
%c=sum(abs(S));c=c/max(c);C=spdiags(1./c(:),0,8,8);
one=ones(size(S,2),1);
C=spdiags([-one one],[0 1],size(S,2)-1,size(S,2));
L=C'*C;
% dBg=invshiftcdp(S,dR,L,D,1,200,1,0.8);
dBg=cglscdp(S,dR,30,L,D);
INV.auto=inv;
Bg(:)=Bg(:).*exp(dBg);
di=Bg(2:end)./Bg(1:end-1);
fi=find(di<1);
di(fi)=1./di(fi);
fac=1.2;
fi=find(di<fac);
while ~isempty(fi),
    mi=min(fi);
    ma=mi+1;
    while find(fi==ma), ma=ma+1; end
    Bg(mi:ma)=exp(mean(log(Bg(mi:ma))));
    di(mi:ma)=2;
    fi=find(di<fac);
end
Bg=rndig(Bg);
if isstruct(Mod),
    out=Mod;
    out.Bg=Bg(:)';
    if iscell(Mod.M), % para model
        dM=[];
        for k=1:length(out.M),
            out.M{k}(:)=Bg(k);
            nc=numel(out.M{k});
            dM=[dM;ones(nc,1)*dBg(k)];
        end
    else % grid model
        dM=Mod.M;
        for k=1:size(Mod.M,3), 
            out.M(:,:,k)=Bg(k); 
            dM(:,:,k)=dBg(k);
        end
    end
else
    out=Bg;
end

function Soll=eindsens(Z,N)
% computes 1d-sensitivity
amq=[];bmq=amq;anq=amq;bnq=amq;
data=length(N.r);
for l = 1:data,
    aaa=N.elec(N.a(l),:);
    mmm=N.elec(N.m(l),:);
    if N.b(l)>0, bbb=N.elec(N.b(l),:); else bbb=[Inf Inf 0]; end
    if N.n(l)>0, nnn=N.elec(N.n(l),:); else nnn=[-Inf -Inf 0]; end
    amq=[amq;sum((aaa-mmm).^2)];    
    anq=[anq;sum((aaa-nnn).^2)];    
    bmq=[bmq;sum((bbb-mmm).^2)];    
    bnq=[bnq;sum((bbb-nnn).^2)];    
end
%amq=0.25*amq;bmq=0.25*bmq;anq=0.25*anq;bnq=0.25*bnq;
bnq(isnan(bnq))=Inf;

Soll=[];
scal=1./sqrt(amq)-1./sqrt(anq)-1./sqrt(bmq)+1./sqrt(bnq);
ss=scal;
for kk = 2:length(Z),
  oss=ss;
%   zq=repmat(4*Z(kk)^2,data,1);
%   ss=1./sqrt(amq+zzq)-1./sqrt(anq+zzq)-1./sqrt(bmq+zzq)+1./sqrt(bnq+zzq);
  zq=Z(kk)^2;
  ss=1./sqrt(amq+zq)-1./sqrt(anq+zq)-1./sqrt(bmq+zq)+1./sqrt(bnq+zq);
  soll=(oss-ss)./scal; % obere Grenze-untere Grenze
  Soll=[Soll soll];
end
Soll=[Soll 1-sum(Soll,2)]; % Rest für z(end) bis inf
