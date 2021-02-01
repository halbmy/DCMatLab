function [Rm,RM,RD,s,lambda]=getrm(S,N,dR)

% GETRM - Get Model Resolution from SVD
% [Rm,RM,RD,s] = getrm(S,N,lambda)  OR
% [Rm,RM,RD,s] = getrm(S,N,dR)

if nargin<3, dR=1; end
data=length(N.r);
if isfield(N,'err'),
    err=N.err;
else
    Umin=100e-6;I=100e-3;proz=1;
    err=abs(N.k./N.r)*Umin/I+proz/100;
    fprintf('Error min=%.1f%% max=%.1f%% mean=%.1f%%\n',...
        min(err)*100,max(err)*100,round(mean(err)*100));
end
D=spdiags(1./log(1+err),0,data,data);
[U,W,V]=svd(D*S);
s=diag(W);
clear W
if length(dR)==1,
    lambda=sqrt(dR); 
else
    lambda=gcv(U,s,dR);
    fprintf('GCV estimation for lambda=%g\n',lambda);
end
r=max(find(s>eps*s(1)));
s=s(1:r);V(:,r+1:end)=[];
f=s.^2./(s.^2+lambda^2);
m=size(V,1);
f=f(:)';
Rm=zeros(m,1);
for i=1:m,
    Rm(i)=sum(V(i,:).^2.*f);
end
if nargout>2,
    RM=V(:,1:r)*diag(f)*V(:,1:r)';
end
if nargout>3,
    RD=U(:,1:r)*diag(f)*U(:,1:r)';
end
fprintf('Information content = %.1f (%d%% model) (%d%% data)\n',...
    sum(Rm),round(mean(Rm)*100),round(sum(Rm)/data*100));
