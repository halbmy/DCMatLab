function [x,rho,eta] = gtik(U,sm,X,b,lambda)
% GTIK - generalized Tikhonov regularization
% x = gtik(U,sm,X,b,lambda)

x=zeros(size(X,1),length(lambda));
utb=(b'*U)';
    si=sm(:,1);
if size(sm)>1,
    my=sm(:,2);
else
    my=ones(size(si));
end
gs=si./my;
p=length(si);
for i=1:length(lambda)
    %fi=gs.^2./(gs.^2+lambda(i)^2);
%     x(:,i)=X(:,1:p)*(utb(1:p)./sm(:,1).*fi);
    x(:,i)=X(:,1:p)*(utb(1:p).*si./(si.^2+lambda(i)^2*my.^2));
    %+X(:,p+1:end)*utp(p+1:end);
    rho(i)=lambda(i)^2*norm(utb(1:p)./(gs + lambda(i)^2)); 
    eta(i)=norm(x(:,i));
end
