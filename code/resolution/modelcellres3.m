function x=modelcellres(S,lam,L,D,n)

% MODELCELLRES - Model cell resolution by approximation
% mcr = modelcellres(S,lambda,L,D,n)

if nargin<4, D=1; end
if nargin<3, L=1; end
if nargin<2, error('Two input arguments required!'); end

% x=cglscdp(S,S(:,n),lam,L,D);return
scal=1;
b=S(:,n);
x=zeros(size(S,2),1);
x(n)=1;
P=1;PI=1;
z = D*(b - S*(P*x)); % residuum of unregularized equation
p = (z'*D*S*P)'-PI*(L*x)*lam; % residuum of normal equation
r = p;normr2 = r'*r;
acc=1e-4;
abbr = normr2*acc; % goal for norm(r)^2
j=0;
while(normr2>abbr)
    j=j+1;  
    q = D*(S*(P*p));
    normr2old=normr2;
    Pp=P*p;
    alpha = normr2/(q'*q+Pp'*(L*Pp)*lam);
    x  = x + alpha*p;
    z  = z - alpha*q;
    r = (z'*D*S*P)'-PI*(L*(P*x))*lam;
    normr2 = r'*r;
    beta = normr2/normr2old;
    p = r + beta*p;
end
