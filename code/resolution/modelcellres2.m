function x=modelcellres(S,lam,n,L,D,P)

% MODELCELLRES - Model cell resolution by approximation
% res = modelcellres(S,lam,n,L,D,P)
% S..Sensitivity matrix
% lam..regularization parameter
% n..Cell number
% L..(quadratic) constraint matrix
% D..data weighting matrix
% P..parameter mapping matrix

if nargin<6, P=1; end
if nargin<5, D=1; end
if nargin<4, L=1; end
if nargin<3, error('Three input arguments required!'); end

% x=cglscdp(S,S(:,n),lam,L,D);return
scal=1;
b=S(:,n);
x=zeros(size(S,2),1);x(n)=1;
% if size(P,1)==size(S,2), 
%     x=zeros(size(P,2),1);x(find(P(n,:)))=1;; end
PI = P\speye(size(P,1));
x=PI*x;
z = D*(b - S*(P*x)); % residuum of unregularized equation
p = (z'*D*S*P)'-PI*(L*(P*x))*lam; % residuum of normal equation
r = p;normr2 = r'*r;
acc=1e-4;
abbr = normr2*acc; % goal for norm(r)^2
j=0;
while(normr2>abbr)
    j=j+1;  
    Pp=P*p;
    q = D*(S*Pp);
    normr2old=normr2;
    alpha = normr2/(q'*q+Pp'*(L*Pp)*lam);
    x  = x + alpha*p;
    z  = z - alpha*q;
    r = (z'*D*S*P)'-PI*(L*(P*x))*lam;
    normr2 = r'*r;
    beta = normr2/normr2old;
    p = r + beta*p;
end
