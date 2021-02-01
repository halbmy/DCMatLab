function RM = approxrm(S,lam,C,D)

% APPROXRM - Approximate resolution matrix
% RM=approxrm(S,lam,C,D)

if nargin<4, D=1; end
if nargin<3, C=1; end
if nargin<2, error('Two input arguments required!'); end
P=1;PI=1;L=C'*C;
vk=1;
if vk,
    scal=zeros(size(S,2),1);
    sdi=zeros(size(S,1),1);
    for i=1:size(S,2),
        sdi=D*S(:,i);
        scal(i)=sdi'*sdi;
    end
    scal=sqrt(scal+lam*diag(L));
    for i=1:length(scal),
        S(:,i)=S(:,i)/scal(i);
        C(:,i)=C(:,i)/scal(i);
    end
    L=C'*C;
else
    scal=1;
end
nd=size(S,1);nm=size(S,2);
acc=1e-4;
RM=zeros(nm);
x=zeros(nm,1);
for n=1:size(RM,2),
    b=S(:,n);
    x(n)=1;
    z = D*(b - S*(P*x)); % residuum of unregularized equation
    p = (z'*D*S*P)'-PI*(L*x)*lam; % residuum of normal equation
    r = p;normr2 = r'*r;
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
    fprintf('%d ',j);
    RM(:,n)=x;
    x(n)=0;
end
