function [x,lam]=diskrep(A,b,VD,s,VM)

% DISKREP - discrepancy principle
% [x,lam] = diskrep(A,b[,VD,s,VM])

if nargin<5, [VD,s,VM]=csvd(A); end
u=0; %logs
o=2.5;
d=length(b);
xu=gtik(VD,s,VM,b,exp(u));
xo=gtik(VD,s,VM,b,exp(o));
tu=anorm(A,b,xu,d);
to=anorm(A,b,xo,d);
t=to;
while 1,
    if (tu-1)*(to-1)>0, % beide drüber o. drunter
        if to<1, % obere schon drunter
            u=o;xu=xo;tu=to;
            o=o+1;
            xo=gtik(VD,s,VM,b,exp(o));
            to=anorm(A,b,xo,d);
        else % untere drunter
            o=u;xo=xu;to=tu;
            u=u-1;
            xu=gtik(VD,s,VM,b,exp(u));
            tu=anorm(A,b,xu,d);
        end
    else % drüber/drunter
       %n=u-(o-u)/(to-tu)*(tu-1);
       n=(u+o)/2;
       x=gtik(VD,s,VM,b,exp(n));
       t=anorm(A,b,x,d);
       if t<1, tu=t;xu=x;u=n; else
           to=t;xo=x;o=n; end
    end
    if abs(t-1)<0.01, break; end
end
lam=exp(n)^2;
    
function tt=anorm(A,b,x,d)
tt=sum((b-A*x).^2)/d;
%tt=sum(abs(b-A*x))/d;