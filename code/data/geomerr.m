function Gerr=geomerr(N,dx)

% GEOMERR - get geometrical error
% Gerr = geomerr(N,dx)

if nargin<2, dx=0.01; end
xa=N.elec(N.a,1);
xm=N.elec(N.m,1);
dk=abs(getk(xa+dx,xm)-N.k)+abs(getk(xa,xm+dx)-N.k);
fn=find(N.n);
if ~isempty(fn),
    xn=N.elec(N.n(fn),1);
    dk(fn)=abs(getk(xa(fn)+dx,xm(fn),xn)-N.k)+...
        abs(getk(xa(fn),xm(fn)+dx,xn)-N.k)+...
        abs(getk(xa(fn),xm(fn),xn+dx)-N.k);
end
fb=find(N.n.*N.b);
if ~isempty(fb),
    xn=N.elec(N.n(fb),1);
    xb=N.elec(N.b(fb),1);
    dk(fb)=abs(getk(xa(fb)+dx,xm(fb),xn,xb)-N.k)+...
        abs(getk(xa(fb),xm(fb)+dx,xn,xb)-N.k)+...
        abs(getk(xa(fb),xm(fb),xn+dx,xb)-N.k)+...
        abs(getk(xa(fb),xm(fb),xn,xb+dx)-N.k);
end
Gerr=dk./abs(N.k);
showdata(N,Gerr*100);

function k=getk(xa,xm,xn,xb)
kk=1./abs(xa-xm);
if nargin>2, kk=kk-1./abs(xa-xn); end
if nargin>3, kk=kk-1./abs(xb-xm)+1./abs(xb-xn); end
k=2*pi./kk;