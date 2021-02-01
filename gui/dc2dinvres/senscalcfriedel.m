function s = senscalcfriedel(x1,x2,z1,z2,nw,a,m,n,b)

%SENSCALCFRIEDEL
% 2D sensitivity of a cell within x1,x2,z1,z2
% Usage: S=senscalcfriedel(x1,x2,z1,z2,nw,a,m[,n[,b]])')
% nw number of weights
% a,m,n,b: vectors of electrode positions

if nargin<7
    error('Usage: S=senscalc(x1,x2,z1,z2,nw,a,m[,n[,b]])');
end
if length(a)~=length(m)
    error('length(a) ~= length(m)');
end
[x,wx]=gauleg(x1,x2,nw);
[z,wz]=gauleg(z1,z2,nw);

[Z,X]=meshgrid(z,x);
[WZ,WX]=ndgrid(wz,wx);

s=zeros(size(a,1),1);

for l = 1:size(a,1),
    S=integral(a(l,:),m(l,:),X,Z);
    if (nargin>7)&&(n(l,1)~=999), % A M N
        S=S-integral(a(l,:),n(l,:),X,Z);
        if (nargin>8)&&(b(l,1)~=999),
            S=S-integral(m(l,:),b(l,:),X,Z)+integral(n(l,:),b(l,:),X,Z);
        end
    end
    s(l)=sum(sum(WX.*WZ.*S));
end

if nargout==0,
    ss=S';
    ss(end+1,end+1)=0;
    xx=[x1 (x(1:end-1)+x(2:end))/2 x2];
    zz=[z1 (z(1:end-1)+z(2:end))/2 z2];
    mm=max(abs(S(:)));
    mm=mm/2;
    surf(xx,zz,ss);
    view(0,-90)
    shading faceted
    colormap('b2r');
    caxis([-mm mm]);
    colorbar
end

function I=integral(a,m,X,Z)
neg=[1 0;0 -1];
I=untegral(X,Z,a,m);
if a(2)>0,
    I=I+untegral(X,Z,neg*a(:),m(:));
    if m(2)>0,
        I=I+untegral(X,Z,neg*a(:),neg*m(:));
    end
end
if m(2)>0,
    I=I+untegral(X,Z,a,neg*m(:));
end    
I=I/(4*pi^2);
if a(2)>0, I=I/2; end
if m(2)>0, I=I/2; end

function I=untegral(X,Z,a,m)
A=(X-a(1)).^2+(Z-a(2)).^2;
B=(X-m(1)).^2+(Z-m(2)).^2;
n=find(B>A);
DD=A(n);
A(n)=B(n);
B(n)=DD;
T1=zeros(size(A));T2=T1;
n=find(A-B>1e-12);
[K,E]=ellipke(1-B(n)./A(n));
T2(n)=2./(sqrt(A(n)).*(A(n)-B(n)).^2);
T1(n)=T2(n)./B(n).*((A(n)+B(n)).*E(:)-2*B(n).*K(:));
T2(n)=T2(n).*((A(n)+B(n)).*K(:)-2*A(n).*E(:));

n=find(abs(A-B)<=1e-10);
B(n)=sqrt(A(n));
T1(n)=3*pi./(8*A(n).*A(n).*B(n));
T2(n)=pi./(8.*A(n).*B(n));
I=T1.*((X-a(1)).*(X-m(1))+(Z-a(2)).*(Z-m(2)))+T2;
