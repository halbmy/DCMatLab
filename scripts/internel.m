function newN=internel(N,x1,x2,y1,y2)

% INTERNEL - Limitation of electrode positions
% newN = internel(N,x1,x2,y1,y2)
%        newN,N..Datasets
%        newN lies completely within (x1,y1)-(x2,y2)

if nargin<5, y2=[]; end
if nargin<4, y1=[]; end
if nargin<3, x2=[]; end
if isempty(x1), x1=min(N.elec(:,1)); end
if isempty(x2), x2=max(N.elec(:,1)); end
if isempty(y1), y1=min(N.elec(:,2)); end
if isempty(y2), y2=max(N.elec(:,2)); end
newN=[];
fb=find(N.b);fn=find(N.n);
drin=(N.elec(:,1)>=x1).*(N.elec(:,1)<=x2)...
    .*(N.elec(:,2)>=y1).*(N.elec(:,2)<=y2);
fd=find(drin);ld=length(fd);
drin(fd)=(1:ld);
ld=ld+1;
drin(find(drin==0))=ld;
newN.elec=N.elec(fd,:);
newN.a=drin(N.a);
newN.m=drin(N.m);
newN.b=N.b(:);newN.n=N.n(:);
newN.b(fb)=drin(N.b(fb));
newN.n(fn)=drin(N.n(fn));
all=(newN.a<ld).*(newN.b<ld).*(newN.m<ld).*(newN.n<ld);
fa=find(all);
newN.a=newN.a(fa);
newN.b=newN.b(fa);
newN.m=newN.m(fa);
newN.n=newN.n(fa);
newN.r=N.r(fa);
newN.k=N.k(fa);
 
% x1a=ones(size(N.a));
% x2a=x1a;y1a=x1a;y2a=x1a;
% x1b=x1a;x2b=x1a;y1b=x1a;y2b=x1a;
% x1m=x1a;x2m=x1a;y1m=x1a;y2m=x1a;
% x1n=x1a;x2n=x1a;y1n=x1a;y2n=x1a;
% x1a(:)=(N.elec(N.a,1)>x1);
% x1b(fb)=(N.elec(N.b(fb),1)>=x1);
% x1m(:)=(N.elec(N.m,1)>=x1);
% x1n(fn)=(N.elec(N.n(fn),1)>=x1);
% x2a(:)=(N.elec(N.a,2)<=x2);
% x2b(fb)=(N.elec(N.b(fb),2)<=x2);
% x2m(:)=(N.elec(N.m,2)<=x2);
% x2n(fn)=(N.elec(N.n(fn),2)<=x2);
% y1a(:)=(N.elec(N.a,1)>=y1);
% y1b(fb)=(N.elec(N.b(fb),1)>=y1);
% y1m(:)=(N.elec(N.m,1)>=y1);
% y1n(fn)=(N.elec(N.n(fn),1)>=y1);
% y2a(:)=(N.elec(N.a,2)<=y2);
% y2b(fb)=(N.elec(N.b(fb),2)<=y2);
% y2m(:)=(N.elec(N.m,2)<=y2);
% y2n(fn)=(N.elec(N.n(fn),2)<=y2);
% all=x1a.*x2a.*y1a.*y2a.*x1b.*x2b.*y1b.*y2b.*...
%     x1n.*x2n.*y1n.*y2n.*x1n.*x2n.*y1n.*y2n;
% 
% [x1a;x2a;y1a;y2a;x1m;x2m;y1m;y2m;all]
% 
% new.a=N.a(all);
% new.b=N.b(all);
% new.m=N.m(all);
% new.n=N.n(all);
% new.r=N.r(all);
% new.k=N.k(all);
