function [cmin,cmax]=showdata3d(N,feld,mal)

% SHOWDATA3D show datum points of 3d data
% showdata3d(N[,field[,options]])
% N..Structure of electrode numbers(a,b,m,n), 
%    k-factors(k) and measurements(r)
%    elec-Electrode Positions
% field - plot field instead of N.r

if nargin<2, feld=N.r; end
if nargin<3, mal=struct('xdir',0,'ydir',0,'vie',0); end
if ~isfield(mal,'vie'), mal.vie=1; end
if isempty(feld), feld=N.r; end
if mal.vie,
    if ~isfield(N,'zweid'), N=getpseudos(N); end
    [cmin,cmax]=plotprofiles(N,mal.vie-1,feld);
    if nargout<2, cmin=[cmin cmax]; end
    return
end
konfs=unique(sort(abs(N.k)));
[aa,bb]=meshgrid(abs(N.k),konfs);
[ik,jk]=find((aa-bb)==0);
midpoint=(N.elec(N.a,1:2)+N.elec(N.m,1:2))/2; % 2-Punkt-Anordnungen
fn=find(N.n>0);
if ~isempty(fn),  % 3-Punkt-Anordnungen
    %midpoint(fn,:)=(2*N.elec(N.a(fn),1:2)+N.elec(N.m(fn),1:2)+...
    %    N.elec(N.n(fn),1:2))/4;
    midpoint(fn,:)=(N.elec(N.m(fn),1:2)+N.elec(N.n(fn),1:2))/2;
end
fb=find(N.b>0);
if ~isempty(fb), % 4-Punkt-Anordnungen
    midpoint(fb,:)=(N.elec(N.a(fb),1:2)+N.elec(N.b(fb),1:2)+...
        N.elec(N.m(fb),1:2)+N.elec(N.n(fb),1:2))/4;
end
midx=unique(sort(midpoint(:,1)));
midy=unique(sort(midpoint(:,2)));
dmx=min(diff(midx));
midx=min(midx):dmx:max(midx);
dmy=min(diff(midy));
midy=min(midy):dmy:max(midy);
[aa,bb]=meshgrid(midpoint(:,1),midx);
[ix,jx]=find((aa-bb)==0);
[aa,bb]=meshgrid(midpoint(:,2),midy);
[iy,jy]=find((aa-bb)==0);
datums=ones(length(midx),length(midy),length(konfs))*NaN;
for l = 1:length(feld),
    datums(ix(l),iy(l),ik(l))=feld(l);
end
mal.elec=[];
mal.cont=[];
mal.cauto=1;
mal.alpha=0;
%mal.xy=0;
%mal.vie=2;
if min(feld(:))<0,
    mal.log=0;
    mal.cmap=2;
    mal.cmax=max(abs(feld(:)));
    mal.cauto=0;
    mal.cmin=-mal.cmax;
end
xx=zeros(length(midx)+1,1);
yy=zeros(length(midy)+1,1);
xx(2:end-1)=(midx(1:end-1)+midx(2:end))/2;
yy(2:end-1)=(midy(1:end-1)+midy(2:end))/2;
xx(1)=2*midx(1)-xx(2);
xx(end)=2*midx(end)-xx(end-1);
yy(1)=2*midy(1)-yy(2);
yy(end)=2*midy(end)-yy(end-1);
%draw3dgridmodel(datums,xx,yy,[0;1;konfs(:)],mal);
[cmin,cmax]=draw3dgridmodel(datums,xx,yy,[0 1 2:(length(konfs)+1)],mal);
if nargout<2, cmin=[cmin cmax]; end