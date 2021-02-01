function clog=tetvis(x,y,z,v,dd,val,sx,sy,sz)

% TETVIS - 3D Tetraeder visualization
% tetvis(x,y,z,v,dd,sx,sy,sz,val)
% x,y,z .. vectors of arbitrary points in 3d-space
% v     .. value vector accompanied to x,yz
% dd    .. grid size
% sx/sy/sz .. slice positions for x/y/z (can be vectors)
% val   .. isovalue for isosurface

if nargin==1,
    v=x;x=1:length(v);y=x;z=x;
elseif nargin<4,
    x=[0 0 0 0 1 1 1 1];
    y=[0 0 1 1 0 0 1 1];
    z=[0 1 0 1 0 1 0 1];
    v=sin(x).*cos(y)+sin(z.^2);
end
if nargin<5, dd=1; end
if length(dd)<2, d(2)=dd(1); end
if length(dd)<3, d(3)=dd(2); end
if nargin<6, val=mean(v); end
if nargin<7, sx=999; end
if nargin<8, sy=999; end
if nargin<9, sz=999; end
if min(v)>0, clog=1; else clog=0; end

if clog, v=log10(v);val=log10(val); end    
dx=dd(1);dy=dd(2);dz=dd(3);
xx=min(x):dx:max(x);yy=min(y):dy:max(y);zz=min(z):dz:max(z);
if sx==999, sx=max(xx); end
if sy==999, sy=max(yy); end
if sz==999, sz=max(zz); end
[X,Y,Z]=meshgrid(xx,yy,zz);
fprintf('Triangulating...');
tic;V=griddata3(x,y,z,v,X,Y,Z);toc
fprintf('ready\n');
sl=slice(X,Y,Z,V,sx,sy,sz);
% caxis(log10([90 250]))
for l=1:length(sl), set(sl(l),'EdgeColor','none'); end
cmap=colormap;
cach=caxis;
for l=1:length(val),
  	cind=1+round((val(l)-cach(1))/(cach(2)-cach(1))*(length(cmap)-1));
	if cind>length(cmap), cind=length(cmap); end
	if cind<1, cind=1; end
	col=cmap(cind,:);
    iso=patch(isosurface(X,Y,Z,V,val(l)),'FaceColor',col,'EdgeColor','black',...
        'SpecularExponent',5,'SpecularColorReflectance',0.05);
    isonormals(X,Y,Z,V,iso);
end
set(gca,'XLim',[min(xx) max(xx)]);
set(gca,'YLim',[min(yy) max(yy)]);
set(gca,'ZLim',[min(zz) max(zz)]);
set(gca,'ZDir','reverse');
set(gca,'DataAspectRatio',[1 1 1]);
xlabel('x in m')
ylabel('y in m')
zlabel('z in m')