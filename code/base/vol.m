function vol(M,Xm,Ym,Zm,val,sx,sy,sz,caps)
%% Isosurface/slice visualization
%% vol(M,xm,ym,zm,iso_value,xslices,yslices,zslices,caps)
global MAL
sq=100;
Im=length(Xm);Jm=length(Ym);Km=length(Zm);
figure(4);
clf
%% Variablen ...
if nargin<5, val=120; end
if nargin<6, sx=[]; end
if nargin<7, sy=[]; end
if nargin<8, sz=[]; end
if nargin<9, caps=0; end
%% Mesh-Grid
%zz=[0 Zm(find(Zm<100))];
zz=Zm;
zm=length(zz);
if MAL.xy==0,
  W=zeros(Im,Jm,Km+1);
  W(:,:,2:Km+1)=reshape(M,Im,Jm,Km);
  [x,y,z]=meshgrid(Ym,Xm,zz);
  xlab='x in m';
  ylab='y in m';
else
  W=zeros(Jm,Im,Km);
  W(:,:,2:Km+1)=reshape(rot90(M),Jm,Im,Km);
  [x,y,z]=meshgrid(Xm,Ym,zz);
  xlab='y in m';
  ylab='x in m';
end
W(:,:,1)=1/sq;
%W=smooth3(W,'box',5);
if MAL.log==1,
  W=log10(W);
  val=log10(val);
end
%% Slices
sl=slice(x,y,z,W,sx,sy,sz);
set(sl,'EdgeColor','none','FaceColor','interp');
colormap default
switch MAL.cmap
  case 1, colormap default
  case 2, colormap(b2r);
  case 3, colormap hot
  case 4, colormap gray
  case 5, colormap jet
  case 6, colormap cool
  otherwise, colormap default
end
if MAL.cauto==1,
    %if MAL.log==1,
    %  [NN,VV]=hist(log10(W(:)),100); 
    %else
      [NN,VV]=hist(W(:),100);
    %end
    CN=cumsum(NN);CN=CN/max(CN);
    imin=max(find(CN<0.01));
    imax=min(find(CN>0.99));
    if isempty(imin), imin=1; end
    if isempty(imax), imax=length(VV); end
    cmin=VV(imin);
    cmax=VV(imax);
    caxis([cmin cmax]);
else
    if MAL.log==1, 
        caxis(log10([MAL.cmin MAL.cmax])); 
    else 
        caxis([MAL.cmin MAL.cmax]); 
    end
end

%% Isosurface
cmap=colormap;
cach=caxis;
for l = 1:length(val),
	cind=1+round((val(l)-cach(1))/(cach(2)-cach(1))*(length(cmap)-1));
	if cind>length(cmap), cind=length(cmap); end
	if cind<1, cind=1; end
	col=cmap(cind,:);
    iso=patch(isosurface(x,y,z,W,val(l)),'FaceColor',col,'EdgeColor','none',...
        'SpecularExponent',5,'SpecularColorReflectance',0.05);
      %isonormals(W,iso)
    if caps>0,
      patch(isocaps(x,y,z,W,val(l),caps),'FaceColor','interp','EdgeColor','none')
    end
end
%% Einstellungen
set(gca,'ZDir','reverse');
%zz=get(gca,'zlim');
%zz(1)=0;
%set(gca,'zlim',zz);
if MAL.xy==0,
  set(gca,'ylim',[min(Xm) max(Xm)],'xlim',[min(Ym) max(Ym)]);
else
  set(gca,'xlim',[min(Xm) max(Xm)],'ylim',[min(Ym) max(Ym)]);
end
set(gca,'zlim',[min(zz) max(zz)]);
xdir='normal';
ydir='normal';
if MAL.xy==0,
  if MAL.xdir==1, xdir='reverse';  end
  if MAL.ydir==1, ydir='reverse';  end
else
  if MAL.xdir==1, ydir='reverse';  end
  if MAL.ydir==1, xdir='reverse';  end
end
set(gca,'XDir',xdir);
set(gca,'XDir',ydir);
xlabel(xlab);
ylabel(ylab);
zlabel('z in m');
set(gca,'DataAspectRatio',[1 1 1]);
drawnow;
axis vis3d
zoom(1.5)
camlight headlight
set(gcf,'Renderer','zbuffer');
lighting gouraud
azn=200;
el=25;
view(azn,el);
drawnow
