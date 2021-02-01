function [mi,ma]=draw3dmodel(Mod,MAL,Field,Alph)

% DRAW3DMODEL - draw 3d model
% [mi,ma]=draw3dmodel(Model,MAL,Field,Alpha)
% Model - Model structure with x,y,z and M
% OPT - structure of possible fields
%   cauto - automatic coloring(1)
%   cmin/cmax - maximum/minimum color
%   cmap - colormapping(0)
%   clog - use logarithmic colortable
%   nu/nv - draw (nu x nv) subplots (3/2)
%   startwith - start with slice number n
%   xy - swap x/y (or x/z or y/z) Orientation
%   xdir/ydir - x/y Direction normal/reverse
%   cont - draw contour lines
%   vie - draw xy(0=default), xz(1) or yz(2) slices
% Field - draw Field instead of M
% Alph  - use Alph values for alpha shading

mal=struct('cauto',1,'cmin',100,'cmax',500,'cmap',0,'log',1,'xdir',0,'ydir',0,'clog',1,...
        'xy',0,'nu',0,'nv',2,'elec',1,'startwith',1,'cont',[],'vie',0,'alpha',0);  
%default struct
namedir={'normal','reverse'};
if nargin<4, Alph=[]; end
if nargin<2, MAL=mal; end
if ~iscell(Mod.M), % grid model
    MM=Mod.M;
    if (nargin>2)&&(prod(size(Field))==prod(size(MM))), 
        MM(:)=Field(:); 
    end
%     [mi,ma]=draw3dgridmodel(MM,Mod.x,Mod.y,Mod.z,MAL,Mod.Bg,Alph);
    [mi,ma]=patch3dgridmodel(MM,Mod.x,Mod.y,Mod.z,MAL,Mod.Bg,Alph);
    if nargout<2, mi=[mi ma]; end
    return;
end
nobreit=0;
if isempty(MAL), MAL=mal; end
if isfield(MAL,'xy'), xy=MAL.xy; else xy=1; end
if isfield(MAL,'xdir'), xdir=MAL.xdir; else xdir=0; end
if isfield(MAL,'ydir'), ydir=MAL.ydir; else ydir=0; end
if isfield(MAL,'cmax'), cmax=MAL.cmax; else cmax=1; end
if isfield(MAL,'cmin'), cmin=MAL.cmin; else cmin=0; end
if isfield(MAL,'cmap'), cmap=MAL.cmap; else cmap=0; end
if isfield(MAL,'cbar'), colbar=MAL.cbar; else colbar=0; end
if isfield(MAL,'log'), clog=MAL.log; else clog=0; end
if isfield(MAL,'clog'), clog=MAL.clog; end
if isfield(MAL,'cauto'), cauto=MAL.cauto; else cauto=1; end
if isfield(MAL,'nu'), nu=MAL.nu; else nu=0; end
if isfield(MAL,'nv'), nv=MAL.nv; else nv=0; end
if isfield(MAL,'elec'), elec=MAL.elec; else elec=0; end
if isfield(MAL,'startwith'), startwith=MAL.startwith; else startwith=1; end
if isfield(MAL,'cont'), cont=MAL.cont; else cont=[]; end
if isfield(MAL,'vie'), vie=MAL.vie; else vie=0; end
if isfield(MAL,'alpha'), alfa=MAL.alpha; else alfa=0; end
if nargin<4, Alph=[]; end

if (alfa)&(nargin<4), Alph=Field; else
    if (nargin>2)&&isequal(Mod.ncells,length(Field)), 
        Mod=modelupdate(Mod,Field,0); end
end
if alfa, Alp=modelupdate(Mod,Alph,0); end
if min(Mod.M(:))<0, clog=0; end
clf;
K=length(Mod.z)-1;
nx=size(Mod.M{1},1);ny=size(Mod.M{1},2);
xmax=Mod.x0+nx*Mod.dx*Mod.nx(1);
ymax=Mod.y0+ny*Mod.dy*Mod.ny(1);
nnx=Mod.nx(1)*nx;nny=Mod.ny(1)*ny;
if cauto==1,
    mi=1000;ma=0.001;
    for k=1:K,
        mi=min(mi,min(Mod.M{k}(:)));
        ma=max(ma,max(Mod.M{k}(:)));
    end
    if mi>=ma, mi=0.99*ma; end
    if mi*ma<=0, clog=0; end
else
    mi=cmin;ma=cmax;
end
if clog,
    mi=log10(mi);
    ma=log10(ma);
else
    if mi*ma<0,
        mm=max(abs([mi ma]));
        mi=-mm;ma=mm;
    end
end

if vie>0, % xz- or yz-slices
    [M,x,y]=mesch3dmodel(Mod);
    z=Mod.z;
    lm=size(M,2);
    if nu*nv==0,
        %rel=(size(Mod.M{1},1)+2)*Mod.dx/(size(Mod.M{1},3)+2)/Mod.dy;
        rel=(max(x)-min(x))/(max(z)-min(z));
        nv=fix(sqrt(lm*rel)+0.499);
        if nv>lm, nv=lm; end
        nu=round(lm/nv+0.499);
        while (nv-1)*nu>=lm, nv=nv-1; end
    end
    if clog, M=log10(M); end
    for l=1:lm,
        if lm>nu*nv, break; end
        if nu*nv>1, subplot(nv,nu,l); end
        mm=squeeze(M(:,l,:));
        mm(end+1,end+1)=0;
        pcolor(x,z,mm');
        axis equal tight
        set(gca,'YDir','reverse');
        caxis([mi ma]);
    end
    k=l;
else % xy-Slices
    if nu*nv==0,
        lm=length(Mod.M);
        rel=(size(Mod.M{1},1)+2)*Mod.dx/(size(Mod.M{1},2)+2)/Mod.dy;
        nv=fix(sqrt(lm*rel)+0.499);
        if nv>lm, nv=lm; end
        nu=round(lm/nv+0.499);
        while (nv-1)*nu>=lm, nv=nv-1; end
    end
    xmima=[Mod.x0 xmax];
    ymima=[Mod.y0 ymax];
    for k=1:min(K,nu*nv),
        if nu*nv>1, subplot(nv,nu,k); end
        modk=ones(size(Mod.M{k})+2*nobreit)*Mod.Bg(k)*clog;
        modk(1+nobreit:end-nobreit,1+nobreit:end-nobreit)=Mod.M{k};
        nx=size(Mod.M{k},1);ny=size(Mod.M{k},2);
        rx=floor(mod(nnx,nx*Mod.nx(k))/2);
        ry=floor(mod(nny,ny*Mod.ny(k))/2);
        x=((0:nx)*Mod.nx(k)+rx)*Mod.dx+Mod.x0;
        y=((0:ny)*Mod.ny(k)+ry)*Mod.dy+Mod.y0;
        if nobreit,
            x=[Mod.x0 x xmax];
            y=[Mod.y0 y ymax];
        else
            x(1)=Mod.x0;x(end)=xmax;
            y(1)=Mod.y0;y(end)=ymax;
        end
        modk(:,end+1)=1;modk(end+1,:)=1;
        if clog,
            pcolor(x,y,log10(modk)');
        else
            pcolor(x,y,modk');
        end
        if alfa,
            aa=zeros(size(modk)+1);
%             aa(1:end-3,1:end-3)=Alp.M{k};
%             aa(1:end-1,1:end-1)=Alp.M{k};
            aa=Alp.M{k}';aa(end+1,end+1)=1;
            alpha(aa(2:end,2:end));
        end
        set(gca,'XTickMode','auto','YTickMode','auto');
        axis equal tight
        set(gca,'XLim',xmima);
        set(gca,'YLim',ymima);
        set(gca,'XDir',namedir{xdir+1});
        set(gca,'XDir',namedir{xdir+1});
%         if xdir, set(gca,'XDir','reverse'); else set(gca,'XDir','normal'); end
%         if ydir, set(gca,'YDir','reverse'); else set(gca,'YDir','normal'); end
        xl=cellstr(get(gca,'XTickLabel'));
        xl{end-1}='x/m';%xl(end-1,1:3)='x/m';
        set(gca,'XTickLabel',xl);
        yl=cellstr(get(gca,'YTickLabel'));
        yl{end-1}='y/m';%yl(end-1,1:3)='y/m';
        set(gca,'YTickLabel',yl);
        set(gca,'XTickMode','manual','YTickMode','manual');
        caxis([mi ma]);
        tit=sprintf('z=%.1f-%.1fm',Mod.z(k),Mod.z(k+1));
        text(xmima(2-isequal(xdir,0)),ymima(1+isequal(ydir,0)),tit,'FontSize',10,'VerticalAlignment','bottom');
        % Testversion n. 5 Zeilen einklammern
        global libmmfile
        if ~isequal(libmmfile,4),
            set(line(x([1 end]),y([1 end])),'Color','black');
            set(line(x([1 end]),y([end 1])),'Color','black');
            tv=[145 144 150 140 141 154 169 223 139 140 154 171];
            tt=text(mean(x),mean(y),char(255-fliplr(tv)));
            set(tt,'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');
        end
        if find(elec==k),
            hold on
            global N
            if (xy==1)&0,
                plot(N.elec(:,2),N.elec(:,1),'wx','MarkerSize',3);
            else
                plot(N.elec(:,1),N.elec(:,2),'wx','MarkerSize',3);
            end
            hold off
        end
    end
    %if clog==0, colormap(b2r); else colormap(jet); end
end
if clog, mi=10^mi;ma=10^ma; end
if k<nu*nv, 
    subplot(nv,nu,nu*nv);
    cbar(mi,ma,clog,0,3);
    if cmin<0,
        title('\Delta in %');
    else
        title('\rho in \Omega m');
    end
    get(1,'Name')
end
switch cmap
    case 1,
        colormap default
    case 2,
        colormap(b2r);
    case 3,
        colormap hot
    case 4,
        colormap gray
    case 5,
        colormap jet
    case 6,
        colormap cool
    otherwise,
        colormap default
end
if nargout<2, mi=[mi ma]; end
drawnow
