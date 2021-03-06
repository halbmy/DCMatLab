function [cmin,cmax,iscb]=draw3dgridmodel(M,X,Y,Z,MAL,Bg,Alph)

% DRAW3DGRIDMODEL - Draw 3d Model M with X,Y,Z
% [cmin,cmax]=draw3dgridmodel(M,X,Y,Z,OPT,Bg,Alph)
% X,Y,Z - Coordinate Vektors
% M - Model to draw (length(X)-1,length(Y)-1,length(Z)-1)
% OPT structure of possible fields
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
% Bg - Background values (displayed in title) (length(Bg)=length(z))
% Alph  - use Alph values for alpha shading

if nargin<2, X=0:size(M,1); end
if nargin<3, Y=0:size(M,2); end
if nargin<4, Z=0:size(M,3); end
if nargin<5,
%     MAL=struct('cauto',1,'cmin',100,'cmax',500,'cmap',0,'log',0,'xdir',0,...
%         'xy',1,'nu',0,'nv',0,'elec',0,'startwith',1,'cont',0,'vie',0,'alpha',0);      
  MAL=struct('cauto',1);
end
if nargin<7,
    Alph=[];
end
if isfield(MAL,'xy'), xy=MAL.xy; else xy=1; end
if isfield(MAL,'xdir'), xdir=MAL.xdir; else xdir=0; end
if isfield(MAL,'ydir'), ydir=MAL.ydir; else ydir=0; end
if isfield(MAL,'cmax'), cmax=MAL.cmax; else cmax=1; end
if isfield(MAL,'cmin'), cmin=MAL.cmin; else cmin=0; end
if isfield(MAL,'cmap'), cmap=MAL.cmap; else cmap=0; end
if isfield(MAL,'cbar'), colbar=MAL.cbar; else colbar=0; end
if isfield(MAL,'log'), clog=MAL.log; else clog=0; end
if isfield(MAL,'clog'), clog=MAL.clog; else clog=0; end
if isfield(MAL,'cauto'), cauto=MAL.cauto; else cauto=1; end
if isfield(MAL,'nu'), nu=MAL.nu; else nu=0; end
if isfield(MAL,'nv'), nv=MAL.nv; else nv=0; end
if isfield(MAL,'elec'), elec=MAL.elec; else elec=0; end
if isfield(MAL,'startwith'), startwith=MAL.startwith; else startwith=1; end
if isfield(MAL,'cont'), cont=MAL.cont; else cont=[]; end
if isfield(MAL,'vie'), vie=MAL.vie; else vie=0; end
if isfield(MAL,'alpha'), alfa=MAL.alpha; else alfa=0; end
if isfield(MAL,'fs'), fs=MAL.fs; else fs=10; end

iscb=0; % no colorbar included
clf;
Im=length(X)-1; Jm=length(Y)-1; Km=length(Z)-1;
if cauto==1,
    if length(unique(M(:)))<10,
        cmin=min(M(:));
        cmax=max(M(:));
        if clog, cmin=log10(cmin);cmax=log10(cmax); end
        if cmin==cmax, cmin=cmax*0.99; end
    else
        if clog==1, 
            [NN,VV]=hist(log10(M(find(M~=NaN))),100); 
        else 
            [NN,VV]=hist(M(find(M~=NaN)),100); 
        end
        CN=cumsum(NN);CN=CN/max(CN);
        imin=max(find(CN<0.01));
        imax=min(find(CN>0.99));
        if isempty(imin), imin=1; end
        if isempty(imax), imax=length(VV); end
        cmin=VV(imin);
        cmax=VV(imax);
    end
    if cmax<=cmin, 
        cmax=max(M(:));
        cmin=min(M(:));
    end
    if cmax<=cmin, cmax=cmin*1.01;cmin=cmin*0.99; end
else
    if clog==1, 
        cmax=log10(cmax);
        cmin=log10(cmin);
    end
end
if alfa>0,
%     global S
%     if ~isempty(S)
%         Alph=M;
%         Alph(:)=log10(sum(abs(S)));
%         Alph=log10(Alph);
%         [nn,hh]=hist(Alph(:),50);
%         nnn=cumsum(nn)/prod(size(Alph));
%         mi=hh(min(find(nnn>0.1)));
%         ma=hh(max(find(nnn<0.7)));
%         mi=-1;ma=-0.5;
        mi=0.1;ma=1;
        Alph=(Alph-mi)/(ma-mi);
        Alph(find(Alph<0))=0;
        Alph(find(Alph>1))=1;
%     end
else
    Alph=[];
end
%Orientierung
mm=M;
if ndims(mm)<3, mm=reshape(mm,Im,Jm,Km); end
if clog==1,
    mm=log10(mm);
end
if vie==1, % xz
    yy=X;
    xx=Z;
    mm=permute(mm,[3 1 2]);
    Alph=permute(Alph,[3 1 2]);
    xlab='x';
    ylab='z';
    ydir=1;
elseif vie==2, %yz
    yy=Y;
    xx=Z;
    xlab='y';
    ylab='z';
    mm=permute(mm,[3 2 1]);
    Alph=permute(Alph,[3 2 1]);
    xdir=ydir;
    ydir=1;
else %xy
    if xy==1,
        xx=Y;
        yy=X;
        xlab='x';
        ylab='y';
        mm=permute(mm,[2 1 3]);
        Alph=permute(Alph,[2 1 3]);
    else
        xx=X;
        yy=Y;
        xlab='y';
        ylab='x';
    end
end
if xdir==0, xdir='normal'; else xdir='reverse'; end
if ydir==0, ydir='normal'; else ydir='reverse'; end

% Alle Scheiben malen
lk=size(mm,3);
if nu*nv==0,
    rel=(max(yy)-min(yy))/(max(xx)-min(xx));
    nv=fix(sqrt(lk*rel)+0.499);
    if nv>lk, nv=lk; end
    nu=round(lk/nv+0.499);
    while (nv-1)*nu>=lk, nv=nv-1; end
end
nr=0;
if clog, cont=log10(cont); end
while (nr<nu*nv)&(nr+startwith-1<lk),
    nr=nr+1;
    subplot(nv,nu,nr);
    nn=nr+startwith-1;
    pp=mm(:,:,nn);
    pp(end+1,:)=1;
    pp(:,end+1)=1;
    if isempty(Alph),
        pcolor(yy,xx,pp);
%         contourf(yy(1:end-1),xx(1:end-1),pp(1:end-1,1:end-1),30);
    else
        cc=Alph(:,:,nn);
        vv=version;
        if str2double(vv(1:3))<6.5,
            pp(:,end+1)=1;
            pp(2:end+1,:)=pp(1:end,:);
            pcolor([yy(1);yy(:)],[xx(1);xx(:)],pp);
            %pcolor(yy,xx,pp);
            %cc(:,2:end+1)=cc;
            cc(end+1,end+1)=0;
            alpha(cc);
        else
            pcolor(yy,xx,pp);
            cc(end+1,end+1)=1;
            alpha(cc(2:end,2:end));
        end
    end
    if (max(length(X),length(Y))>20)&(length(unique(pp(:)))>1), 
        shading flat; 
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
    caxis([cmin cmax]); 
    if ~isempty(cont),
        hold on
        contour(yy,xx,pp,cont,'k-')
        hold off
    end
    if (~isempty(elec))&(vie==0),
        if find(elec==nn),
            global N
            if ~isempty(N)&&isfield(N,'elec'),
                hold on
                if xy==0,
                    plot(N.elec(:,2),N.elec(:,1),'wx','MarkerSize',3);
                else
                    plot(N.elec(:,1),N.elec(:,2),'wx','MarkerSize',3);
                end
                hold off
            end
            %else
            %grid on
        end
        %    else
        %grid on
    end
    set(gca,'XTickMode','auto','YTickMode','auto');
    axis equal tight
    set(gca,'XDir',xdir,'YDir',ydir);
    if vie==0, tit=sprintf('z=%.1f-%.1fm',Z(nn),Z(1+nn)); end
    if vie==1, tit=sprintf('y=%.1f-%.1fm',Y(nn),Y(1+nn)); end
    if vie==2, tit=sprintf('x=%.1f-%.1fm',X(nn),X(1+nn)); end
%     if (nargin>5)&&(vie==0)&&(nn<=length(Bg))&&(Bg(nn)>0),
%         tit=sprintf('%s: %d\\Omegam',tit,round(Bg(nn)));  end
    %title(tit,'FontSize',10);
    xmima=[min(X) max(X)];ymima=[min(Y) max(Y)];
    lx=length(X);ly=length(Y);lz=length(Z);
    if  vie==1, xmima=[min(X) max(X)];ymima=[min(Z) max(Z)];
        hold on;set(line(X([1 lx lx 1 1]),Z([1 1 lz lz 1])),'Color','black');hold off;
    elseif vie==2, xmima=[min(Y) max(Y)];ymima=[min(Z) max(Z)];
        hold on;set(line(Y([1 ly ly 1 1]),Z([1 1 lz lz 1])),'Color','black');hold off; 
    else xmima=[min(X) max(X)];ymima=[min(Y) max(Y)]; 
        hold on;set(line(X([1 lx lx 1 1]),Y([1 1 ly ly 1])),'Color','black');hold off; end
    if fs, 
        set(gca,'FontSize',fs);
        text(xmima(2-isequal(xdir,'normal')),ymima(1+isequal(ydir,'normal')),tit,'FontSize',fs,'VerticalAlignment','bottom'); 
    end
    xl=get(gca,'XTickLabel');
    if ischar(xl), xl=cellstr(xl); end
    xl{end-1}=[xlab '/m'];
    set(gca,'XTickLabel',xl);
    yl=get(gca,'YTickLabel');
    if ischar(yl), yl=cellstr(yl); end
    yl{end-1}=[ylab '/m'];
    set(gca,'YTickLabel',yl);
    set(gca,'XTickMode','manual','YTickMode','manual');
    % Testversion n�chste 5 Z. einklammern
    global libmmfile
    if ~isequal(libmmfile,4)
        set(line(X([1 end]),Y([1 end])),'Color','black');
        set(line(X([1 end]),Y([end 1])),'Color','black');
        tv=[145 144 150 140 141 154 169 223 139 140 154 171];
        tt=text(mean(X),mean(Y),char(255-fliplr(tv)));
        set(tt,'FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle');
    end
    if colbar>0,
        if (vie==0)&(colbar==1),
            hc=colorbar;
            clabel=10.^str2num(get(hc,'YTickLabel'));
        else
            hc=colorbar('horiz');
            %get(hc,'XTick');
            %set(hc,'XTick',log10([100 200 400]));
            clabel=10.^str2num(get(hc,'XTickLabel'));
        end
        if(clog==1)
            fi=find(clabel>10);
            clabel(fi)=round(clabel(fi));
            fi=find(clabel<=10);
            clabel(fi)=round(10*clabel(fi))/10;
            if (vie==0)&(colbar==1),
                set(hc,'YTickLabel',num2str(clabel));
            else
                set(hc,'XTickLabel',num2str(clabel));
            end
        end
    end
end % of subplots
if clog, cmin=10^cmin;cmax=10^cmax; end
if nr<nu*nv, % one subplot left
    subplot(nv,nu,nu*nv);
    iscb=1;
    cbar(cmin,cmax,clog,0,3);
    if strcmp(get(gcf,'Name'),'Sensitivity'),
        title('sensitivity');
    else
        title('\rho in \Omega m');
    end
end
if nargout<2, cmin=[cmin cmax]; end
