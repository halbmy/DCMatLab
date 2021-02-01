function cbar(cmin,cmax,lolo,dir,anz,lab)

% CBAR - Draw single color bar (for exporting)
% cbar(cmin,cmax[,lolo,dir,num,label])
% cmin/cmax - minimum/maximum color
% lolo - logarithmic spacing
% dir - direction (0=horizontal,1=vertical)
% num - number of ticks
% label - text label

if nargin<6, lab='\rho in \Omega\cdotm'; end
if nargin<5, anz=0; end
if nargin<4, dir=0; end %horiz
if nargin<3, lolo=0; end
if nargin<1, cmin=0; end
if nargin<2, cmax=cmin+1; end
if lolo, cmax=log10(cmax);cmin=log10(cmin); end
n=size(colormap,1);
if anz==0, anz=11; end
if dir==0,
    image([0 n 0],[0 n/20],(0:n));
    axis tight
    set(gca,'DataAspectRatio',[1 1 1]);
    set(gca,'YTick',[],'XDir','normal');
    set(gca,'XTick',linspace(0,n,anz));
    yt=get(gca,'XTick');
    ytl=yt/n*(cmax-cmin)+cmin;
    if lolo==1, ytl=10.^ytl; end
    fi=find(ytl>=9.5);ytl(fi)=round(ytl(fi));
    fi=find(ytl<9.5);
    for i=1:length(fi),
        l=0;yy=ytl(fi(i));
        if yy~=0,
            while abs(yy)<9.5, yy=yy*10;l=l+1; end
            yy=round(yy);
            for ll=1:l, yy=yy/10; end
        end
        ytl(fi(i))=yy;
    end
%     if ~isempty(lab), ytl(end-1,1:length(lab))=lab; end
    set(gca,'XTickLabel',ytl);
%     po=get(gca,'Position');
%     po(4)=po(4)/2;
%     po(2)=po(2)*2;
%     set(gca,'Position',po);
else
    image([0 n/10],[1 n 1],(1:n)');
    set(gca,'DataAspectRatio',[1 1 1]);
    axis tight
    set(gca,'XTick',[],'YDir','normal');
    set(gca,'YTick',linspace(1,n,anz))
    yt=get(gca,'YTick');
    ytl=yt/n*(cmax-cmin)+cmin;
    if lolo==1, ytl=10.^ytl; end
    fi=find(ytl>=9.5);ytl(fi)=round(ytl(fi));
    fi=find(ytl<9.5);
    for i=1:length(fi),
        l=0;yy=ytl(fi(i));
        if yy~=0,
            while abs(yy)<9.5, yy=yy*10;l=l+1; end
            yy=round(yy);
            for ll=1:l, yy=yy/10; end
        end
        ytl(fi(i))=yy;
    end
%     if ~isempty(lab), ytl(end-1,1:length(lab))=lab; end
    set(gca,'YTickLabel',ytl,'YAxisLocation','right');
end
if cmin>0,
%     title('\rho in \Omega\cdotm');
title(lab);
end