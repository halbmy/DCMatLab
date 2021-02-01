function [cmin,cmax]=showpseudoring(N,feld,cmin,cmax)

cla reset;
if nargin<2, feld=N.r; end
islog=(min(feld)>0);
if nargin<3, cmin=min(feld); end
if nargin<4, cmax=max(feld); end
if islog, feld=log10(feld);cmin=log10(cmin);cmax=log10(cmax); end
nel=size(N.elec,1);
xm=mean(N.elec(:,1));
ym=mean(N.elec(:,2));
phi1=atan2(N.elec(1,2)-ym,N.elec(1,1)-xm);
phi2=atan2(N.elec(2,2)-ym,N.elec(2,1)-xm);
di=(phi2>phi1)*2-1;
aa=double(N.a-1);bb=double(N.b-1);mm=double(N.m-1);nn=double(N.n-1);
% fi=find(bb<aa);bb(fi)=bb(fi)+nel;
% fi=find(mm<bb);mm(fi)=mm(fi)+nel;
% fi=find(nn<mm);nn(fi)=nn(fi)+nel;
% sep=mm-bb;%fi=find(sep<0);sep(fi)=sep(fi)+nel;
ab=(aa+bb)/2;ab(abs(aa-bb)>nel/2)=-0.5;
mn=(mm+nn)/2;mn(abs(mm-nn)>nel/2)=-0.5;
sep=abs(ab-mn)-1;
sep=min(sep,nel-2-sep);
% mid=mod(bb+sep/2,nel)+1;
mid=(bb+mm)/2;fi=find(abs(bb-mm)>nel/2);mid(fi)=mid(fi)-nel/2;
fi=find(mid<0);mid(fi)=mid(fi)+nel;
fi=find(mm-bb>(nel-2)/2);
mid(fi)=mod(aa(fi)-sep(fi)/2,nel)+1;
r=max(sep)-sep+1;
r=r;
phi=(mid-1)/nel*2*pi;
rr=sqrt((N.elec(:,1)-xm).^2+(N.elec(:,2)-ym).^2);
xx=r.*cos(phi1+phi);
yy=r.*sin(phi1+phi);
x=linspace(min(xx),max(xx),100);
y=linspace(min(yy),max(yy),100);
[X,Y]=meshgrid(x,y);
Z=griddata1(xx,yy,feld,X,Y);
XY=X.^2+Y.^2;Z(find(XY<0.8))=NaN;
imagesc(x,y,Z);
cmap=colormap(jet);
caxis([cmin-(cmax-cmin)/length(cmap)*2 cmax]);
cmap(1,:)=1;
colormap(cmap);
axis equal tight
set(gca,'YDir','normal','XAxisLocation','top');
% xlabel('relative x');
% ylabel('relative y');
hold on
plot(xx,yy,'k.','MarkerSize',1);
rrr=rr*max(r)/max(rr);
if 0,
    phiel=(0:nel-1)'*di/nel*2*pi;
    plot(rrr.*cos(phi1+phiel),rrr.*sin(phi1+phiel),'kx-');
else
    phiel=atan2(N.elec(:,2)-ym,N.elec(:,1)-xm);
    plot(rrr.*cos(phiel),rrr.*sin(phiel),'kx-');
end
plot(rrr(1).*cos(phi1),rrr(1).*sin(phi1),'ro');
hold off
hc=colorbar('horiz');
dar=get(hc,'DataAspectRatio');
set(hc,'DataAspectRatio',dar.*[1 32 1]);
% set(hc,'DataAspectRatio',[1 32 1]);
xl=get(hc,'XLim');xl(1)=xl(1)+diff(xl)/64;
set(hc,'XLim',xl);
xt=get(hc,'XTick');xtl=num2strcell(xt);
if islog,
    if min(feld(:))>1,
        xtl=num2strcell(round(10.^xt));
    else
        xtl=num2strcell(round(10.^xt*10)/10);
    end
    set(hc,'XTickMode','manual','XTickLabel',xtl);
end

    
function ce=num2strcell(vec)

for i=1:length(vec),
    ce{i}=num2str(vec(i));
end