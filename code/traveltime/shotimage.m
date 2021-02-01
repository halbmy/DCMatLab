function [A,B] = shotimage(Shot,showrez,field)

% SHOTIMAGE - Plot apparent velocity image of a shot
% shotimage(Shot)
% [VA,REC] = shotimage(Shot)
% VA ... apparent velocity matrix (in m/s)
% REC .. reciprocity matrix (in %) 
% [...]=shotimage(Shot,1) shows reciprocity

if nargin<1, error('specify shot!'); end
A=ones(size(Shot.pos,1))*NaN;
su=0;
for i=1:length(Shot.ns),
    for j=1:length(Shot.nx{i}),
        dl=norm(Shot.pos(Shot.ns{i},:)-Shot.pos(Shot.nx{i}(j),:));
        if nargin>2, aa=field(su+j);
        else aa=dl/Shot.tt{i}(j)*1000; end
        A(Shot.ns{i},Shot.nx{i}(j))=aa;
    end
    su=su+length(Shot.nx{i});
end
B=abs((A-A')./(A+A')*200);
if (nargin>1)&&(showrez>0), 
    imagesc(B);
else
    imagesc(A); 
end
if nargin==11, colormap([1 1 1;jet(64)]);
else colormap(jet(64)); end
set(gca,'XTickMode','auto','YTickMode','auto','XTickLabelMode','auto','YTickLabelMode','auto');
axis ij equal tight;
xt=get(gca,'XTick');xtl=num2strcell(xt);
yt=get(gca,'YTick');ytl=num2strcell(yt);
ytl{end-1}='shot';xtl{end-1}='re';
set(gca,'XTickMode','manual','YTickMode','manual','XTickLabelMode','manual','YTickLabelMode','manual');
set(gca,'XTickLabel',xtl,'YTickLabel',ytl,'XAxisLocation','top');
cb=colorbar;
set(cb,'YTickLabel',num2strcell(rndig(get(cb,'YTick'),3)));
% imagesc(B);axis equal tight;colorbar
