function finex(fig,name,fontsize)

% FINEX - fine export of figures
% finex(filename)
% finex(figure,filename)
% filename is without extension

if nargin<3, fontsize=0; end
if nargin<2, if ischar(fig), name=0; end; end
if nargin<1, error('No filename present'); end
if ischar(fig), 
    fontsize=name;
    name=fig; 
end
if ~ishandle(fig), fig=gcf; end 
oldfs=get(gca,'FontSize');
if strfind(name,'.eps'),
    efile=name;
else
    efile=[name '.eps'];
end
xt=get(gca,'XTick');xtl=get(gca,'XTickLabel');
yt=get(gca,'YTick');ytl=get(gca,'YTickLabel');
if fontsize>0,
    ch=get(gcf,'Children');
    for i=1:length(ch),
        if isequal(get(ch(i),'Type'),'axes'),
            set(ch(i),'FontSize',fontsize);
        end
    end
    set(gca,'XTick',xt,'XTickLabel',xtl);
    set(gca,'YTick',yt,'YTickLabel',ytl);
    xl=get(gca,'XLabel');
    set(xl,'FontSize',fontsize);
    yl=get(gca,'YLabel');
    set(yl,'FontSize',fontsize);
    tit=get(gca,'Title');
    if ishandle(tit), set(tit,'FontSize',fontsize); end
end
print(gcf,'-depsc2',efile);
dos(['epstopdf ' efile]);
% print(gcf,'-dpng','-r150',strrep(efile,'.eps','.png'));
exportpng(strrep(efile,'.eps','.png'));
if fontsize>0,
    set(xl,'FontSize',oldfs);
    set(yl,'FontSize',oldfs);
    for i=1:length(ch),
        if isequal(get(ch(i),'Type'),'axes'),
            set(ch(i),'FontSize',oldfs);
        end
    end
end
