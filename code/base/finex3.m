function finex(fig,name)

% FINEX - fine export of figures
% finex(filename)
% finex(figure,filename)
% filename is without extension

if nargin<1, error('No filename present'); end
if ischar(fig), name=fig; end
if ~ishandle(fig), fig=gcf; end 
oldfs=get(gca,'FontSize');
if strfind(name,'.eps'),
    efile=name;
else
    efile=[name '.eps'];
end
fontsize=14;
ch=get(gcf,'Children');
for i=1:length(ch),
    if isequal(get(ch(i),'Type'),'axes'),
        set(ch(i),'FontSize',fontsize);
    end
end
xl=get(gca,'XLabel');
set(xl,'FontSize',fontsize);
yl=get(gca,'YLabel');
set(yl,'FontSize',fontsize);
zl=get(gca,'ZLabel');
set(zl,'FontSize',fontsize);
print(gcf,'-depsc2',efile);
dos(['epstopdf ' efile]);
print(gcf,'-dpng','-r150',strrep(efile,'.eps','.png'));
set(xl,'FontSize',oldfs);
set(yl,'FontSize',oldfs);
set(zl,'FontSize',oldfs);
for i=1:length(ch),
    if isequal(get(ch(i),'Type'),'axes'),
        set(ch(i),'FontSize',oldfs);
    end
end
