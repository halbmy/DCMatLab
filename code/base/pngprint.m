function pngprint(fig,outfile)

% PNGPRINT - print figure into eps graphics file
% pngprint(figure_handle,filename)

if (nargin==1)&&ischar(fig), outfile=fig;fig=gcf; end
set(fig,'PaperPositionMode','manual');
set(fig,'Units',get(fig,'PaperUnits'));
po=get(fig,'PaperPosition');
newpo=po([4 3])-po([2 1]);
if isequal(get(fig,'PaperPositionMode'),'manual')&&(newpo(1)>0)&(newpo(2)>0),
    set(fig,'PaperSize',po([4 3])-po([2 1]));
    set(fig,'PaperPositionMode','auto');
end
print(fig,'-dpng',outfile);
