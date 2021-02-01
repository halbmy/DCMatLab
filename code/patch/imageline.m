function imageline(tnew,nRHOA,cax,canot);

if nargin<4, canot='Ohmm'; end
if nargin<3, cax=interperc(nRHOA,[2 98]); end
islog=(min(nRHOA(:))>0);
if islog,
    imagesc(tnew,1:size(nRHOA,2),log10(nRHOA)');caxis(log10(cax));
else
    imagesc(tnew,1:size(nRHOA,2),nRHOA');caxis(cax);
end
alpha(1-isnan(nRHOA'));
datetick;axis tight;grid on;
% xt=get(gca,'XTick');xtl={};for i=1:length(xt), xtl{i}=datestr(xt,15); end
% set(gca,'XTickLabel',xtl);
xlabel('time');ylabel('channel');
cb=colorbar('horiz');
if islog, set(cb,'XTickLabel',rndig(10.^get(cb,'XTick'))); end
set(cb,'YTick',mean(get(cb,'YLim')),'YTickLabel',canot,'DataAspectRatio',[1 256 1]);
