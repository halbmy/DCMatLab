function filfak(s,f)

%FILFAK - Plot Filter Factors
% filfak(s,f) or
% filfak(s,lambda)

if (nargin<2), f=1; end
if(length(f)~=length(s)), f=s.^2./(s.^2+f(1)); end

r=max(find(s>10*eps));
ax=plotyy(1:r,f,1:r,cumsum(f))
xlabel('singular values','FontSize',12)
ylabel('filter factors','FontSize',12)
set(get(ax(2),'Ylabel'),'String','cumulative information','FontSize',12)
set(ax(2),'FontSize',12)
set(ax(1),'FontSize',12)
set(ax(2),'YTick',0:100:400)
grid on
hold on
tf=f./s;
plot(tf/max(tf),'r--')
hold off
le=legend('filter factors','transform function');
po=get(le,'Position');
po(2)=0.55;
set(le,'Position',po)
le=legend(ax(2),'cumulative information')
set(le,'Color',[1 1 1]);
%print(gcf,'-dpng','-r100','filfak.png');