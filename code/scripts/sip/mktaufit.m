t=logspace(-3.5,0,10);
f=logspace(-1,3,100);
clf;set(gca,'FontSize',14);
for i=1:length(t),
    wt=2*pi*f*t(i);
    ri=wt./(wt.^2+1);
    rr=wt.^2./(wt.^2+1);    
    semilogx(f,ri,'-');
    th=text(1/t(i)/2/pi,0.5,num2str(rndig(t(i)),'%g'));
    set(th,'FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','center');
%     loglog(f,atan(ri./rr),'x-');
    hold on;
end
set(text(f(1),0.5,' \tau='),'FontSize',14,'VerticalAlignment','bottom','HorizontalAlignment','center');
hold off;
grid on;
xlabel('f in Hz');
ylabel('amplitude');
set(gca,'YTick',[]);