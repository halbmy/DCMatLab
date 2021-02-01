clear;clc;
xmb=(10:0.5:30)';
topo=[0 0;15 0;20 4;25 4;40 0];
plot(topo(:,1),topo(:,2));
rp=4;
xz=mbm2xz(xmb,topo,1,rp);
hold on;plot(xz(:,1),xz(:,2),'x');hold off
hold on;plot(topo(rp,1),topo(rp,2),'ro');hold off
ylim([-1 6]);grid on;
di=sqrt(sum(diff(xz).^2,2))';