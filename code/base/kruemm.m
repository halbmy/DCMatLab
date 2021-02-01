function cc=kruemm(x,y,t)

if nargin<3, t=1:length(x); end
xp=diff(x)./diff(t);
yp=diff(y)./diff(t);
xp(2:end-1)=(xp(1:end-2)+xp(3:end)+xp(2:end-1))/3;
yp(2:end-1)=(yp(1:end-2)+yp(3:end)+yp(2:end-1))/3;
tt=t(1:end-1)+diff(t)/2;
xpp=diff(xp)./diff(tt);
ypp=diff(yp)./diff(tt);
xp=(xp(1:end-1)+xp(2:end))/2;
yp=(yp(1:end-1)+yp(2:end))/2;
cc=(xp.*ypp-xpp.*yp)./((xp.^2+yp.^2).^1.5);
cc=[0;cc(:);0];