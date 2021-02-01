function curv=curvature(x,y,n,t)
% Curvature of parametric function x(t), y(t)
% curv=curvature(x,y,n)
% n-number of sampling points
nn=min([length(x) length(y)]);
x=x(1:nn);
y=y(1:nn);
if nargin<4, t=reshape(1:nn,size(x)); end
if nargin<3, n=0; end
vor=n;nach=n;
if n==0,
    ys=(y(3:nn)-y(1:nn-2))./(x(3:nn)-x(1:nn-2));
    dydx=diff(y)./diff(x);
    yss=diff(dydx)./(x(3:nn)-x(1:nn-2))*2;
    curv=yss./((1+ys.^2).^1.5);
    curv=[curv(1);curv(:);curv(length(curv))];
else
  if n==1,
      xp=(x(3:nn)-x(1:nn-2))./(t(3:nn)-t(1:nn-2));
      yp=(y(3:nn)-y(1:nn-2))./(t(3:nn)-t(1:nn-2));
      %xpp=x(1:nn-2)+x(3:nn)-2*x(2:nn-1);
      %ypp=y(1:nn-2)+y(3:nn)-2*y(2:nn-1);
      xpp=(x(3:nn)-x(2:nn-1))./(t(3:nn)-t(2:nn-1));
      xpp=xpp-(x(2:nn-1)-x(1:nn-2))./(t(2:nn-1)-t(1:nn-2));
      xpp=xpp./(t(3:nn)-t(1:nn-2))*2;
      ypp=(y(3:nn)-y(2:nn-1))./(t(3:nn)-t(2:nn-1));
      ypp=ypp-(y(2:nn-1)-y(1:nn-2))./(t(2:nn-1)-t(1:nn-2));
      ypp=ypp./(t(3:nn)-t(1:nn-2))*2;
  end
  if n==2,
      xp=(x(1:nn-4)-8*x(2:nn-3)+8*x(4:nn-1)-x(5:nn))/12;
      yp=(y(1:nn-4)-8*y(2:nn-3)+8*y(4:nn-1)-y(5:nn))/12;
      xpp=(-x(1:nn-4)+16*x(2:nn-3)-30*x(3:nn-2)+16*x(4:nn-1)-x(5:nn))/12;
      ypp=(-y(1:nn-4)+16*y(2:nn-3)-30*y(3:nn-2)+16*y(4:nn-1)-y(5:nn))/12;
  end
  if n==3,
      xp=-x(1:nn-3)+x(2:nn-2);
      yp=-y(1:nn-3)+y(2:nn-2);    
      xpp=2*x(1:nn-3)-5*x(2:nn-2)+4*x(3:nn-1)-x(4:nn);
      ypp=2*y(1:nn-3)-5*y(2:nn-2)+4*y(3:nn-1)-y(4:nn);
      vor=0;nach=3;
  end
  if n==4,
      xp=(-3*x(2:nn-5)-10*x(3:nn-4)+18*x(4:nn-3)-6*x(5:nn-2)+x(6:nn-1))/12;
      yp=(-3*y(2:nn-5)-10*y(3:nn-4)+18*y(4:nn-3)-6*y(5:nn-2)+y(6:nn-1))/12;
      xpp=(-13*x(1:nn-6)+228*x(2:nn-5)-420*x(3:nn-4)+200*x(4:nn-3)+15*x(5:nn-2)-12*x(6:nn-1)+2*x(7:nn))/180;    
      ypp=(-13*y(1:nn-6)+228*y(2:nn-5)-420*y(3:nn-4)+200*y(4:nn-3)+15*y(5:nn-2)-12*y(6:nn-1)+2*y(7:nn))/180;    
      vor=2;nach=4;
  end
  curv=(xp(:).*ypp(:)-xpp(:).*yp(:))./((xp(:).^2+yp(:).^2).^1.5);
  curv=[curv(1)*ones(vor,1);curv(:);curv(length(curv))*ones(nach,1)];
end
curv=-reshape(curv,size(x));
return;
