function mm=lkurv(rho,eta,ak,n)

%% LKURV - L-Curve Criterion
%%         Estimation of maximum curvature and plotting
%% index=lkurv(rho,eta[,alphas]) 

%global hinv
%cc=-(curvature(log(rho),log(eta),1,ak));
%mm=ecke(rho,eta);

if nargin<3, ak=[]; end
if nargin<4, n=0; end
if isempty(ak),
%   cc=-curvature(rho,eta,1);
    cc=curvature(rho/max(rho),eta/max(eta));
else
%   cc=curvature(rho/max(rho),eta/max(eta),n,ak);
  cc=curvature(rho,eta,n,ak);
end
mm=min(find(cc==max(cc)));
if (mm==1)|(mm==length(cc)-1),
    mm=1;flag=0;
    while((flag==0)&(mm+1<length(rho))),
        mm=mm+1;
        if((cc(mm)>cc(mm-1))&(cc(mm)>cc(mm+1))) flag=1; end
    end
    if flag==0, mm=min(find(cc==max(cc))); end
end
if nargout==0,
  figure(9);
  plot(rho,eta,'bx-',rho(mm),eta(mm),'ro',rho,abs(cc)/max(abs(cc))*max(eta)*0.9,'g+');
  le=legend('L-Curve','Optimized \lambda','Curvature');
  set(le,'FontSize',14);
  xlabel('Solution norm ||D(S\Delta m - \Delta d)||','Fontsize',14);
  ylabel('Model norm ||C \Delta m||','FontSize',14)
  fac=1;ak=round(ak*100)/100;
  if ~isempty(ak),
      text(rho(1)*1.02,eta(1)*0.95,strcat('\lambda=',num2str(ak(1),'%4g')),'FontSize',12)
      text(rho(mm)*1.01,eta(mm)*1.05,strcat('\lambda=',num2str(ak(mm),'%4g')),'FontSize',12)
      text(rho(end)*1.01,eta(end)*1.05,strcat('\lambda=',num2str(ak(end),'%4g')),'FontSize',12)
  end
  %exportfig(gcf,'l-kurve','Format','png','Resolution',300,'Color','rgb');
end

function curv=curvature(x,y,n,t)
% Curvature of parametric function x(t), y(t)
% curv=curvature(x,y,n)

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
      xpp=xpp./(t(3:nn)-t(1:nn-2))/2;
      ypp=(y(3:nn)-y(2:nn-1))./(t(3:nn)-t(2:nn-1));
      xpp=ypp-(y(2:nn-1)-y(1:nn-2))./(t(2:nn-1)-t(1:nn-2));
      xpp=ypp./(t(3:nn)-t(1:nn-2))/2;
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
curv=reshape(curv,size(x));
return;
