function [cmin,cmax,iscb]=plotprofiles(N,field,mal)

% PLOTPROFILES - Plot profiles of a 3D data set
% plotprofiles(N)
% where N is a data structure of electrodes (elec)
% a/b/m/n electrode numbers, k-factors(k) and resistivities(r)

iscb=0;
if (nargin<1)||(~isstruct(N)),
    error('Data structure required as input!');
end
if nargin<3, mal=struct('cauto',1); end
if (nargin<2)||isempty(field), field=N.r; end
nu=0;nv=0;
if ~isfield(mal,'cauto'), mal.cauto=1; end
if isfield(mal,'nu'), nu=mal.nu; end
if isfield(mal,'nv'), nv=mal.nv; end
clog=0;
if isfield(mal,'log'), clog=mal.log; end
if isfield(mal,'clog'), clog=mal.clog; end
if min(field)<0, clog=0; end
% 2d profiles available 
if mal.cauto==1,
  if clog,
    [NN,VV]=hist(log10(field),100);
  else
    [NN,VV]=hist(field,100);
  end
  CN=cumsum(NN);CN=CN/max(CN);
  imin=max(find(CN<0.01));
  if isempty(imin), imin=1; end
  imax=min(find(CN>0.99));
  if isempty(imax), imax=100; end
  cmin=VV(imin);
  cmax=VV(imax);
  if clog, cmin=10^cmin;cmax=10^cmax; end
  mal.cmin=cmin;mal.cmax=cmax;mal.cauto=0;
else
  cmin=mal.cmin;cmax=mal.cmax;
end
if isfield(N,'eind'),
    nplots=length(N.zweid);
    if nu*nv==0,
        nv=round(sqrt(nplots))+1;
        nu=round(nplots/nv+0.499);
        while((nu-1)*nv>=nplots), nu=nu-1; end
        while((nv-1)*nu>=nplots), nv=nv-1; end
    end
    if nu*nv>1, clf; end
    rmin=min(field);rmax=max(field);
    if rmax<=rmin, rmax=rmin*2;rmin=rmax/2; end
    amin=99999;amax=0;
    for n=1:nplots,
        abh=N.eind{n}(:,1);mi=min(abh);ma=max(abh);
        if mi<amin, amin=mi; end
        if ma>amax, amax=ma; end
    end
    amin=round(amin*10)/10;amax=round(amax*10)/10;
    amed=round(sqrt(amin*amax));
    if rmin<0, rmed=(rmin+rmax)/2; else rmed=sqrt(rmin*rmax); end
    if (abs(rmin)>0.85)&&(abs(rmax)>10), rmin=round(rmin);rmax=round(rmax);rmed=round(rmed); end
    xtl=num2strcell([rmin rmed rmax]);ytl=num2strcell([amin amed amax]);
    ytl{2}='AB/2';
    for n=1:nplots,
        if n<=nu*nv,
            if nu*nv>1, subplot(nu,nv,n); end
            abh=N.eind{n}(:,1);
            rhoa=field(N.nr{n});%N.eind{n}(:,2);
            if rmin>0, loglog(rhoa,abh,'x-'); else
                semilogy(rhoa,abh,'x-'); end
            grid on;
            set(gca,'xlim',[rmin rmax],'ylim',[amin amax]);
            set(gca,'YDir','reverse','XAxisLocation','top');
            set(gca,'XMinorGrid','off','YMinorGrid','off')
            if length(get(gca,'XTick')<3), set(gca,'XTick',[rmin rmed rmax],'XTickLabel',xtl); end
            if length(get(gca,'YTick')<3), set(gca,'YTick',[amin amed amax],'YTickLabel',ytl); end
            tt=text(rmax,amax,N.names{n});
            set(tt,'VerticalAlignment','bottom','HorizontalAlignment','right');
        end
    end
    return;
end
if ~isfield(N,'zweid'),
%   [cmin,cmax]=showdata3d(N,field,mal);
  return;
end
nplots=length(N.zweid);
nn=N.zweid{1};
if nu*nv==0,
    nv=round(sqrt(nplots)/size(nn.elec,1)*12)+1;
    nu=round(nplots/nv+0.499);
    while((nu-1)*nv>=nplots), nu=nu-1; end
    while((nv-1)*nu>=nplots), nv=nv-1; end
end
if nu*nv>1, clf; end
%mal.cbar=(nu*nv==nplots);
mal.cbar=0;mal.elec=0;
for n=1:nplots,
    if n<=nu*nv,
        if nu*nv>1, subplot(nu,nv,n); end
        NN=N.zweid{n};
        NN.r=field(N.nr{n});
%         [mids,konfs]=showdata2d(NN,NN.r,mal);
        showdata2d(NN,NN.r,mal);
        yl=get(gca,'Ylim');xl=get(gca,'Xlim');
        nn=strrep(strrep(N.names{n},'_',' '),' ','');
        tt=text(xl(2),yl(2),nn);
        set(tt,'HorizontalAlignment','Right','VerticalAlignment','Bottom');
    end
end
if nu*nv>n, % one subplot left
    iscb=1;
    subplot(nu,nv,nu*nv);
    hc=cbar(mal.cmin,mal.cmax,clog,0,5,1);
    xl=xlim;xl(1)=xl(1)+diff(xl)/64;
    set(gca,'XLim',xl);
    %     xt=get(hc,'XTick');if xt(1)==0, xt(1)=1; end
    if isfield(mal,'canot')&&isstr(mal.canot),
        tit=mal.canot;
    else
        if cmin<0,
            tit='\Delta in %';
        else
            tit='\rho_a in \Omega\cdotm';
        end
    end
    xl=get(gca,'Xlim');yl=get(gca,'Ylim');
    set(text(mean(xl),yl(1),tit),'VerticalAlignment','bottom','HorizontalAlignment','center');
end
if nargout<2, cmin=[cmin cmax]; end
