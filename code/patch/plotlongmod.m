function plotlongmod(Mod,MAL)

% PLOTLONGMOD - Plot long model using serveral stripes (subfigures)
% plotlongmod(Mod[,OPT]);

if nargin<2, MAL=struct('clog',1,'cauto',1); end

% Mod=modelimport2d('c:\halbmy\src\dcmatlab\dc2dinvres\rollalong\hafen8t1rollalong.mod');
% Mod=modelimport2d('c:\halbmy\src\dcmatlab\dc2dinvres\rollalong\teil22a-roll.mod');
% Mod=modelimport2d('c:\halbmy\src\dcmatlab\dc2dinvres\rollalong\teil22c.mod');
% Mod=modelimport2d('c:\halbmy\2d\diverse\ggl\Test_roll_along_z01.mod');
minx=min(Mod.x);maxx=max(Mod.x);
n=round(sqrt((maxx-minx)/max(Mod.z))/1.5);
mal=MAL;mal.high=1;
if ~isfield(mal,'clog'),
    if isfield(mal,'log'), mal.clog=mal.log; else mal.clog=1; end
end
if ~isfield(mal,'cauto')||(mal.cauto>0),
    mal.cauto=0;
    if mal.clog, mm=10.^interperc(log10(Mod.M(:)),[5 95]);
    else mm=interperc(Mod.M(:),[5 95]); end
    mal.cmin=mm(1);mal.cmax=mm(2);
    if mal.cmin==mal.cmax, mal.cmin=0.99*mal.cmax; end
end
mal.cbar=0;
di=(maxx-minx)/n/8;
% figure(3);
clf;
for i=1:n,
    mi=(maxx-minx-di)*(i-1)/n+minx;
    ma=(maxx-minx-di)*i/n+minx+di;
%     fprintf('%d\t%.1f-%.1f=%.1f\n',i,ma,mi,ma-mi);
    imi=max(find(Mod.x<=mi));if isempty(imi), imi=1; end
    ima=min(find(Mod.x>=ma));if isempty(ima), ima=length(Mod.x); end    
    xx=Mod.x(imi:ima);
    if i==1, len=max(xx)-min(xx); else
        if max(xx)-min(xx)>len, xx(1)=xx(1)+max(xx)-min(xx)-len; end
    end
%     xx(end)-xx(1)
    olddx=max(xx)-min(xx);
    MM=Mod.M(imi:ima-1,:);
    subplot(n+1,1,i);
    patch2dmodel(xx,Mod.z,MM,mal);
    set(gca,'XTickMode','manual','XTickLabelMode','manual');
end
subplot(n+1,1,n+1);
cbar(mal.cmin,mal.cmax,mal.clog,0,9);
set(gca,'DataAspectRatio',get(gca,'DataAspectRatio').*[1 2 1]);
title('\rho in \Omegam');