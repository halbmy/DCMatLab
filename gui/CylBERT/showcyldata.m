ndata=length(Data);
d0=readunifile('980i\980i.data');
% [cmin,cmax]=interperc(d0.r);
cmin=7;cmax=70;
nv=3;nu=floor(ndata/nv+0.5);
for i=1:ndata,
%     if  i>1, dd{i}=readunifile(['980i\' num2str(i) '.dat']); end
    Data{i}.r=Data{i}.rho.*d0.k;
    Data{i}.r(Data{i}.r<=0)=NaN;
    subplot(nu,nv,i);
    showpseudoflat(d0,Data{i}.r,cmin,cmax,0);
end
% epsprint(1,'daten-absolut7-70_Ohmm_log',1)
%%
figure(2);clf;
nv=3;nu=floor((ndata-1)/nv+0.5);
cmin=-80;cmax=0;
for i=2:ndata,
   subplot(nu,nv,i-1);
   di=(Data{i}.r./d0.r-1)*100;
   showpseudoflat(d0,di,cmin,cmax,0);
end
subplot(nu,nv,nu*nv);cbar(cmin,cmax,0);
% epsprint(2,'daten-diff_gegen0',1)
%%
figure(3);clf;
nv=3;nu=floor((ndata-1)/nv+0.5);
cmin=-30;cmax=10;
for i=2:ndata,
   subplot(nu,nv,i-1);
   di=(Data{i}.r./Data{i-1}.r-1)*100;
   showpseudoflat(d0,di,cmin,cmax,0);
end
subplot(nu,nv,nu*nv);cbar(cmin,cmax,0);
% epsprint(3,'daten-diff_stepwise',1)