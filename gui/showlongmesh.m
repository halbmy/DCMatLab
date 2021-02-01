mi=min(Mesh.node(:,1));
ma=max(Mesh.node(:,1));
n=4;
mal=struct('clog',1,'cauto',0,'cmin',interperc(Mesh.model,3),'cmax',interperc(Mesh.model,97),'cbar',0,'ylim',[-17 0]);
clf;
for i=1:n,
   mal.xlim=(ma-mi)/n*[i-1 i]+mi;
   subplot(n+1,1,i);
   tripatchmod(Mesh,Mesh.model,mal);
end
subplot(n+1,1,n+1);cbar(mal.cmin,mal.cmax,mal.clog);