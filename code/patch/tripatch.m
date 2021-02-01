Mesh=loadmesh(meshname);
if exist(vecname)==2,
    att=load(vecname);
else
    patch('Vertices',Mesh.node,'Faces',Mesh.cell,'FaceColor',[1 1 1]);
    axis equal tight
    return;
    fprintf(['file ' vecname ' does not exist!']);
    att=1;
end
if ~exist('perc'), perc=5; end
di=max(Mesh.node)-min(Mesh.node);
po=get(gcf,'Position');
po(1)=50;po(3)=800;po(4)=di(2)/di(1)*1.3*800+50;
while po(4)>800, po(3:4)=po(3:4)/2; end
set(gcf,'Position',po);
if length(att)<Mesh.ncells, att(Mesh.ncells)=1; end
if length(att)>Mesh.ncells, att(Mesh.ncells+1:end)=[]; end
islog=(min(att)>0);
if islog, att=log10(att); end
if ~exist('cmin'), 
  [N,X]=hist(att,100);
  C=cumsum(N)/sum(N);
  cmin=X(min(find(C>perc/100)));%min(att);
end
if ~exist('cmax'),
  [N,X]=hist(att,100);
  C=cumsum(N)/sum(N);
  cmax=X(max(find(C<1-perc/100)));%max(att);
end
if ~(cmax>cmin), cmax=cmin+1; end
if islog, colormap(jet); else
  colormap(b2r);
  cmax=max(abs([cmin cmax]));cmin=-cmax;
end
%if (length(att)<100)|(length(unique(att))<10), emap(:)=0; end
amap=ones(size(att));
if exist('sensCov.vector','file'),
  senscov=log(load('sensCov.vector','-ascii'));
  [nn,hh]=hist(senscov,50);
  nnn=cumsum(nn)/length(senscov);
  mi=hh(min(find(nnn>0.02)));
  ma=hh(max(find(nnn<0.5)));
  amap=(senscov-mi)/(ma-mi);
  amap(find(amap<0))=0;
  amap(find(amap>1))=1;
end
att(find(att<cmin))=cmin;
att(find(att>cmax))=cmax;
patch2dmesh(Mesh,att,amap);
caxis([cmin cmax]);
xtl=cellstr(get(gca,'XTickLabel'));
xtl{end-1}='x/m';
set(gca,'XTickLabel',xtl);
ytl=cellstr(get(gca,'YTickLabel'));
ytl{end-1}='z/m';
set(gca,'YTickLabelMode','manual','YTickLabel',ytl);
cb=colorbar('horiz');
set(cb,'DataAspectRatio',[1 64 1]);
if islog,
  xt=get(cb,'XTick');
  xt=rndig(10.^xt,3);
  %fi=find(xt>1);xt(fi)=round(xt(fi)*10)/10;
  %fi=find(xt>10);xt(fi)=round(xt(fi));
  set(cb,'XTickLabel',num2strcell(xt));
end
set(cb,'YTick',mean(get(cb,'Ylim')),'YTickLabel','Ohmm');
elec=[];
if exist('datafile','file'),
  N=readinv2dfile(datafile);
  elec=N.elec;%[N.elec;N.elec(1,:)];
else
  rad=max(abs([xlim ylim]));
  if rad<1, %obviously a tree
    nel=24;
    elec=rad*sin((0:nel)'/nel*2*pi);
    elec(:,2)=rad*cos((0:nel)'/nel*2*pi);
  end
end
if ~isempty(elec),
    hold on;plot(elec(:,1),elec(:,2),'k.-','MarkerSize',1);hold off
end
