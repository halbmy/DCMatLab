function Poly = circles2cyl(N,ztop,zbottom,dd,dz)

zs=unique(N.elec(:,3));
if nargin<2, ztop=max(zs)+(max(zs)-min(zs))/2; end
if nargin<3, zbottom=min(zs)-(max(zs)-min(zs))/2; end
if nargin<4, dd=1; end
if nargin<5, dz=median(diff(zs)); end

%%
nodes=[];
for i=1:length(zs),
    elec=N.elec(N.elec(:,3)==zs(i),:);
    nel=size(elec,1);
    if dd==2,
       me=mean(elec);
       rad=sqrt((elec(:,1)-me(1)).^2+(elec(:,2)-me(2)).^2);
       wi=atan2(elec(:,2)-me(2),elec(:,1)-me(1));
       dwi=mod(diff(wi([1:end 1])),pi)';
       fi=find(abs(dwi)>pi/2);dwi(fi)=dwi(fi)-pi*sign(dwi(fi));
       newwi=wi'+dwi/2;
       newrad=(rad+rad([2:end 1]))'/2;
       newel=reshape([elec(:,1:2)';newrad.*cos(newwi)+me(1);newrad.*sin(newwi)+me(2)],2,nel*2)';
    else
       newel=elec;
    end
    elec(:,3)=zs(i);
    newel(:,3)=zs(i);
    l=size(N.elec,1);
    nodes=[nodes;newel];
end
%%
Poly=[];
nel=nel*dd;
nodes(1:dd:end,4)=-99;
top=nodes(end-nel+1:end,:);top(:,3)=ztop;top(:,4)=1;
bottom=nodes(1:nel,:);bottom(:,3)=zbottom;bottom(:,4)=1;
Poly.node=[bottom;nodes;top];
l=0;
Poly.face={};
for j=1:length(zs)+1,
    for i=1:nel-1,
        l=l+1;
        Poly.face{l}=[1 2 nel+2 nel+1]+i-1+(j-1)*nel;
    end
    l=l+1;
    Poly.face{l}=[nel 1 nel+1 2*nel]+(j-1)*nel;
end
Poly.face{end+1}=1:nel;
Poly.face{end+1}=(1:nel)+(length(zs)+1)*nel;
Poly.region=[mean(elec(:,1:2)) mean(zs) 2 1e-4]; % parameter region
Poly.node(end+1,:)=[mean(elec(:,1:2)) mean(zs)-dz/2 -999]; % reference electrode (arbitrary)
Poly.node(end+1,:)=[mean(elec(:,1:2)) mean(zs)+dz/2 -1000]; % reference point for Neumann problem
