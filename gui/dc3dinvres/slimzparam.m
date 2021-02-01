function zpar=zparam(N,anz,zmax)

% ZPARAM - Computes optimal z-levels from 1D-sensitivity
% z_par = zparam(N[,nint]);
% N=struct of arrays a,b,m,n(electrode numbers), r(rhoa), k(konf)
%             elec..Coordinates of Electrodes
% nint..number of intervals(8)

if nargin<2, anz=0; end
if nargin<3, zmax=0; end
if anz==0,
    anz=length(unique(N.k));
    if anz>8, anz=8; end
end

amq=[];bmq=amq;anq=amq;bnq=amq;
data=length(N.a);
fb=find(N.b);fn=find(N.n);fbn=find(N.b.*N.n);
amq=sum((N.elec(N.a,:)-N.elec(N.m,:)).^2,2);
anq=amq;bmq=amq;bnq=amq;
bmq(fb)=sum((N.elec(N.b(fb),:)-N.elec(N.m(fb),:)).^2,2);
anq(fn)=sum((N.elec(N.a(fn),:)-N.elec(N.n(fn),:)).^2,2);
bmq(fbn)=sum((N.elec(N.b(fbn),:)-N.elec(N.n(fbn),:)).^2,2);
% for l = 1:data,
%     aaa=N.elec(N.a(l),:);
%     mmm=N.elec(N.m(l),:);
%     if N.b(l)>0, bbb=N.elec(N.b(l),:); else bbb=[Inf Inf 0]; end
%     if N.n(l)>0, nnn=N.elec(N.n(l),:); else nnn=[Inf Inf 0]; end
%     amq=[amq sum((aaa-mmm).^2)];    
%     anq=[anq sum((aaa-nnn).^2)];    
%     bmq=[bmq sum((bbb-mmm).^2)];    
%     bnq=[bnq sum((bbb-nnn).^2)];    
% end
%amq=0.25*amq;bmq=0.25*bmq;anq=0.25*anq;bnq=0.25*bnq;
bnq(find(isnan(bnq)))=Inf;
am=sqrt(amq);an=sqrt(anq);bm=sqrt(bmq);bn=sqrt(bnq);
depth=0.1;
ss=0;
zi=[];
si=[];
sens=0;
scal=1./am;scal(fn)=scal(fn)-1./an(fn);
scal(fb)=scal(fb)-1./bm(fb);scal(fbn)=scal(fbn)+1./bn(fbn);
while sens<0.99,
    %zzq=repmat(4*depth.^2,data,1);
    zq=4*depth^2;
    ss=1./sqrt(amq+zq);
    ss(fn)=ss(fn)-1./sqrt(anq(fn)+zq);
    ss(fb)=ss(fb)-1./sqrt(bmq(fb)+zq);
    ss(fbn)=ss(fbn)+1./sqrt(bnq(fbn)+zq);
    %sens=sqrt(sum(((1-ss./scal).^2)/data));
    sens=sum(1-ss./scal)/data;
    %ss=am./sqrt(amq+zzq)-an./sqrt(anq+zzq)-bm./sqrt(bmq+zzq)+bn./sqrt(bnq+zzq);
    %sens=sqrt(sum((ss).^2)/data);
    zi=[zi depth];
    si=[si sens];
    depth=depth*1.1;
end
if zmax==0, % by 1d-sens
    bound=linspace(0,0.95,anz+1);
    zpar=[0];
    for b = bound,
        nb=max(find(si<b));
        zpar=[zpar zi(nb)];
    end
else % logspaced
    zmi=zmax/anz;zma=zmax;
    zpar=[0 exp(linspace(log(zmi),log(zma),anz))];
end
    
% koregel dz darf nicht abnehmen
ddzpar=diff(diff(zpar));
while find(ddzpar<-1e-3),
    ma=max(find(ddzpar<0));
    zpar(1:ma+2)=linspace(zpar(1),zpar(ma+2),ma+2);
    ddzpar=diff(diff(zpar));
end
% koregel dz darf sich nicht mehr als verdoppeln
dz=diff(zpar);
for k=1:length(dz)-1,
  if dz(k+1)>2*dz(k), zpar(k+2)=zpar(k+1)+dz(k);
      dz=diff(zpar); 
  end
end
for l = 2:length(zpar)
    mul=1;
    zb=zpar(l);
    while zb<50, zb=zb*10;mul=mul*10; end
    zb=round(zb)/mul;
    zpar(l)=zb;
end