function zpar=zparam(N,anz,lolo)

% ZPARAM - Computes optimal z-levels from 1D-sensitivity
% z_par = zparam(N[,nlay[,maxdepth]]);
% N=structure of vectors of a,b,m,n..electrode numbers,
% k..konfiguration factors, r..apparent resistivities
% elec..Electrode positions
% nlay..number of layers, otherwise No. of separations
if nargin<3, lolo=0; end 
if nargin<2, anz=0; end
data=length(N.a);
amq=zeros(data,1);bmq=amq;anq=amq;bnq=amq;
for l = 1:data,
    if N.a(l)>0, aaa=N.elec(N.a(l),:); else aaa=[Inf 0]; end
    if N.m(l)>0, mmm=N.elec(N.m(l),:); else mmm=[-Inf 0]; end
    if N.b(l)>0, bbb=N.elec(N.b(l),:); else bbb=[-Inf 0]; end
    if N.n(l)>0, nnn=N.elec(N.n(l),:); else nnn=[Inf 0]; end
    amq(l)=sum((aaa-mmm).^2);    
    anq(l)=sum((aaa-nnn).^2);    
    bmq(l)=sum((bbb-mmm).^2);    
    bnq(l)=sum((bbb-nnn).^2);    
end
% amq=0.25*amq;bmq=0.25*bmq;anq=0.25*anq;bnq=0.25*bnq;
bnq(find(isnan(bnq)))=Inf;
am=sqrt(amq);an=sqrt(anq);bm=sqrt(bmq);bn=sqrt(bnq);
depth=0.1;
ss=0;
zi=[];
si=[];
sens=0;
scal=1./am-1./an-1./bm+1./bn;
while sens<0.99,
    zzq=repmat(4*depth.^2,data,1);
    ss=1./sqrt(amq+zzq)-1./sqrt(anq+zzq)-1./sqrt(bmq+zzq)+1./sqrt(bnq+zzq);
    %sens=sqrt(sum(((1-ss./scal).^2)/data));
    sens=sum(1-ss./scal)/data;
    %ss=am./sqrt(amq+zzq)-an./sqrt(anq+zzq)-bm./sqrt(bmq+zzq)+bn./sqrt(bnq+zzq);
    %sens=sqrt(sum((ss).^2)/data);
    zi=[zi depth];
    si=[si sens];
    depth=depth*1.1;
end
if nargout==0, plot(zi,si,'bx-'); end
if anz==0, 
    anz=length(unique(N.k));
    if anz>12, anz=10; end
end
if lolo==0,
    bound=linspace(0,0.95,anz+1);
    zpar=[0];
    for b = bound,
        nb=max(find(si<b));
        zpar=[zpar zi(nb)];
    end
else
%     zmi=zi(max(find(si<0.95/(anz+1))));
%     zma=zi(max(find(si<0.95)));
    zmi=lolo/(anz);zma=lolo;
    zpar=[0 exp(linspace(log(zmi),log(zma),anz))];
end
% koregel dz darf nicht abnehmen
ddzpar=diff(diff(zpar));
while find(ddzpar<-1e-3),
    ma=max(find(ddzpar<0));
    zpar(1:ma+2)=linspace(zpar(1),zpar(ma+2),ma+2);
    ddzpar=diff(diff(zpar));
end
% koregel dz darf sich nicht mehr als 1.5-fachen
dz=diff(zpar);
for k=1:length(dz)-1,
  if dz(k+1)>1.5*dz(k),
      zpar(k+2)=zpar(k+1)+1.5*dz(k);
      dz=diff(zpar); 
  end
end
% runden auf vernünftige Werte
for l = 2:length(zpar)
    mul=1;
    zb=zpar(l);
    while zb<30, zb=zb*10;mul=mul*10; end
    zb=round(zb)/mul;
    zpar(l)=zb;
end