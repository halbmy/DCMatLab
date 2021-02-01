function S=modsens(Mod,N)

% MODSENS - Integrate Model Sensitivity
% S = modsens(Mod,N)
% Mod - Model structure with z,x0,dx,...
% N - Data structure with a,b,m,n,r,k,elec

if ~iscell(Mod.M),
    S=calcsens3d(Mod.x,Mod.y,Mod.z,N);
    return;
end    
K=length(Mod.z)-1;
for k=1:K,
    kmod(k)=prod(size(Mod.M{k}));        
end
nmod=sum(kmod);
data=length(N.r);
nobreit=0;
S=zeros(data,nmod);
n=5;
[p,w]=gauleg(0,1,n);
nx=size(Mod.M{1},1);ny=size(Mod.M{1},2);
xmax=Mod.x0+nx*Mod.dx*Mod.nx(1);
ymax=Mod.y0+ny*Mod.dy*Mod.ny(1);
nnx=Mod.nx(1)*size(Mod.M{1},1);
nny=Mod.ny(1)*size(Mod.M{1},2);
start=0;
dz=diff(Mod.z);
wb=waitbar(0,'Integrating Sensitivity matrix...');
for k=1:K,
    modk=Mod.M{k};
    nx=size(modk,1);ny=size(modk,2);
    rx=floor(mod(nnx,nx*Mod.nx(k))/2);
    ry=floor(mod(nny,ny*Mod.ny(k))/2);
    x=((0:nx)*Mod.nx(k)+rx)*Mod.dx+Mod.x0;
    y=((0:ny)*Mod.ny(k)+ry)*Mod.dy+Mod.y0;
    if ~nobreit,
        x(1)=Mod.x0;x(end)=xmax;
        y(1)=Mod.y0;y(end)=ymax;
    end
    aller=fix(data/10);  % show waitbar every aller times
    mal=aller;
    for l = 1:data,
        sens=sens3dplane(x,y,Mod.z(k),Mod.z(k+1),N.elec(N.a(l),:),N.elec(N.m(l),:));
        if N.n(l)>0, 
            sens=sens-...
           sens3dplane(x,y,Mod.z(k),Mod.z(k+1),N.elec(N.a(l),:),N.elec(N.n(l),:)); end
        if N.b(l)>0, sens=sens-...
           sens3dplane(x,y,Mod.z(k),Mod.z(k+1),N.elec(N.b(l),:),N.elec(N.m(l),:));
           if N.n(l)>0, sens=sens+...
              sens3dplane(x,y,Mod.z(k),Mod.z(k+1),N.elec(N.b(l),:),N.elec(N.n(l),:)); end
        end
        S(l,start+1:start+kmod(k))=sens'*N.k(l);
        mal=mal-1;
        if mal<1,
            waitbar((start+l/data*kmod(k))/nmod,wb);
            mal=aller;
        end
    end
    start=start+kmod(k);
end
close(wb);  
