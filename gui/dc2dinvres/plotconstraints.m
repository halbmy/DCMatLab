function plotconstraints(Mod,FIX,XX,ZZ,N)

if nargin<5, N=[]; end
if nargin<4, ZZ=[]; end
if nargin<3, XX=[]; end
if nargin<2, FIX=[]; end

xz=zeros(size(Mod.x));zz=Mod.z;
if isfield(N,'topo'),
    xz=interp1(N.topo(:,1),N.topo(:,2),Mod.x,'linear','extrap');
    zz=-Mod.z;
end

hold on
if isequal(size(FIX),size(Mod.M)),
    [i,j]=find(FIX==-1);
    plot((Mod.x(i)+Mod.x(i+1))/2,(zz(j)+zz(j+1))/2+xz(i),'wx');
    [i,j]=find(FIX>0);
    plot((Mod.x(i)+Mod.x(i+1))/2,(zz(j)+zz(j+1))/2+xz(i),'wo');
end
if isequal(size(XX)+[1 0],size(Mod.M)),
    [i,j]=find(XX>0);
    plot(Mod.x(i+1),(zz(j)+zz(j+1))/2+xz(i),'w+');
end
if isequal(size(ZZ)+[0 1],size(Mod.M)),
    [i,j]=find(ZZ>0);
    plot((Mod.x(i)+Mod.x(i+1))/2,zz(j+1)+xz(i),'w+');    
end
hold off
