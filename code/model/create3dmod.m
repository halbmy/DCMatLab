function Mod=create3dmod(N,nz,dx,dy)

% CREATEMOD - Create model
% Model = create3dmod(N[,nz,dx,dy])
% create parametric model from data struct N
% nz(number of layers),
% dx/dy (smallest spacing) are optional

if nargin<2, nz=[]; end
if nargin<3, dx=[]; end
if nargin<4, dy=dx; end

if isempty(nz),
    nz=length(unique(sort(round(N.k))));
    if nz>8, nz=8; end
end
Mod.z=zparam(N,nz);

aa=diff(unique(sort(N.elec(:,1))));
dxel=min(aa(find(aa)));
aa=diff(unique(sort(N.elec(:,2))));
dyel=min(aa(find(aa)));
minx=min(N.elec(:,1))-dxel;
maxx=max(N.elec(:,1))+dyel;
miny=min(N.elec(:,2))-dyel;
maxy=max(N.elec(:,2))+dyel;
Mod.x0=minx;
Mod.y0=miny;
if isempty(dy), dy=dyel; end %or /2 ?
if isempty(dx), dx=dxel; end
dx=round(dx*20)/20;
dy=round(dy*20)/20;
if nargin<3,
    refx=round(dx/2/(Mod.z(2)-Mod.z(1)))*2;
else
    refx=1;
end
if refx<1, refx=1; end
Mod.dx=dx/refx;
if nargin<4,
    refy=round(dx/2/(Mod.z(2)-Mod.z(1)))*2;
else
    refy=1;
end
if refy<1, refy=1; end
Mod.dy=dy/refy;

nnx=floor((maxx-minx)/Mod.dx);
nny=floor((maxy-miny)/Mod.dy);

%message(sprintf('Min(x)=%g Max(x)=%g Min(y)=%g Max(y)=%g D=%g',minx,maxx,miny,maxy,D)); 
zusatz=1;
K=length(Mod.z)-1;
% Model estimation
minKonf=min(abs(N.k)); 
% Sigma_q=Mittelwert der geringsten Eindringtiefe o. Mittel
rq=median(N.r(find(abs(N.k)==minKonf)));
%rq=mean(N.r);
if rq>20, rq=round(rq/10)*10; end
Mod.Bg=ones(size(Mod.z))*rq;
mode=1;
switch mode,
case 2, % ascending
    Mod.nx=(1:K)+refx-1;Mod.ny=(1:K)+refy-1;
case 1, % complicating
    Mod.nx=round(min(diff(Mod.z)/Mod.dx+1,Mod.z(2:end)/Mod.dx));
    Mod.ny=round(min(diff(Mod.z)/Mod.dy+1,Mod.z(2:end)/Mod.dy));
    Mod.nx=max(Mod.nx,round(refx.*((0:K-1)*0.13+1)));
    Mod.ny=max(Mod.ny,round(refy.*((0:K-1)*0.13+1)));    
    for k=1:length(Mod.nx)-1,
        if Mod.nx(k+1)>Mod.nx(k)+1, Mod.nx(k+1)=Mod.nx(k)+1; end
        if Mod.ny(k+1)>Mod.ny(k)+1, Mod.ny(k+1)=Mod.ny(k)+1; end
    end
case 3, % quadratic
    Mod.nx=round(diff(Mod.z)/Mod.dx);
    Mod.ny=round(diff(Mod.z)/Mod.dx);
otherwise, %logarithmically ascending
    Mod.nx=round(refx.*((0:K-1)*0.13+1));
    Mod.ny=round(refx.*((0:K-1)*0.13+1));
end
Mod.ncells=0;
for k=1:K,
    Mod.M{k}=ones(floor(nnx/Mod.nx(k)),floor(nny/Mod.ny(k)))*Mod.Bg(k);
    Mod.ncells=Mod.ncells+prod(size(Mod.M{k}));
end
message(sprintf('Creating para model: %d cells, rho=%.1f',Mod.ncells,rq));
message(sprintf('dx=%.1f (%.1f-%.1f) dy=%.1f (%.1f-%.1f)',Mod.dx,Mod.nx(1)*Mod.dx,...
    Mod.nx(end)*Mod.dx,Mod.dy,Mod.ny(1)*Mod.dx,Mod.ny(end)*Mod.dy));