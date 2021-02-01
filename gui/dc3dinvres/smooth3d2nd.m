function C=smooth3d2nd(x,y,z,rb,nor,az)

% SMOOTH3D2nd - Second order smoothness matrix for 3d models
% C = smooth3d2nd(x,y,z[,bc[,nor[,alfaz]]])
% x,y,z - vectors of grid lines
% bc    - boundary condition index (default=0)
%        (0=dirichlet, 1=neumann, 2=modified neumann)
% nor   - normalize vectors (default=0)
% alfaz - weighting for vertical part

if nargin<6, az=1; end
if nargin<5, nor=0; end
if nargin<4, rb=0; end
if nargin<2, % for testing
    x=[0 1 2 4];%x=0:3;
    y=0:4;z=0:2;
end
if nor, x=1:length(x);y=1:length(y);z=1:length(z); end
% midpoints representing model cells
xq=(x(1:end-1)+x(2:end))/2;
yq=(y(1:end-1)+y(2:end))/2;
zq=(z(1:end-1)+z(2:end))/2;
% scaling because of comparability of regularization parameter
xmed=median([diff(xq(:));diff(yq(:))]);
xq=xq/xmed;
yq=yq/xmed;
zq=zq/xmed;
dx=diff(xq(:));
dx=[dx(1);dx;dx(end)];
lx=length(dx)-1;
dy=diff(yq(:));
dy=[dy(1);dy;dy(end)];
ly=length(dy)-1;
dz=diff(zq(:));
dz=[dz(1);dz;dz(end)];
lz=length(dz)-1;
[DX,DY,DZ]=ndgrid(dx,dy,dz);
DX(:,end,:)=[];DX(:,:,end)=[];
DY(end,:,:)=[];DY(:,:,end)=[];
DZ(end,:,:)=[];DZ(:,end,:)=[];
DXQ=DX(1:end-1,:,:)+DX(2:end,:,:);
DYQ=DY(:,1:end-1,:)+DY(:,2:end,:);
DZQ=DZ(:,:,1:end-1)+DZ(:,:,2:end);
im=2./DX(1:end-1,:,:)./DXQ;
ip=2./DX(2:end,:,:)./DXQ;
jm=2./DY(:,1:end-1,:)./DYQ;
jp=2./DY(:,2:end,:)./DYQ;
km=az*2./DZ(:,:,1:end-1)./DZQ;
kp=az*2./DZ(:,:,2:end)./DZQ;
if rb==0, ik=-im-ip-jm-jp-km-kp; end
im(1,:,:)=0;
ip(end,:,:)=0;
jm(:,1,:)=0;
jp(:,end,:)=0;
km(:,:,1)=0;
kp(:,:,end)=0;
if rb>1,
    im(end,:,:)=0;
    ip(1,:,:)=0;
    jm(:,end,:)=0;
    jp(:,1,:)=0;
    km(:,:,end)=0;
    kp(:,:,1)=0;
end
if rb>0, ik=-im-ip-jm-jp-km-kp; end
lxy=lx*ly;lxyz=lxy*lz;
C=spdiags([km(:) jm(:) im(:) ik(:) ip(:) jp(:) kp(:)],...
    [lxy lx 1 0 -1 -lx -lxy],lxyz,lxyz)';