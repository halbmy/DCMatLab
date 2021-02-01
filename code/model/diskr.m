function C=diskr(x,y,z,sigma,type,bound)

% DISKR - apply discretization scheme to obtain matrix C
% C = diskr(x,y,z,sigma,type,bound)
% x,y,z - grid node positions
% sigmac - cell conductivity
% type - discretization scheme (0-default)
%    0 = Dey & Morrison (volume discr.) (1979)
%    1 = Brewitt Taylor & Weaver (1976)
%    2 = Zhang, Mackie & Madden (1995)
%    3 = Wurmstich & Morgan (1994)
%    4 = Dey & Morrison 1 (point discr.)
% bound - use mixed boundary conditions(1) or dirichlet(0)

% if nargin==0, testfd;return; end
if nargin<4, error('Four input arguments are required!'); end
if nargin<5, type=0; end % default
if nargin<6, bound=0; end % default

I=length(x);
J=length(y);
K=length(z);
IJK=I*J*K;
%C=spalloc(IJK,IJK,7*IJK);
dx=diff(x);dx=[dx(1);dx(:);dx(end)];
dy=diff(y);dy=[dy(1);dy(:);dy(end)];
dz=diff(z);dz=[dz(1);dz(:);dz(end)];
[DX,DY,DZ]=ndgrid(dx,dy,dz);
S=DX.*DY.*DZ.*sigma;
SIGMA=S(1:end-1,1:end-1,1:end-1)+S(2:end,1:end-1,1:end-1)+...
    S(1:end-1,2:end,1:end-1)+S(1:end-1,1:end-1,2:end)+...
    S(2:end,2:end,1:end-1)+S(2:end,1:end-1,2:end)+...
    S(1:end-1,2:end,2:end)+S(2:end,2:end,2:end);
Vp=(DX(1:end-1,2:end,2:end)+DX(2:end,2:end,2:end)).*...
    (DY(2:end,1:end-1,2:end)+DY(2:end,2:end,2:end)).*...
    (DZ(2:end,2:end,1:end-1)+DZ(2:end,2:end,2:end));
SIGMA=SIGMA./Vp;
switch type,
    case 1, % B,T & W
    case 2, % Z,M & M
        sxyz=sigma.*DX.*DY.*DZ;
    case 3, % W&M
    case 4,
    otherwise, % D&M or B,T&W
        sxyz=sigma.*DX.*DY.*DZ;
        zC=(sxyz(1:I,1:J,1:K)+sxyz(2:I+1,1:J,1:K)+sxyz(1:I,2:J+1,1:K)+sxyz(2:I+1,2:J+1,1:K))./DZ(1:I,1:J,1:K).^2*(-0.25);
        Cz=(sxyz(1:I,1:J,2:K+1)+sxyz(2:I+1,1:J,2:K+1)+sxyz(1:I,2:J+1,2:K+1)+sxyz(2:I+1,2:J+1,2:K+1))./DZ(1:I,1:J,2:K+1).^2*(-0.25);
        yC=(sxyz(1:I,1:J,1:K)+sxyz(2:I+1,1:J,1:K)+sxyz(1:I,1:J,2:K+1)+sxyz(2:I+1,1:J,2:K+1))./DY(1:I,1:J,1:K).^2*(-0.25);
        Cy=(sxyz(1:I,2:J+1,1:K)+sxyz(2:I+1,2:J+1,1:K)+sxyz(1:I,2:J+1,2:K+1)+sxyz(2:I+1,2:J+1,2:K+1))./DY(1:I,2:J+1,1:K).^2*(-0.25);
        xC=(sxyz(1:I,1:J,1:K)+sxyz(1:I,2:J+1,1:K)+sxyz(1:I,1:J,2:K+1)+sxyz(1:I,2:J+1,2:K+1))./DX(1:I,1:J,1:K).^2*(-0.25);
        Cx=(sxyz(2:I+1,1:J,1:K)+sxyz(2:I+1,2:J+1,1:K)+sxyz(2:I+1,1:J,2:K+1)+sxyz(2:I+1,2:J+1,2:K+1))./DX(2:I+1,1:J,1:K).^2*(-0.25);
end
Cs=-(Cx+Cy+Cz+xC+yC+zC);

xC(1,:,:)=0.0;Cx(I,:,:)=0.0;
yC(:,1,:)=0.0;Cy(:,J,:)=0.0;
Cz(:,:,K)=0.0;

if bound,
    % Mixed boundary conditions
    %% C_0 = C_0 - C_n*cos(n,r)*dx/r
    % Take Midpoint(surface) as source
    xm=(x(1)+x(I))/2;
    ym=(y(1)+y(J))/2;
    zm=z(1);
    % cos(n1,n2)=(n1.n2)/(abs(n1)*abs(n2));  n2=[1 0 0];
    % cos(r,e)/r = r.e / (r^2)
    RR=[reshape(repmat(x-xm,1,J*K),1,IJK);reshape(repmat(repmat(y-ym,1,K),I,1),1,IJK);reshape(repmat(z-zm,I*J,1),1,IJK)];
    rq=reshape(sum(RR.^2),I,J,K);
    % Left boundary
    i=1;xx=xm-x(i);
    Cs(1,:,:)=Cs(1,:,:)-dx(1)*Cx(1,:,:)*(xm-x(1))./rq(1,:,:);
    % Right boundary
    Cs(end,:,:)=Cs(end,:,:)-dx(end)*xC(end,:,:)*(x(end)-xm)./rq(end,:,:);
    % Front boundary
    Cs(:,1,:)=Cs(:,1,:)-dy(1)*Cy(:,1,:)*(ym-y(1))./rq(:,1,:);
    %Back boundary
    Cs(:,end,:)=Cs(:,end,:)-dy(end)*yC(:,end,:)*(y(end)-ym)./rq(:,end,:);
    % Bottom boundary
    k=K;zz=z(k)-zm;
    Cs(:,:,end)=Cs(:,:,end)-dz(end)*zC(:,:,end)*(z(end)-zm)./rq(:,:,end);
end    % end of boundary conditions

C=spdiags([Cz(:) Cy(:) Cx(:) Cs(:) xC(:) yC(:) zC(:)],[-I*J -I -1 0 1 I I*J],IJK,IJK);