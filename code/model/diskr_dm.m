function C = diskr_dm(X,Y,Z,SIGMA,rb)
% DISKR_DM - Discretization after Dey&Morrison(1979)
% C = diskr_dm(x,y,z,sigma)
% x,y,z - grid nodes
% sigma - conductivities (lx+1).(ly+1).(lz+1)
% with sigma(:,:,1)=0;
if nargin<5, rb=0; end
if nargin<3, 
    SIGMA=ones(length(X)+1,length(Y)+1,length(Z)+1);
    SIGMA(:,:,1)=0;
end
I=length(X);
J=length(Y);
K=length(Z);
IJK=I*J*K;
C=spalloc(IJK,IJK,7*IJK);
Cs=zeros(IJK,1);Cx=reshape(Cs,I,J,K);
Cy=Cx;Cz=Cx;xC=Cx;yC=Cx;zC=Cx;

dX=X(2:I)-X(1:I-1);
dY=Y(2:J)-Y(1:J-1);
dZ=Z(2:K)-Z(1:K-1);
dX=[dX(1) dX dX(I-1)];
dY=[dY(1) dY dY(J-1)];
dZ=[dZ(1) dZ dZ(K-1)];

dX=reshape(dX,I+1,1,1);dY=reshape(dY,1,J+1,1);dZ=reshape(dZ,1,1,K+1);

DX=reshape(repmat(dX,1,J*K),I+1,J,K);
DY=reshape(repmat(repmat(dY,1,K),I,1),I,J+1,K);
DZ=reshape(repmat(dZ,I*J,1),I,J,K+1);
SIGMA(:,:,1)=0; % for CQ !

%% TOP (-z)
   zC=SIGMA(1:I,2:J+1,1:K).*DX(1:I,:,:).*DY(:,2:J+1,:);
zC=zC+SIGMA(2:I+1,2:J+1,1:K).*DX(2:I+1,:,:).*DY(:,2:J+1,:);
zC=zC+SIGMA(1:I,1:J,1:K).*DX(1:I,:,:).*DY(:,1:J,:);
zC=zC+SIGMA(2:I+1,1:J,1:K).*DX(2:I+1,:,:).*DY(:,1:J,:);
zC=-0.25*zC./DZ(:,:,1:K);
%% BOTTOM (+z)
   Cz=SIGMA(1:I,2:J+1,2:K+1).*DX(1:I,:,:).*DY(:,2:J+1,:);
Cz=Cz+SIGMA(2:I+1,2:J+1,2:K+1).*DX(2:I+1,:,:).*DY(:,2:J+1,:);
Cz=Cz+SIGMA(1:I,1:J,2:K+1).*DX(1:I,:,:).*DY(:,1:J,:);
Cz=Cz+SIGMA(2:I+1,1:J,2:K+1).*DX(2:I+1,:,:).*DY(:,1:J,:);
Cz=-0.25*Cz./DZ(:,:,2:K+1);
%% LEFT (-x)
   xC=SIGMA(1:I,1:J,1:K).*DY(:,1:J,:).*DZ(:,:,1:K);
xC=xC+SIGMA(1:I,2:J+1,1:K).*DY(:,2:J+1,:).*DZ(:,:,1:K);
xC=xC+SIGMA(1:I,1:J,2:K+1).*DY(:,1:J,:).*DZ(:,:,2:K+1);
xC=xC+SIGMA(1:I,2:J+1,2:K+1).*DY(:,2:J+1,:).*DZ(:,:,2:K+1);
xC=-0.25*xC./DX(1:I,:,:);
%% RIGHT (+x)
   Cx=SIGMA(2:I+1,1:J,1:K).*DY(:,1:J,:).*DZ(:,:,1:K);
Cx=Cx+SIGMA(2:I+1,2:J+1,1:K).*DY(:,2:J+1,:).*DZ(:,:,1:K);
Cx=Cx+SIGMA(2:I+1,1:J,2:K+1).*DY(:,1:J,:).*DZ(:,:,2:K+1);
Cx=Cx+SIGMA(2:I+1,2:J+1,2:K+1).*DY(:,2:J+1,:).*DZ(:,:,2:K+1);
Cx=-0.25*Cx./DX(2:I+1,:,:);
%% FRONT (-y)
   yC=SIGMA(1:I,1:J,1:K).*DX(1:I,:,:).*DZ(:,:,1:K);
yC=yC+SIGMA(2:I+1,1:J,1:K).*DX(2:I+1,:,:).*DZ(:,:,1:K);
yC=yC+SIGMA(1:I,1:J,2:K+1).*DX(1:I,:,:).*DZ(:,:,2:K+1);
yC=yC+SIGMA(2:I+1,1:J,2:K+1).*DX(2:I+1,:,:).*DZ(:,:,2:K+1);
yC=-0.25*yC./DY(:,1:J,:);
%% BACK (+y)
   Cy=SIGMA(1:I,2:J+1,1:K).*DX(1:I,:,:).*DZ(:,:,1:K);
Cy=Cy+SIGMA(2:I+1,2:J+1,1:K).*DX(2:I+1,:,:).*DZ(:,:,1:K);
Cy=Cy+SIGMA(1:I,2:J+1,2:K+1).*DX(1:I,:,:).*DZ(:,:,2:K+1);
Cy=Cy+SIGMA(2:I+1,2:J+1,2:K+1).*DX(2:I+1,:,:).*DZ(:,:,2:K+1);
Cy=-0.25*Cy./DY(:,2:J+1,:);

Cs=-(Cx+Cy+Cz+xC+yC+zC);

xC(1,:,:)=0.0;Cx(I,:,:)=0.0;
yC(:,1,:)=0.0;Cy(:,J,:)=0.0;
Cz(:,:,K)=0.0;

if rb,
% Mixed boundary conditions
%% C_0 = C_0 - C_n*cos(n,r)*dx/r
% Take Midpoint(surface) as source
xm=median(X);%(X(1)+X(I))/2;
ym=median(Y);%(Y(1)+Y(J))/2;
zm=Z(1);
% cos(n1,n2)=(n1.n2)/(abs(n1)*abs(n2));  n2=[1 0 0];
% cos(r,e)/r = r.e / (r^2)
RR=[reshape(repmat(X-xm,1,J*K),1,IJK);reshape(repmat(repmat(Y-ym,1,K),I,1),1,IJK);reshape(repmat(Z-zm,I*J,1),1,IJK)];
rq=reshape(sum(RR.^2),I,J,K);
% Left boundary
i=1;xx=xm-X(i);
Cs(i,:,:)=Cs(i,:,:)-dX(i)*Cx(i,:,:)*xx./rq(i,:,:);
% Right boundary
i=I;xx=X(i)-xm;
Cs(i,:,:)=Cs(i,:,:)-dX(i)*xC(i,:,:)*xx./rq(i,:,:);
% Front boundary
j=1;yy=ym-Y(j);
Cs(:,j,:)=Cs(:,j,:)-dY(j)*Cy(:,j,:)*yy./rq(:,j,:);
%Back boundary
j=J;yy=Y(j)-ym;
Cs(:,j,:)=Cs(:,j,:)-dY(j)*yC(:,j,:)*yy./rq(:,j,:);
% Bottom boundary
k=K;zz=Z(K)-zm;
Cs(:,:,k)=Cs(:,:,k)-dZ(k)*zC(:,:,k)*zz./rq(:,:,k);
% end of boundary conditions
end
C=spdiags([Cz(:) Cy(:) Cx(:) Cs(:) xC(:) yC(:) zC(:)],[-I*J -I -1 0 1 I I*J],IJK,IJK);