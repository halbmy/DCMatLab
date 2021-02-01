function S=asens(x,y,z,N,M,S,rho0,R)

% ASENS - Sensitivity by anomaleous potential
% S = asens(x,y,z,N,M,S,rho_0,R)
% x/y/z - Cell coordinates
% N - vector of electrode numbers(a,b,m,n)
%               and konfiguration factors(k)
%               elec - electrode positions
% M - Model parameter
% S - (old) sensitivity matrix
% rho_0 - Primary Potential Resistivity
% R - Forward calculation

rhoBg=getbg(x,y,M(:,:,1),N.elec);
data=length(N.r);
ijk=prod(size(M));
if nargin<5, error('Not Enough Input arguments!');end
if nargin<6, S=zeros(data,ijk); end
if nargin<7, rho0=M(1,1,1); end
M1=M;
if nargin<8, % nicht logarithmisch
    R=ones(size(N.r));
    M1=ones(size(M));
end
data=length(N.a);
global PHI
[X,Y,Z]=ndgrid(x,y,z);
[DX,DY,DZ]=ndgrid(diff(x),diff(y),diff(z));
V=DX.*DY.*DZ*rho0/(32*pi)./M.^2; % 4 phis 4 rs und 2 pi
wb=waitbar(0,'Calculating anomalous sensitivities...');
aller=fix(data/25);  % show waitbar every aller times
mal=aller;
phi=zeros(length(x),length(y),length(z));
psi=phi;
dr=min(diff(z))/5;
% R1=Phi_p, R2=Psi_p (ohne 1/2pisigma)
warning('off','MATLAB:divideByZero');
for l = 1:data,
    phi(:)=PHI(:,N.a(l));
    RA=N.elec(N.a(l),:);
    X1=X-RA(1);Y1=Y-RA(2);Z1=Z-RA(3);
    R1=1./sqrt(X1.^2+Y1.^2+Z1.^2);
    R1(find(isinf(R1)))=1/dr;
    if N.b(l)>0,
        phi(:)=phi(:)-PHI(:,N.b(l));
        RB=N.elec(N.b(l),:);
        X1=X-RB(1);Y1=Y-RB(2);Z1=Z-RB(3);
        R1=R1-1./sqrt(X1.^2+Y1.^2+Z1.^2);
        R1(find(isinf(R1)))=-1/dr;
    end
    psi(:)=PHI(:,N.m(l));
    RM=N.elec(N.m(l),:);
    X2=X-RM(1);Y2=Y-RM(2);Z2=Z-RM(3);
    R2=1./sqrt(X2.^2+Y2.^2+Z2.^2);
    R2(find(isinf(R2)))=1/dr;
    if N.n(l)>0,
        psi(:)=psi(:)-PHI(:,N.n(l));
        RN=N.elec(N.n(l),:);
        X2=X-RN(1);Y2=Y-RN(2);Z2=Z-RN(3);
        R2=R2-1./sqrt(X2.^2+Y2.^2+Z2.^2);
        R2(find(isinf(R2)))=-1/dr;
    end
    % phix = Flächenmittel der x-Ableitung von phi
    phix=diff(phi(:,1:end-1,1:end-1)+phi(:,1:end-1,2:end)+phi(:,2:end,1:end-1)+phi(:,2:end,2:end),1,1)./DX;
    phiy=diff(phi(1:end-1,:,1:end-1)+phi(1:end-1,:,2:end)+phi(2:end,:,1:end-1)+phi(2:end,:,2:end),1,2)./DY;
    phiz=diff(phi(1:end-1,1:end-1,:)+phi(1:end-1,2:end,:)+phi(2:end,1:end-1,:)+phi(2:end,2:end,:),1,3)./DZ;
    psix=diff(psi(:,1:end-1,1:end-1)+psi(:,1:end-1,2:end)+psi(:,2:end,1:end-1)+psi(:,2:end,2:end),1,1)./DX;
    psiy=diff(psi(1:end-1,:,1:end-1)+psi(1:end-1,:,2:end)+psi(2:end,:,1:end-1)+psi(2:end,:,2:end),1,2)./DY;
    psiz=diff(psi(1:end-1,1:end-1,:)+psi(1:end-1,2:end,:)+psi(2:end,1:end-1,:)+psi(2:end,2:end,:),1,3)./DZ;
    % X1 = Summe der R1 über benachbarte y,z usw.
    X1=R1(:,1:end-1,1:end-1)+R1(:,1:end-1,2:end)+R1(:,2:end,1:end-1)+R1(:,2:end,2:end);
    Y1=R1(1:end-1,:,1:end-1)+R1(1:end-1,:,2:end)+R1(2:end,:,1:end-1)+R1(2:end,:,2:end);
    Z1=R1(1:end-1,1:end-1,:)+R1(1:end-1,2:end,:)+R1(2:end,1:end-1,:)+R1(2:end,2:end,:);
    X2=R2(:,1:end-1,1:end-1)+R2(:,1:end-1,2:end)+R2(:,2:end,1:end-1)+R2(:,2:end,2:end);
    Y2=R2(1:end-1,:,1:end-1)+R2(1:end-1,:,2:end)+R2(2:end,:,1:end-1)+R2(2:end,:,2:end);
    Z2=R2(1:end-1,1:end-1,:)+R2(1:end-1,2:end,:)+R2(2:end,1:end-1,:)+R2(2:end,2:end,:);
    sens1=diff(X1,1,1).*psix./DX.^2+diff(Y1,1,2).*psiy./DY.^2+diff(Z1,1,3).*psiz./DZ.^2;
    sens2=diff(X2,1,1).*phix./DX.^2+diff(Y2,1,2).*phiy./DY.^2+diff(Z2,1,3).*phiz./DZ.^2;
%     X1=diff(R1,1,1).*diff(psi,1,1);
%     Y1=diff(R1,1,2).*diff(psi,1,2);
%     Z1=diff(R1,1,3).*diff(psi,1,3);
%     X2=diff(R2,1,1).*diff(phi,1,1);
%     Y2=diff(R2,1,2).*diff(phi,1,2);
%     Z2=diff(R2,1,3).*diff(phi,1,3);
%     sens1=(X1(:,1:end-1,1:end-1)+X1(:,2:end,1:end-1)+X1(:,1:end-1,2:end)+X1(:,2:end,2:end))./DX.^2+...
%     (Y1(1:end-1,:,1:end-1)+Y1(2:end,:,1:end-1)+Y1(1:end-1,:,2:end)+Y1(2:end,:,2:end))./DY.^2+...
%     (Z1(1:end-1,1:end-1,:)+Z1(2:end,1:end-1,:)+Z1(1:end-1,2:end,:)+Z1(2:end,2:end,:))./DZ.^2;
%     sens2=(X2(:,1:end-1,1:end-1)+X2(:,2:end,1:end-1)+X2(:,1:end-1,2:end)+X2(:,2:end,2:end))./DX.^2+...
%     (Y2(1:end-1,:,1:end-1)+Y2(2:end,:,1:end-1)+Y2(1:end-1,:,2:end)+Y2(2:end,:,2:end))./DY.^2+...
%     (Z2(1:end-1,1:end-1,:)+Z2(2:end,1:end-1,:)+Z2(1:end-1,2:end,:)+Z2(2:end,2:end,:))./DZ.^2;
    sens3=(phix.*psix+phiy.*psiy+phiz.*psiz)*2*pi/rho0;
    sens=4*(sens1+sens2)+sens3;   
    sens=sens1+sens2+sens3;   
    S(l,:)=(S(l,:)-(sens(:).*V(:)*N.k(l))').*M1(:)'/R(l);
    mal=mal-1;
    if mal<1,
        waitbar(l/data,wb);
        mal=aller;
    end
end
warning('off','MATLAB:divideByZero');
close(wb);
