function erg=getdoi(Rhoa,S,rho,fak,INV,err);

% DOI - Depth of investigation(after Oldenburgh&Li)
% Dependance from strting model by two inverse steps
% erg = getdoi(Rhoa,S,rho1[,factor[,INV]]);
% Rhoa - Measured Resistivities
% S - Sensitivity
% rho1 - starting rho
% factor - for second rho(1.1)
% INV - Inversion settings
% err - data errors(std)

if nargin<5,
    INV=struct('redu',0,'mitschicht',0,'method',2,'auto',1,'start',1);  
end
if nargin<6, err=0.01; end
if nargin<4, fak=1.1; end
inv=INV;
inv.method=0;
inv.lam=100;
inv.const=1;
MOD.rand=3;
M=ones(size(S,2))*rho;
R=ones(size(Rhoa))*rho;
dR=log(Rhoa(:))-log(R(:));
dM1=invers(S,dR,inv,err);
rho=rho*fak;
M=ones(size(S,2))*rho;
R=ones(size(Rhoa))*rho;
dR=log(Rhoa(:))-log(R(:));
dM2=invers(S,dR,inv,err);
erg=abs(log(fak)+(dM2-dM1))/log(fak);
