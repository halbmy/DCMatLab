function S = brutesens2d(x,z,M,Lay,N,FOR)

% BRUTESENS2D - Brute force 2d sensitivity calculation
% S = brutesens2d(x,z,M,Lay,N,FOR)
% x,z .. grid nodes in x/z direction
% M   .. model resistivity  matrix (length(x)-1)*(length(z)-1)
% Lay .. background resistivities
% N   .. data structure with elec=electrodes and a/b/m/n numbers
% FOR .. forward calculation options
if nargin<6,
    FOR=struct('fillup',0,'rand',3,'prolong',5,'zref',1);
end
FOR.silent=1;
fak=1.01; % 1% perturbation
R=dcfwd2d(x,z,M,Lay,N,FOR);
nd=length(N.r);
nm=numel(M);
S=zeros(nd,nm);
aller=fix(nm/25);
if aller==0, aller=100; end
mal=aller;
wb=waitbar(0,'Sensitivity calculation');
for l=1:nm,
    mm=M(l);
    M(l)=mm*fak;
    R1=dcfwd2d(x,z,M,Lay,N,FOR);
    S(:,l)=(log(R1)-log(R))/log(fak);
    M(l)=mm;
    mal=mal-1;
    if mal==0,
        waitbar(l/nm,wb);
        mal=aller;
    end
end
close(wb);