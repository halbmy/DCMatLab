function S = sens2dbrute(M,x,z,Lay,N,FOR)

% SENS2DBRUTE - 2D Sensitivity by brute force (sorry computer!)
% S = sens2dbrute(M,x,z,Lay,N,FOR)

optfor=FOR;optfor.silent=1;
modcells=prod(size(M));
S=zeros(length(N.r),modcells);
R0=dcfwd2d(x,z,M,Lay,N,optfor);
fac=1.05;
wb=waitbar(0,'Brute force Sensitivity...');
aller=fix(modcells/25);
mal=aller;
for m=1:modcells,
    mm=M(m); 
    M(m)=mm*fac;
    R=dcfwd2d(x,z,M,Lay,N,optfor);
    S(:,m)=(log(R(:))-log(R0(:)))/log(fac);
    mal=mal-1;
    if mal==0,
        waitbar(m/modcells);
        mal=aller;
    end
    M(m)=mm;
end
close(wb);