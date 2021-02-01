function S=spcalcsens2dt(x,z,N,sp)

% SPCALCSENS3D - Calculate (sparse) sensitivity matrix by gauss-legendre integration
% S = spcalcsens2dt(x,y,z,N)
% x,z - limits of model cells
% N - data structure with electrode positions(elec),
%     numbers(a,b,m,n) and k-factors(k)
% bases on fortran routine sens2dxz.f (Thomas Günther)
% it uses - dependent on electrode distance - from 2^3 to 10^3 points
% to integrate the sensitivity function over hexahedral domains

if nargin<5, sp=5e-4; end
lx=length(x)-1;lz=length(z)-1;
ik=lx*lz;
data=length(N.a);
S=spalloc(ik,data,round(data*ik/10)); %5));
wb=waitbar(0,'Integrating Sensitivity matrix (sparse) ...');
aller=fix(data/25);  % show waitbar every aller times
mal=aller;
for l = 1:data,
    sens=sens2dxz(x,z,N.elec(N.a(l),:),N.elec(N.m(l),:));
    if N.n(l)>0,
        sens=sens-sens2dxz(x,z,N.elec(N.a(l),:),N.elec(N.n(l),:));
        if N.b(l)>0,
            sens=sens+sens2dxz(x,z,N.elec(N.b(l),:),N.elec(N.n(l),:));
        end
    end
    if N.b(l)>0,
        sens=sens-sens2dxz(x,z,N.elec(N.b(l),:),N.elec(N.m(l),:));
    end
    sens=sens*N.k(l);
    fi=find(abs(sens)>sp);
    S(fi,l)=sens(fi);
    mal=mal-1;
    if mal<1,
        waitbar(l/data,wb);
        mal=aller;
    end
end
close(wb);
