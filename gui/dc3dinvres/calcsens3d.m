function S=calcsens3d(x,y,z,N)

% CALCSENS3D - Calculate sensitivity matrix by gauss-legendre integration
% S = calcsens3d(x,y,z,N)
% x,y,z - limits of model cells
% N - data structur with electrode positions(elec),
%     numbers(a,b,m,n) and k-factors(k)
% bases on fortran routine sens3dfull.f (Thomas Günther)
% it uses - dependent on electrode distance - from 2^3 to 10^3 points

lx=length(x)-1;ly=length(y)-1;lz=length(z)-1;
ijk=lx*ly*lz;
data=length(N.a);
S=zeros(data,ijk);
wb=waitbar(0,'Integrating Sensitivity matrix...');
aller=fix(data/25);  % show waitbar every aller times
mal=aller;
for l = 1:data,
    sens=sens3dfull(x,y,z,N.elec(N.a(l),:),N.elec(N.m(l),:));
    if N.n(l)>0,
        sens=sens-sens3dfull(x,y,z,N.elec(N.a(l),:),N.elec(N.n(l),:));
        if N.b(l)>0,
            sens=sens+sens3dfull(x,y,z,N.elec(N.b(l),:),N.elec(N.n(l),:));
        end
    end
    if N.b(l)>0,
        sens=sens-sens3dfull(x,y,z,N.elec(N.b(l),:),N.elec(N.m(l),:));
    end
    S(l,:)=sens'*N.k(l);
    mal=mal-1;
    if mal<1,
        waitbar(l/data,wb);
        mal=aller;
    end
end
close(wb);