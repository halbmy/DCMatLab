function S=calcsens3d(x,y,z,N,sp)

% CALCSENS3D - Calculate sensitivity matrix by gauss-legendre integration
% S = calcsens3d(x,y,z,N)
% x,y,z - limits of model cells
% N - data structur with electrode positions(elec),
%     numbers(a,b,m,n) and k-factors(k)
% bases on fortran routine sens3dfull.f (Thomas Günther)
% it uses - dependent on electrode distance - from 2^3 to 10^3 points

if nargin<5, sp=1e-3; end
lx=length(x)-1;ly=length(y)-1;lz=length(z)-1;
ijk=lx*ly*lz;
data=length(N.a);
S=spalloc(data,ijk,round(data*ijk/15));
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
    sens=sens'*N.k(l);
    tol=1e-5*max(sens);ngt=length(sens)*sp;
    asens=abs(sens);
    while length(find(asens>tol))>ngt, tol=tol*1.3; end
    sens(asens<tol)=0;
    S(l,:)=sparse(sens);
    mal=mal-1;
    if mal<1,
        waitbar(l/data,wb);
        mal=aller;
    end
end
close(wb);
