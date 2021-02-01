function S=spcalcsens3dt(x,y,z,N,sp,spratio)

% CALCSENS3D - Calculate sensitivity matrix by gauss-legendre integration
% S = calcsens3d(x,y,z,N)
% x,y,z - limits of model cells
% N - data structur with electrode positions(elec),
%     numbers(a,b,m,n) and k-factors(k)
% bases on fortran routine sens3dfull.f (Thomas Günther)
% it uses - dependent on electrode distance - from 2^3 to 10^3 points
% to integrate the sensitivity function over hexahedral domains

if nargin<6, spratio=10; end
if nargin<5, sp=5e-4; end
lx=length(x)-1;ly=length(y)-1;lz=length(z)-1;
ijk=lx*ly*lz;
data=length(N.a);
S=spalloc(ijk,data,round(data*ijk*spratio/100)); %5));
wb=waitbar(0,'Integrating Sensitivity matrix (sparse) ...');
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
    sens=sens*N.k(l);
%     sens(abs(sens)<sp*max(sens))=0;
%    sens(abs(sens)<sp*sum(sens))=0;
%    S(l,:)=sparse(sens);
%     mo=sp*sum(abs(sens));
%     fi=find(abs(sens)>mo);
     fi=find(abs(sens)>sp);
     S(fi,l)=sens(fi);
    mal=mal-1;
    if mal<1,
        waitbar(l/data,wb);
        if nargin>6, fprintf('%5d data: %d/%d=%.1f%% nnz\n',l,nnz(S),l*ijk,nnz(S)/l/ijk*100); end
        mal=aller;
    end
end
close(wb);
