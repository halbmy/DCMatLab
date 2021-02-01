datafile='ro/ro_dd.dat';
% datafile='rotest/rothschoen.dat';
% datafile='soos/innen/mof-innen.dat';
N=read3dfile(datafile);
N.err=estimateerror(N,3);
[x,y,z,M]=modelfromdata3d(N);
rho=M(1,1,1);
S=calcsens3d(x,y,z,N);
C=smooth3d1st(x,y,z,1);
D=spdiags(1./log(1+N.err),0,length(N.err),length(N.err));
dR=log(N.r)-log(rho);
dM=cglscdp(S,dR,30,C,D);
tol=7e-4;
S1=spcalcsens3d(x,y,z,N,tol);
fprintf('Elements: full=%d sparse=%d pack=%.1f%%\n',...
  prod(size(S)),nnz(S1),nnz(S1)/prod(size(S))*100);
dM1=cglscdp(S1,dR,30,C,D);
plot(dM,dM1,'.');
S2=spcalcsens3dnew(x,y,z,N,nnz(S1)/prod(size(S1)));
fprintf('Elements: full=%d sparse=%d pack=%.1f%%\n',...
  prod(size(S)),nnz(S2),nnz(S2)/prod(size(S))*100);
dM2=cglscdp(S2,dR,30,C,D);
plot(dM,dM2,'.');
[norm(dM1./dM-1) norm(dM2./dM-1)]
[norm(dM1-dM) norm(dM2-dM)]