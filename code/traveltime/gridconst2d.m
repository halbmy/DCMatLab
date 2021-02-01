function C = gridconst2d(M)

% GRIDCONST2D - Grid constraint matrix
% C = gridconst(M)

% M=ones(3,2);
alphaz=0.3;
nx=size(M,1);nz=size(M,2);
hor=(nx-1)*nz;ver=nx*(nz-1);
% C=spalloc(hor+ver,nx*nz,(hor+ver)*2);
l=0;
mm=[];nn=[];vv=[];
ox=ones(nx-1,1);oz=ones(nz-1,1)*alphaz;
ii=(1:nx-1)';kk=(1:nz-1)';
for k=1:nz,
    ik=(k-1)*nx+ii;
    mm=[mm;l+ii];
    nn=[nn;ik];
    vv=[vv;ox];
    mm=[mm;l+ii];
    nn=[nn;ik+1];
    vv=[vv;-ox];
    l=l+nx-1;
end
for i=1:nx,
    ik=(kk-1)*nx+i;
    mm=[mm;l+kk];
    nn=[nn;ik];
    vv=[vv;oz];
    mm=[mm;l+kk];
    nn=[nn;ik+nx];
    vv=[vv;-oz];
    l=l+nz-1;
end
C=sparse(mm,nn,vv);