function P = pmatrix2d(x,z,S,mincov)

% PMATRIX - Form p-matrix from sensitivity by minimum coverage
% P = pmatrix2d(x,z)
% or P = pmatrix2d(x,z,S,mincov)

if nargin<2, error('2 input arguments required!'); end
if nargin<4, mincov=0.5; end
lx=length(x)-1;
lz=length(z)-1;
dx=mean(diff(x));
dz=diff(z);
P=speye(lx*lz);
for k=lz:-1:2,
  start=lx*(k-1);
  nx=round(dz(k)/dx);
  if nx>1,
    kx=fix(lx/nx);
    nn=nx+lx-nx*kx;
    for kk=kx:-1:1,
      first=start+(kk-1)*nx+1;
      P(:,first)=sum(P(:,first:first+nn-1),2);
      P(:,first+1:first+nn-1)=[];
      nn=nx;
    end
  end
end
if nargin>2, %mincov given
    COV=sum(abs(S))*P./sum(P);
    P(:,COV<mincov)=[];
%     COV=sum(abs(S));
%     P(:,find((COV<mincov)*P))=[];
end
% if size(P,2)<size(P,1),
%     message(sprintf('reduced parameter from %d to %d',size(P,1),size(P,2)));
% end
