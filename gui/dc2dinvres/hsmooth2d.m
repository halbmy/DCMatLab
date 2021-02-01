function C=hsmooth2d(M,n)

% HSMOOTH2D - 2D HORIZONTAL SMOOTHNESS CONSTRAINT MATRIX
% C = hsmooth2d(M,n)

if nargin<2, n=1; end
m=prod(size(M));
spalte=ones(m,1);
C=spdiags([-spalte spalte],[0 1],m-1,m);
nx=size(M,1);
C(nx:nx:end,:)=[];