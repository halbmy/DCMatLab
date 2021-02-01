function [L,Cx,Cy,Cz] = smooth3d1st(x,y,z,nor,az)

% SMOOTHMAT2D - Computes smoothness matrix
% L = smooth3dst(x,y,z[,nor,az])
% [L,Cx,Cy,Cz] = smooth3d1st(...
% L = Cx'*Cx + Cy'*Cy + Cz'*Cz
% x/y/z .. model grid vectors
% nor   .. normalize vector (1 2 3...)
% az    .. weight for z-direction derivatives

argmin=3;
if nargin<argmin, error('Too less input arguments!'); end
if nargin<4, nor=0; end
if nargin<5, az=1; end
if nor,
  x=1:length(x);
  y=1:length(y);
  z=1:length(z);
end
nx=length(x)-1;ny=length(y)-1;nz=length(z)-1;
x=(x(1:end-1)+x(2:end))/2; % Midpoint
y=(y(1:end-1)+y(2:end))/2;
z=(z(1:end-1)+z(2:end))/2;
nn=nx*ny*nz;
one=repmat([1./diff(x(:));0],ny*nz,1);
Cx=spdiags([-one one],[0 1],nn-1,nn);
Cx(nx:nx:end,:)=[];
two=repmat([repmat(1./diff(y(:)),nx,1);zeros(nx,1)],nz,1);
Cy=spdiags([-two(:) two(:)],[0 nx],nn-nx,nn);
su=sum(abs(Cy),2);
Cy(find(su==0),:)=[];
%Cy(nx*ny:nx*ny:end,:)=[];
three=repmat(1./diff(z(:)),nx*ny,1);
Cz=spdiags([-three(:) three(:)],[0 nx*ny],nn-nx*ny,nn);
L=Cx'*Cx+Cy'*Cy+az*Cz'*Cz;
%if nargout==0, % only for testing
    %clf;
    %spy(L);
    %L=1;
%end