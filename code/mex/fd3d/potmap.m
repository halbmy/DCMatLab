function C=potmap(elec,x,y,z)

% POTMAP - Matrix to extract potential data
% C = potmap(elec,x,y,z)

I=length(x);J=length(y);K=length(z);
IJK=I*J*K;
m=size(elec,1);
A=zeros(I,J,K);
C=spalloc(m,IJK,8*m);
for n=1:m,
  xx=elec(n,1);
  yy=elec(n,2);
  zz=elec(n,3);
  i=max(find(x<=xx));
  j=max(find(y<=yy));
  k=max(find(z<=zz));
  dv=(x(i+1)-x(i))*(y(j+1)-y(j))*(z(k+1)-z(k));
  dxi=xx-x(i);dx1=x(i+1)-xx;
  dyi=yy-y(j);dy1=y(j+1)-yy;
  dzi=zz-z(k);dz1=z(k+1)-zz;
  A(:)=0;
  A(i,j,k)=dx1*dy1*dz1/dv;
  A(i+1,j,k)=dxi*dy1*dz1/dv;
  A(i,j+1,k)=dx1*dyi*dz1/dv;
  A(i+1,j+1,k)=dxi*dyi*dz1/dv;
  A(i,j,k+1)=dx1*dy1*dzi/dv;
  A(i+1,j,k+1)=dxi*dy1*dzi/dv;
  A(i,j+1,k+1)=dx1*dyi*dzi/dv;
  A(i+1,j+1,k+1)=dxi*dyi*dzi/dv;
  C(n,:)=A(:)';
end

%if nargout==0,
    %spy(C(:,1:max(find(sum(C)))));
%end
