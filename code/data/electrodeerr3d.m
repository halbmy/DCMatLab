function [maxerr,E,A]=electrodeerr3d(N,dx,dy)

% ELECTRODEERR3D - Error caused by electrode displacement (3d)
% [A,maxerr] = electrodeerr3d(N,dx)
% N  - data structure
% dx - electrode displacement (default 0.01m)
% A - data covariance matrix
% maxerr - maximum err for each datum

nel=size(N.elec,1);
%E=zeros(length(N.a),nel);
E=spalloc(length(N.a),nel*2,length(N.a)*8);
E=spalloc(length(N.a),nel*2,length(N.a)*8);
sumerr=zeros(length(N.a),1);
if nargin<2, dx=0.01; end
if nargin<3, dy=dx; end
for e=1:nel,
    elx=N.elec(e,1);
    N.elec(e,1)=elx+dx;
    newk=getkonf(N);
    E(:,e)=(newk-N.k)./N.k;
    N.elec(e,1)=elx;
    ely=N.elec(e,2);
    N.elec(e,2)=ely+dx;
    newk=getkonf(N);
    E(:,nel+e)=(newk-N.k)./N.k;
    N.elec(e,2)=ely;
end
maxerr=sum(abs(E),2);
E=E';
if nargout>2,
    A=E'*E;
end
