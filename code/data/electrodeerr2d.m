function [maxerr,A,E]=electrodeerr2d(N,dx)

% ELECTRODEERR2D - Error caused by electrode displacement (2d)
% [maxerr,A,E] = electrodeerr2d(N,dx)
% N  - data structure
% dx - electrode displacement (default 0.01m)
% maxerr - maximum err for each datum
% A - data covariance matrix
% E - Electrode-Data dependency matrix

nel=size(N.elec,1);
%E=zeros(length(N.a),nel);
E=spalloc(length(N.a),nel,length(N.a)*4);
% sumerr=zeros(length(N.a),1);
if nargin<2, dx=0.01; end
oldk=getkonf(N);
for e=1:nel,
    elx=N.elec(e,1);
    N.elec(e,1)=elx+dx;
    newk=getkonf(N);
    E(:,e)=(newk-oldk)./oldk;
    N.elec(e,1)=elx;
end
maxerr=sum(abs(E),2);
if nargout>1, A=E*E'; end
if nargout>2, E=E'; end
