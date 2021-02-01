function [ii,jj,kk]=antiindex(mm,I,J)

% ANTIINDEX - determines three indices out of a single
% [i,j,k] = antiindex(n,I,J)
 
if nargin<3,
  global X Y
  I=length(X);
  J=length(Y);
end
ii=mod(mm,I);
nn=floor((mm-ii)/I);
jj=mod(nn,J)+1;
kk=floor(nn/J)+1;
return;
