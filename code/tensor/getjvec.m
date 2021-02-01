function [jx,jy]=getjvec(A,P,B)

% GETJVEC - Get current density vector at point P
% j = getjvec(A,P)
% j = getjvec(A,P,B)
% [jx,jy] = getjvec(A,P[,B]);
% P - (x,y)-row-vector of point (can also be matrix of many)
% A/B - vector of positive/negative current injection (can also be matrix)

dx=P(:,1)-A(:,1);
dy=P(:,2)-A(:,2);
r=sqrt(dx.^2+dy.^2);
jx=dx./r.^3;
jy=dy./r.^3;
if nargin>2, % B given
    dx=P(:,1)-B(:,1);
    dy=P(:,2)-B(:,2);
    r=sqrt(dx.^2+dy.^2);
    jx=jx-dx./r.^3;
    jy=jy-dy./r.^3;    
end
jx=jx/2/pi;
jy=jy/2/pi;
if nargout<2,
    jx=[jx jy];
end