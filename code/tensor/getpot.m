function u=getjpot(A,P,B)

% GETJVEC - Get current density vector at point P
% j = getjvec(A,P)
% j = getjvec(A,P,B)
% [jx,jy] = getjvec(A,P[,B]);
% P - (x,y)-row-vector of point (can also be matrix of many)
% A/B - vector of positive/negative current injection (can also be matrix)

dx=A(:,1)-P(:,1);
dy=A(:,2)-P(:,2);
r=sqrt(dx.^2+dy.^2);
u=1./r;
if nargin>2, % B given
    dx=B(:,1)-P(:,1);
    dy=B(:,2)-P(:,2);
    r=sqrt(dx.^2+dy.^2);
    u=u-1./r;
end
u=u/2/pi;