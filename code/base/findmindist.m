function [ind,midi,di]=findmindist(pl,pp)

% FINDMINDIST - Find minimum distance from point to list of points
% [ind,mindist,distances] = findmindist(pointlist,point)
% size(pointlist) = [N,D]
% size(point) = [1,D]

if nargin<2, pp=zeros(1,size(pl,2)); end
di=sqrt(sum((pl-repmat(pp,size(pl,1),1)).^2,2));
[midi,ind]=min(di);