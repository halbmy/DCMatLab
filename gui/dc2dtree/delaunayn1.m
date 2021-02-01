function t = delaunayn(x)
%DELAUNAYN  N-D Delaunay tessellation.
%   T = DELAUNAYN(X) returns a set of simplices such that no data
%   points of X are contained in any circumspheres of the simplices.
%   The set of simplices forms the Delaunay tessellation. X is an
%   m-by-n array representing m points in n-D space. T is a
%   numt-by-(n+1) array where each row is the indices into X of the
%   vertices of the corresponding simplex. When the simplices cannot
%   be computed (such as when X is degenerate, or X is empty), an
%   empty matrix is returned.
%
%   DELAUNAYN is based on Qhull.  
%
%   See also QHULL, VORONOIN, CONVHULLN, DELAUNAY, DELAUNAY3.

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.1.1.1 $ $Date: 2005/06/01 10:29:16 $

if isempty(x), t = []; return; end

[x,idx,jdx] = unique(x,'rows');
[m,n] = size(x);

if m < n+1,
  error('Not enough unique points to do tessellation.');
end
if any(isinf(x(:)) | isnan(x(:)))
  error('Data containing Inf or NaN cannot be tessellated.');
end
if m == n+1
  t = 1:n+1;
  return;
end

% qhull needs the sum of squares at the end of the points
x = [x sum(x.*x,2)];
t = qhullmx(x', 'd ', 'QJ ', 'Pp');

x(:,end) = [];

% try to get rid of zero volume simplices. They are generated
% because of the fuzzy jiggling.

[m, n] = size(t);
v = true(m,1);

seps = eps^(4/5)*max(abs(x(:)));
for i=1:m
  if abs(det(x(t(i,1:n-1),:)-x(t(i,n)*ones(n-1,1),:))) < seps
    v(i) = 0;
  end
end

t = t(v,:);
t = idx(t);
