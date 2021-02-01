function dm=tantransdiff(p,l,u)

% TANTRANS - Tangens transformation with lower and upper bound
% m = tantrans(p,lower_bound,upperbound)
% both lower and upper bound are required
% See also ITANTRANS

if nargin<3, error('Specify model, lower, and upper bound!'); end
dm = pi / ( u - l ) * ( 1 + tantrans(p,l,u).^2 );