function p=itantrans(m,l,u)

% ITANTRANS - Inverse tangens transformation with lower and upper bound
% p = itantrans(m,lower_bound,upperbound)
% both lower and upper bound are required
% See also TANTRANS

if nargin<3, error('Specify model, lower, and upper bound!'); end
p=atan(m)/pi*(u-l)+(l+u)/2;