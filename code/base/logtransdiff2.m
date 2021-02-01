function m=logtransdiff(p,l,u)

% LOGTRANS - Logarithmic transformation with lower/upper bound
% m = logtrans(p,[lower_bound[,upper_bound]])
% if upperbound is not specified or zero only lowerbound is used
% if lowerbound is not specified it is set to zero
% see also ILOGTRANS

if nargin<2, l=0; end
if nargin<3, u=l; end
m=1./(p-l);
if u>l, m=m+1./(u-p); end
