function x=fillgap(d,a,b)

% FILLGAP - Fill gap of length d with final spacing a
%           and (optional) initial starting b
% x = fillgap(d,a[b])

if nargin<3, b=a/2; end
n=floor(d/a)+1;
if n==1,
    x=0;
else
    ddx=(d-n*b)/sum(1:n-1);
    x=[0 cumsum([0:n-2]*ddx+b)];
end