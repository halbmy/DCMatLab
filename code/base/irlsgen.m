function weight = irlsGEN(a,p,locut,hicut)

% IRLS - iteratively reweighted least squares function
% weight = irls(vector[,locut,hicut])
% locut/hicut restrict the values to lie above/below bounds
%            |a_i|^p / || a ||_p^p            sum(a^2) * |a|^p
% weight_i = ------------------- => weight =  -----------------
%            a_i^2 / || a ||_2^2              a^2 * sum(|a|^p)

if nargin<1, error('Specify vector'); end
if nargin<2, p=1; end
if nargin<3, locut=0; end
if nargin<4, hicut=1; end

absp=abs(a).^p;
weight=ones(size(a));
fi = find(absp);
weight(fi) = absp(fi)./(a(fi).^2) * ( a(:)'*a(:) ) / sum(absp);
weight( weight < locut ) = locut;
if hicut>0, weight( weight > hicut ) = hicut; end
