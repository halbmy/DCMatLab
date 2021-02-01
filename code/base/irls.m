function weight = irls(a,locut,hicut)

% IRLS - iteratively reweighted least squares function
% weight = irls(vector[,locut,hicut])
% locut/hicut restrict the values to lie above/below bounds
%            |a_i| / || a ||_1^1                 sum(a^2)
% weight_i = ------------------- => weight =  --------------
%            a_i^2 / || a ||_2^2              |a| * sum(|a|)

if nargin<1, error('Specify vector'); end
if nargin<2, locut=0; end
if nargin<3, hicut=1; end

absa=abs(a);
weight=ones(size(a));
fi = find(isfinite(a)&(absa>0));
if isempty(fi), return; end
weight(fi) = ( sum(a(fi).^2) ) / sum(absa(fi)) ./ absa(fi) ;
weight( weight < locut ) = locut;
if hicut>0, weight( weight > hicut ) = hicut; end
