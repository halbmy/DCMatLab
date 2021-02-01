function b=logdrop(a,mina)
if nargin<2, mina=max(1e-6,min(abs(a))); end
lma=log10(mina);
b=(max(log10(abs(a)),lma)-lma).*sign(a);