function P = pmincov(S,mincov)

% PMINCOV - For p-matrix from sensitivity by minimum coverage
% P = pmincov(S,mincov)
% mincov - minimum coverage (default 0.4)
if nargin<2, mincov=0.4; end
P=speye(size(S,2));
COV=ones(size(S,2),1);
for i=1:length(Cov), COV(i)=sum(abs(S(:,i))); end
P(:,find(COV<mincov))=[];
if size(P,2)<size(P,1),
    message(sprintf('Reduced parameters from %d to %d',size(P,1),size(P,2)));
end