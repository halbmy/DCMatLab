function newM=blockify2dmodel(M)

% BLOCKIFYMODEL - Blockify model (logarithmical)
% newModel = blockifymodel(oldModel)

part=[1 1.5 2.5 4 6.5];
dlog=log10(max(M(:)))-log10(min(M(:)));
if dlog>3, part=[1 2 5]; end
if dlog<0.5, part=[1 1.4 1.9 2.7 3.7 5.2 7.2]; end
row=part;
while row(end)<5000,
    part=part*10;
    row=[row part];
end
row(end+1)=part(1)*10;
lrow=log(row);
lM=log(M(:));
newM=M;
for m=1:length(lM),
    [x,im]=min(abs(lM(m)-lrow));
    newM(m)=row(im(1));
end