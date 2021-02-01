function nM=cluster2dmodel(M,M2)

% CLUSTER2DMODEL - Cluster 2d model
% clustered_model = cluster2dmodel(M1[,M2]);

xm=(1:size(M,1))';zm=(1:size(M,2))';
[X,Z]=ndgrid(xm,zm);
A=[X(:) Z(:) log10(M(:))];
if (nargin>1)&&isequal(size(M),size(M2)), A=[A M2(:)]; end
for i=1:size(A,2),
    mi=min(A(:,i));ma=max(A(:,i));
    A(:,i)=(A(:,i)-mi)/(ma-mi);
end
A(:,3)=A(:,3)*2;
distA = pdist(A,'euclid');
linkA = linkage(distA,'complete');
[mx,my] = size(linkA);
plot(flipud(linkA(mx-20:mx,3)),'o-');
grid on
snum=inputdlg('Number of clusters?');
numclust=str2double(snum{1});
CM = cluster(linkA,numclust);
nM=M;cl=[];
for i=1:numclust,
    fi=find(CM==i);
    cl(i)=median(M(fi));
    nM(fi)=cl(i);
end
