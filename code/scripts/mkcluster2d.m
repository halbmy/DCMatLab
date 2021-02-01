global x z M
xm=x(1:end-1)+diff(x)/2;
zm=z(1:end-1)+diff(z)/2;
xm=(1:size(M,1))';zm=(1:size(M,2))';
[X,Z]=ndgrid(xm,zm);
A=[X(:) Z(:) log10(M(:))];
% A=log10(M(:));
for i=1:size(A,2),
    mi=min(A(:,i));ma=max(A(:,i));
    A(:,i)=(A(:,i)-mi)/(ma-mi);
end
A(:,3)=A(:,3)*2;
distA = pdist(A,'euclid');
% Compute Clusters based on Ward...
linkA = linkage(distA,'complete');
% [H,T] = dendrogram(linkA,0);
% Plot Distances to decide # of clusters...
fprintf('Plotting distances...\n\n')
figure(1);clf;
subplot(2,1,1);
mal=struct('clog',1,'cauto',0,'cmin',min(M(:)),'cmax',max(M(:)));
patch2dmodel(x,z,M,mal);
subplot(2,1,2);
[mx,my] = size(linkA);
plot(linkA(mx-50:mx,3),'o-')
snum=inputdlg('Wieviele?');
numclust=str2num(snum{1});
% for numclust=8,%1:10,
    CM = cluster(linkA,numclust);
    nM=M;cl=[];
    for i=1:numclust,
        fi=find(CM==i);
        cl(i)=median(M(fi));
        nM(fi)=cl(i);
    end
    patch2dmodel(x,z,nM,mal);
%     pause(1.0);
% end