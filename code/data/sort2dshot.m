function NN=sort2dshot(N,islinear)

% SORT2DSHOT - sort 2d electrode points clock-wise
% Nnew = sort2dshot(Nold) sorts by angle to midpoint
% sort2dshot(N,1) sorts by x

if nargin<2, islinear=0; end

if islinear,
    wi=N.pos(:,1);
else
    mi=mean(N.pos);mi(end)=mi(end)+0.001;
    wi=atan2(N.pos(:,2)-mi(2),N.pos(:,1)-mi(1));
end
[swi,ii,jj]=unique(wi);
% jj=zeros(size(ii));jj(ii)=1:length(jj);
NN=N;
NN.pos=N.pos(ii,:);
for i=1:length(N.ns),
    NN.ns{i}=jj(N.ns{i});
    NN.nx{i}=jj(N.nx{i});
end
