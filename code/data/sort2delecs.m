function NN=sort2delecs(N,islinear)

% SORT2DELECS - sort 2d electrode points clock-wise
% Nnew = sort2delecs(Nold) sorts by angle to midpoint
% sort2delec(N,1) sorts by x

if nargin<2, islinear=0; end

if islinear,
    wi=N.elec(:,1);
else
    mi=mean(N.elec);mi(end)=mi(end)+0.001;
    wi=atan2(N.elec(:,2)-mi(2),N.elec(:,1)-mi(1));
end
[swi,ii]=sort(wi);
NN=N;
NN.elec=N.elec(ii,:);
jj=zeros(size(ii));jj(ii)=1:length(jj);
fi=find(N.a);NN.a(fi)=jj(N.a(fi));
fi=find(N.b);NN.b(fi)=jj(N.b(fi));
fi=find(N.m);NN.m(fi)=jj(N.m(fi));
fi=find(N.n);NN.n(fi)=jj(N.n(fi));
