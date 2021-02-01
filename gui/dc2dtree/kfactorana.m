function k = kfactorana(N,ana)

% KFACTORANA - calculate k(configuration) factor
%              analytically (after Weidelt&Weller)
% k = kfactorana(N)

nel=size(N.elec,1);
ana1=[ana(1:end-1);flipud(ana(:))];
anan=interp1(1:length(ana1),ana1,(1:nel-1)'*length(ana)*2/nel);
am=max(N.a,N.m)-min(N.a,N.m);
an=max(N.a,N.n)-min(N.a,N.n);
bm=max(N.b,N.m)-min(N.b,N.m);
bn=max(N.b,N.n)-min(N.b,N.n);
mid=mean(N.elec);
rad=mean(sqrt(sum((N.elec-repmat(mid,nel,1)).^2,2)));
k=rad./(anan(am)-anan(an)-anan(bm)+anan(bn));