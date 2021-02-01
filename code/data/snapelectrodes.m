% function NN = snapelectrodes(N,ic)
NN=N;
dx=2;
% relec=round(N.elec/dx)*dx;
% NN.elec=unique(relec,'rows');
distA = pdist(N.elec,'euclidean');
linkA = linkage(distA,'average');
[mx,my] = size(linkA);
% plot(linkA(:,3),'o');grid on
ic=750;
if ic>0,
    CM = cluster(linkA,ic);
    NN.elec=zeros(ic,size(N.elec,2));
    for i=1:ic, NN.elec(i,:)=mean(N.elec(CM==i,:),1); end
    dist=sqrt(sum((N.elec-NN.elec(CM,:)).^2,2));
%     hist(dist,30)
end
% return
% else
%     for ic=size(N.elec,1):-1:1;
fi=find(N.a);NN.a(fi)=CM(N.a(fi));
fi=find(N.b);NN.b(fi)=CM(N.b(fi));
fi=find(N.m);NN.m(fi)=CM(N.m(fi));
fi=find(N.n);NN.n(fi)=CM(N.n(fi));
plot(N.elec(:,1),N.elec(:,2),'b.',NN.elec(:,1),NN.elec(:,2),'ro');
[pp,ff,ee]=fileparts(datfile);
saveinv3dfile(strrep(datfile,ee,['snap' num2str(ic) '.dat']),NN);