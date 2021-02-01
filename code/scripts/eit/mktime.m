% data=
if ~exist('isload')
    data=getToInvertData;
    invertData(data(1,:));
    Mesh=loadmesh('tmp/meshPara.bms');
    Ne=readinv2dfile('newDD16.dat');
    N=readinv2dfile('tmp.data');
    S=loadsens;
    isload=1;
end
abmn1=N.a*999+N.b*99+N.m*9+N.n;
abmn2=Ne.a*999+Ne.b*99+Ne.m*9+Ne.n;
[tf,loc]=ismember(abmn1,abmn2);
dat=data(:,loc)';
refdat=dat(:,1);refmod=load(fullfile('tmp','model_iter.7.vector'));
for i=1:size(dat,2), dat(:,i)=dat(:,i)./refdat; end
C=meshsmoothness(Mesh);
return
lam=0.1;
GI=inv(S'*S+lam*C)*S';
dmods=GI*log(dat);
for i=1:size(dmods,2),
    clf;
%     tripatchmod(Mesh,dmods(:,i),N,-0.05,0.05);
    model=refmod.*exp(dmods(:,i));
    clf;tripatchmod(Mesh,model,N,1,30);
    title(sprintf('Frame %3d: %5.2fs',i-1,(i-1)/13));
    Mov(i)=getframe(gcf);
end
% movie(Mov);
% movie2avi(Mov,'lungheart3.avi','fps',13);
return
pos=[-0.02 0.05];
cellmids=zeros(Mesh.ncells,2);
for i=1:Mesh.ncells, cellmids(i,:)=mean(Mesh.node(Mesh.cell(i,:),:)); end
di=sqrt(sum((cellmids-repmat(pos,Mesh.ncells,1)).^2,2));
[mi,nmi]=min(di);
plot(dmods(nmi,:));
return
dos('dc2dtreerun -B tmp.data')
model=load(fullfile('tmp','model_iter.6.vector'));
