Mesh=loadmesh('ray\ray3.bms'); %four electrode mesh
Mesh.node(:,end)=Mesh.node(:,end)/3;
cellmids=zeros(Mesh.ncells,2);
for i=1:Mesh.ncells, cellmids(i,:)=mean(Mesh.node(Mesh.cell(i,:),:)); end
cellrad=sqrt(sum(cellmids.^2,2));
Mesh.cellattr=ones(Mesh.ncells,1);m0=Mesh.cellattr;
Mesh.cellattr(cellrad<0.5)=0.5;msynth=Mesh.cellattr;
tripatchmod(Mesh,msynth);pause(0.5);
Ws=waymatrix(Mesh);
ts=Ws*Mesh.cellattr.*(1+randn(size(t))*0.03);
Mesh.cellattr=ones(Mesh.ncells,1);
W=waymatrix(Mesh);t=W*Mesh.cellattr;
appsl=ts./sum(W,2);
coverage=sum(W);
C=meshsmoothness(Mesh);
D=0.03;lam=1;it=0;
fprintf('It. %d: RMS=%.2f\n',it,rms(ts,t));
for it=1:5,
    deltaData=ts-t;
    deltaModel=0;Mesh.cellattr-m0;
    dm=cglscdp(W,deltaData,0.01,C,D,1,deltaModel);
    Mesh.cellattr=Mesh.cellattr+dm;
    W=waymatrix(Mesh);    
    t=W*Mesh.cellattr;
    fprintf('It. %d: RMS=%.2f\n',it,rms(ts,t));
    tripatchmod(Mesh);pause(0.5);
end