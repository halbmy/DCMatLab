function [newMesh,sumvol]=extractvtk(Mesh,att,xl,yl,zl)

if nargin<2, att=100; end
if nargin<3, xl=-Inf; end
if nargin<4, yl=-Inf; end
if nargin<5, zl=-Inf; end
newMesh=Mesh;
cellmids=zeros(Mesh.ncells,3);
for i=1:Mesh.ncells, cellmids(i,:)=mean(Mesh.node(Mesh.cell(i,:),:)); end
lala=(Mesh.cellattr>att)&(cellmids(:,1)>xl)&(cellmids(:,2)>zl)&(cellmids(:,3)>zl);
if length(att)>1, lala=lala&(Mesh.cellattr<att(2)); end
if length(xl)>1, lala=lala&(cellmids(:,1)<xl(2)); end
if length(yl)>1, lala=lala&(cellmids(:,2)<yl(2)); end
if length(zl)>1, lala=lala&(cellmids(:,3)<zl(2)); end
fi=find(lala);    
A=zeros(3);
for i=1:Mesh.ncells,
    nodes=Mesh.cell(i,:);
    A(1,:)=diff(Mesh.node(nodes([1 2]),:));
    A(2,:)=diff(Mesh.node(nodes([1 3]),:));
    A(3,:)=diff(Mesh.node(nodes([1 4]),:));
    volume(i)=det(A)/6;
end
sumvol=sum(volume(fi));
newMesh.cell=Mesh.cell(fi,:);
newMesh.cellattr=log10(Mesh.cellattr(fi,:));
newMesh.ncells=size(newMesh.cell,1);
savevtkmesh(newMesh,'test.vtk');