function Mesh = hex2tet(Mod)

% HEX2TET - Convert 3d grid model to 
% Mesh = hex2tet(Mod)
% Mod - 3d grid model containing x,y,z and M
% Mesh - 3d tetraedral model containing node,cell and cellattr

Mesh=[];Mesh.dim=3;
[X,Y,Z]=ndgrid(Mod.x,Mod.y,Mod.z);
Mesh.node=[X(:) Y(:) Z(:)];
Mesh.cell=[];
nx=length(Mod.x);ny=length(Mod.y);nz=length(Mod.z);
nxy=nx*ny;
% uninode=[1 2 7 3;1 3 7 4;1 7 8 4;6 7 1 2;6 5 1 7;5 8 1 7];
uninode=[1 2 4 5;2 3 4 7;2 4 5 7;6 5 7 2;5 8 7 4];
% uninode=uninode(:,[1 2 4 3]);
nt=size(uninode,1);
aa=repmat(Mod.M(:)',nt,1);
Mesh.cellattr=aa(:);
Mesh.cell=zeros(numel(Mod.M)*size(uninode,1),4);
for k=1:nz-1,
    for j=1:ny-1,
        for i=1:nx-1,
            cellxy=[i i+1 i+nx+1 i+nx ]+ny*(j-1);
            cellp=[nxy*(k-1)+cellxy nxy*k+cellxy];
            nc=i-1+(nx-1)*(j-1)+(nx-1)*(ny-1)*(k-1);
            Mesh.cell(nc*nt+(1:nt),:)=cellp(uninode);
        end
    end
end
Mesh.nnodes=size(Mesh.node,1);
Mesh.ncells=size(Mesh.cell,1);
Mesh.nbounds=0;
Mesh.nodemarker=zeros(Mesh.nnodes,1);
Mesh.cellnodes=ones(Mesh.ncells,1)*4;