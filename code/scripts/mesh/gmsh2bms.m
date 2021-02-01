%% gmsh2bms - fwagner@GFZ (04/2012)
% reads in a gmsh-file (*.msh) using load_gmsh.m and exports it to the binary format used in
% BERT (*.bms) using savemesh.m (DCMatlab required)

function gmsh2bms(gmesh_name,bmesh_name,dim)

if isempty(strfind(gmesh_name,'.msh')), gmesh_name=[gmesh_name '.msh']; end

% load .msh-file
GMesh=load_gmsh(gmesh_name);
BMesh=[];

switch dim
    case 2
        % Nodes
        BMesh.dim=dim;
        BMesh.node=GMesh.POS;
        BMesh.node(:,3)=[];
        BMesh.nnodes=size(BMesh.node,1);
        
        % Cells
        BMesh.cell=GMesh.TRIANGLES(1:GMesh.nbTriangles,1:dim+1);
        BMesh.ncells=size(BMesh.cell,1);
        BMesh.cellnodes=ones(BMesh.ncells,1)*size(BMesh.cell,2);
        
        % Boundaries and regions
        % Electrode definition
        BMesh.nodemarker=zeros(BMesh.nnodes,1);
        for i = 1:GMesh.nbPoints, BMesh.nodemarker(GMesh.POINTS(i,1),1)=GMesh.POINTS(i,2)*(-1); end
        % Parameter domain
        BMesh.cellattr=GMesh.TRIANGLES(1:GMesh.nbTriangles,dim+2);
        % Boundary conditions
        BMesh.nbounds=GMesh.nbLines;
        BMesh.boundnodes=ones(BMesh.nbounds,1)*dim;
        BMesh.bound=GMesh.LINES(1:BMesh.nbounds,1:dim);
        BMesh.boundmarker=GMesh.LINES(1:BMesh.nbounds,dim+1);
        BMesh.boundmarker(find(BMesh.boundmarker==2))=BMesh.boundmarker(find(BMesh.boundmarker==2))*-1;
        BMesh.boundmarker(find(BMesh.boundmarker==3))=BMesh.boundmarker(find(BMesh.boundmarker==3))./-3;
        BMesh.boundleft=zeros(BMesh.nbounds,1);
        BMesh.boundright=zeros(BMesh.nbounds,1);
            
    case 3
        % Nodes
        BMesh.dim=dim;
        BMesh.node=GMesh.POS;
        BMesh.nnodes=size(BMesh.node,1);
        
        % Cells
        BMesh.cell=GMesh.TETS(1:GMesh.nbTets,1:dim+1);
        BMesh.ncells=size(BMesh.cell,1);
        BMesh.cellnodes=ones(BMesh.ncells,1)*size(BMesh.cell,2);
        
        % Boundaries and regions
        % Electrode definition
        BMesh.nodemarker=zeros(BMesh.nnodes,1);
        for i = 1:GMesh.nbPoints, BMesh.nodemarker(GMesh.POINTS(i,1),1)=GMesh.POINTS(i,2)*(-1); end
        % Parameter domain
        BMesh.cellattr=GMesh.TETS(1:GMesh.nbTets,dim+2);
        % Boundary conditions
        BMesh.nbounds=GMesh.nbTriangles;
        BMesh.boundnodes=ones(BMesh.nbounds,1)*dim;
        BMesh.bound=GMesh.TRIANGLES(1:BMesh.nbounds,1:dim);
        BMesh.boundmarker=GMesh.TRIANGLES(1:BMesh.nbounds,dim+1);
        BMesh.boundmarker(find(BMesh.boundmarker==2))=BMesh.boundmarker(find(BMesh.boundmarker==2))*-1;
        BMesh.boundmarker(find(BMesh.boundmarker==3))=BMesh.boundmarker(find(BMesh.boundmarker==3))./-3;
        BMesh.boundleft=zeros(BMesh.nbounds,1);
        BMesh.boundright=zeros(BMesh.nbounds,1);
        
    otherwise 
        disp('ERROR: Dimension not specified.')
end

%%
savemesh(BMesh,[bmesh_name '.bms']);
savevtkmesh(BMesh,[bmesh_name '.vtk']);
% mesh2vtk('Gitter1.vtk',BMesh); % ganz ähnlich
