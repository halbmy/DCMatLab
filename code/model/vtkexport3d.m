function vtkexport3d(outfile,Mod,x,y,z,name,islog)

% VTKEXPORT3d - Export 3d model to VTK file
%   visualization toolkit - rectilinear grid
% vtkexport3d(modfile,Model[,name]) - grid or para model
% vtkexport3d(modfile,M,x,y,z,[name])

if nargin<7, islog=0; end
if nargin<6, name='Resistivity'; end
if nargin<2, error('Filename and Model array must be specified!'); end
if isstruct(Mod), % model structure
    if nargin>3, islog=y; end
    if nargin>2, name=x; end
    if isfield(Mod,'M'),
       if iscell(Mod.M), % Para Model
           [M,x,y,z]=mesch3dmodel(Mod);
       else % grid model
           M=Mod.M;x=Mod.x;y=Mod.y;z=Mod.z;
       end
    else
        error('Model corrupt!');
    end    
else % real grid model
    M=Mod;
    if nargin<5, z=0:size(M,3); end
    if nargin<4, y=0:size(M,2); end
    if nargin<3, x=0:size(M,1); end
end
lx=length(x);
ly=length(y);
lz=length(z);
if (lx~=size(M,1)+1)||(ly~=size(M,2)+1)||(lz~=size(M,3)+1),
  error('Array and Vector sizes do not match!'); end
fid=fopen(outfile,'w'); 
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'DC3dInvRes 3d grid model\n');
fprintf(fid,'ASCII\n');
fprintf(fid,'DATASET RECTILINEAR_GRID\n');
fprintf(fid,'DIMENSIONS %d %d %d\n',lx,ly,lz);
fprintf(fid,'X_COORDINATES %d float\n',lx);
fprintf(fid,'%g',x(1));
for i=2:lx, fprintf(fid,' %g',x(i)); end
fprintf(fid,'\n');
fprintf(fid,'Y_COORDINATES %d float\n',ly);
fprintf(fid,'%g',y(1));
for i=2:ly, fprintf(fid,' %g',y(i)); end
fprintf(fid,'\n');
fprintf(fid,'Z_COORDINATES %d float\n',lz);
fprintf(fid,'%g',-z(1));
for i=2:lz, fprintf(fid,' %g',-z(i)); end
fprintf(fid,'\n');
fprintf(fid,'CELL_DATA %d\n',prod(size(M)));
if islog,
    fprintf(fid,'SCALARS %s(log10) float\nLOOKUP_TABLE default\n',name);
    fprintf(fid,'%g\n',log10(M(:)));
else
    fprintf(fid,'SCALARS %s float\nLOOKUP_TABLE default\n',name);
    fprintf(fid,'%g\n',M(:));
end
fclose(fid);
