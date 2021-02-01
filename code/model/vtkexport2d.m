function vtkexport2d(outfile,M,x,z,topo,coverage)

% VTKEXPORT2d - Export 3d model to VTK file
%   visualization toolkit - rectilinear grid
% vtkexport(modfile,M,x,z)
% vtkexport(modfile,M,x,z,topography)

if nargin<2, error('Filename and Model array must be specified!'); end
if nargin<7, y=0; end
if nargin<6, coverage=[]; end
if nargin<5, topo=[]; end
if nargin<4, z=0:size(M,2); end
if nargin<3, x=0:size(M,1); end
if (nargin==5)&&(size(topo,2)~=2), coverage=topo;topo=[]; end
lx=length(x);lz=length(z);ly=length(y);
z=-z;
if (lx~=size(M,1)+1)||(lz~=size(M,2)+1),
  error('Array and Vector sizes do not match!'); end
% y=z;z=0;ly=lz;lz=1;
fid=fopen(outfile,'w'); 
fprintf(fid,'# vtk DataFile Version 3.0\n');
fprintf(fid,'DC2dInvRes 2d grid model\n');
fprintf(fid,'ASCII\n');
if isempty(topo)||(size(topo,2)<2)||(length(unique(topo(:,2)))<2), % RECTILINEAR GRID
    if (nargin>4)&&(size(topo,2)>1), z=z+topo(1,2); end %add constant topo
    fprintf(fid,'DATASET RECTILINEAR_GRID\n');
    fprintf(fid,'DIMENSIONS %d %d %d\n',lx,ly,lz);
    fprintf(fid,'X_COORDINATES %d double\n',lx);
    fprintf(fid,'%g',x(1));
    for i=2:lx, fprintf(fid,' %g',x(i)); end
    fprintf(fid,'\n');
    fprintf(fid,'Y_COORDINATES %d double\n',ly);
    fprintf(fid,'%g',y(1));
    for i=2:ly, fprintf(fid,' %g',y(i)); end
    fprintf(fid,'\n');
    fprintf(fid,'Z_COORDINATES %d double\n',lz);
    fprintf(fid,'%g',z(1));
    for i=2:lz, fprintf(fid,' %g',z(i)); end
    fprintf(fid,'\n');
    if ly==1, ncells=(lx-1)*(lz-1); else ncells=(lx-1)*(ly-1)*(lz-1); end
else %topography->unstructured grid
    ncells=numel(M);npoints=lx*lz;
    fprintf(fid,'DATASET UNSTRUCTURED_GRID\n');
    fprintf(fid,'POINTS %d double\n',npoints);
    xz=interp1(topo(:,1),topo(:,2),x,'linear','extrap');
    for k=1:lz,
        for i=1:lx,
            fprintf(fid,'%g %g %g\n',x(i),0,xz(i)+z(k));
        end
    end
    fprintf(fid,'CELLS %d %d\n',ncells,ncells*5);
    for k=1:lz-1,
        for i=1:lx-1,
            ii=(k-1)*lx+i-1;
            fprintf(fid,'4 %d %d %d %d\n',ii,ii+1,ii+lx+1,ii+lx);
        end
    end
    fprintf(fid,'CELL_TYPES %d\n',ncells);
    fprintf(fid,'%d ',ones(1,ncells)*9); % VTK QUADS
end
fprintf(fid,'CELL_DATA %d\n',ncells);
fprintf(fid,'SCALARS Resistivity(log10) double 1\n');
fprintf(fid,'LOOKUP_TABLE default\n');
fprintf(fid,'%g ',log10(M(:)));
fprintf(fid,'\n');
if ~isempty(coverage)&&(numel(coverage)==numel(M)),
    fprintf(fid,'SCALARS Coverage double 1\n');
    fprintf(fid,'LOOKUP_TABLE default\n');
    fprintf(fid,'%g ',coverage(:));
    fprintf(fid,'\n');
end
fclose(fid);
