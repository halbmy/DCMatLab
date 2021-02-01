function Mod=loadhexvtk(filename)

% LOADHEXVTK - Load hex vtk file
% Mod = loadhexvtk(filename)

Mod=[];
fid=fopen(filename);
for i=1:5, zeile=fgetl(fid); end
dims=sscanf(zeile,'%*s%f%f%f');
zeile=fgetl(fid);
xdim=sscanf(zeile,'%*s%d%*s');
Mod.x=fscanf(fid,'%f',xdim);
zeile=fgetl(fid);
zeile=fgetl(fid);
ydim=sscanf(zeile,'%*s%d%*s');
Mod.y=fscanf(fid,'%f',ydim);
zeile=fgetl(fid);
zeile=fgetl(fid);
zdim=sscanf(zeile,'%*s%d%*s');
Mod.z=fscanf(fid,'%f',zdim);
Mod.M=zeros(xdim-1,ydim-1,zdim-1);
zeile=fgetl(fid);
zeile=fgetl(fid);
ncells=sscanf(zeile,'%*s%d');
zeile=fgetl(fid);
zeile=fgetl(fid);
data=fscanf(fid,'%f',ncells);
fclose(fid);
if numel(Mod.M)==length(data), Mod.M(:)=data(:); end