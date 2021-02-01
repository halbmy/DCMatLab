function savemesh(Mesh,meshname)

% SAVEMESH - Save Mesh to binary mesh (.bms) file
% savemesh(Mesh,meshfile)

if nargin<2, meshname='tmp\meshPara.bms'; end
if ~strfind(meshname,'.bms'), meshname=[meshname '.bms']; end

zahl='int32';wert='double';
fid=fopen(meshname,'w');
if fid<0, error('Could not open file!'); end
try,
    fwrite(fid,Mesh.dim,zahl);
    fwrite(fid,zeros(127,1),zahl);%vertinfo
    fwrite(fid,Mesh.nnodes,zahl);
    fwrite(fid,Mesh.node',wert);
    fwrite(fid,Mesh.nodemarker,zahl);
    fwrite(fid,zeros(127,1),zahl);%cellinfo
    fwrite(fid,Mesh.ncells,zahl);
    fwrite(fid,Mesh.cellnodes,zahl);
    fwrite(fid,Mesh.cell'-1,zahl);
    fwrite(fid,Mesh.cellattr,wert);
    fwrite(fid,zeros(127,1),zahl);%boundinfo
    fwrite(fid,Mesh.nbounds,zahl);
    fwrite(fid,Mesh.boundnodes,zahl);
    fwrite(fid,Mesh.bound'-1,zahl);
    fwrite(fid,Mesh.boundmarker,zahl);
    fwrite(fid,Mesh.boundleft-1,zahl);
    fwrite(fid,Mesh.boundright-1,zahl);
    fclose(fid);
catch
    display(lasterr);
    fclose(fid);
end