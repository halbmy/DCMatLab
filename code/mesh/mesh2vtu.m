
name='resistivity';
fid=fopen('test.vtk','w');
fprintf(fid,'<?xml version="1.0"?>\n');
fprintf(fid,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid,'  <UnstructuredGrid>\n');
fprintf(fid,'    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n',Mesh.nnodes,Mesh.ncells);
fprintf(fid,'      <CellData Scalars="%s">\n',name);
fprintf(fid,'        <DataArray type="Float32" Name="%s" format="binary">',name);
fwrite(fid,resS,'float32');
fprintf(fid,'        </DataArray>\n');
fprintf(fid,'      </CellData>\n');
fprintf(fid,'');
fprintf(fid,'');
fprintf(fid,'  </Piece>\n');
fprintf(fid,'  </UnstructuredGrid>\n');
fprintf(fid,'</VTKFile>\n');