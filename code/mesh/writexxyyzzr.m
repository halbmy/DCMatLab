function writexxyyzzr(filename,Mesh,res)

A=zeros(Mesh.ncells,7);
A(:,7)=res;
A(:,1:2)=Mesh.node(Mesh.cell(:,1),1:2);
A(:,3:4)=Mesh.node(Mesh.cell(:,2),1:2);
A(:,5:6)=Mesh.node(Mesh.cell(:,3),1:2);
fid=fopen(filename,'w');
fprintf(fid,'%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\r\n',A');
fclose(fid);