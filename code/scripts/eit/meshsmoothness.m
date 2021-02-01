function C=meshsmoothness(Mesh)

% MESHSMOOTHNESS - Calculate smoothness matrix

fb=find(Mesh.boundmarker==0);
C=spalloc(Mesh.ncells,Mesh.ncells,length(fb)*2);
for n=1:length(fb),
    i=Mesh.boundleft(fb(n));
    j=Mesh.boundright(fb(n));
    C(i,i)=C(i,i)+1;
    C(j,j)=C(j,j)+1;
    C(i,j)=-1;
    C(j,i)=-1;
end