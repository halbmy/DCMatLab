function C=meshsmoothness(Mesh,po)

% MESHSMOOTHNESS - Calculate smoothness matrix

if nargin<2, po=0; end
% fb=find(Mesh.boundmarker==0);
if Mesh.dim>2, error('Dimension 3 not yet implemented!'); end
fb=find((Mesh.boundleft>0)&(Mesh.boundright>0));
C=spalloc(Mesh.ncells,Mesh.ncells,length(fb)*2);
for n=1:length(fb),
    i=Mesh.boundleft(fb(n));
    j=Mesh.boundright(fb(n));
    ve=diff(Mesh.node(Mesh.bound(fb(n),:),:));
    val=(1-abs(ve(1))/norm(ve))^po;
%     fprintf('%.3f ',[ve val]);fprintf('\n');
    C(i,i)=C(i,i)+val;
    C(j,j)=C(j,j)+val;
    C(i,j)=-val;
    C(j,i)=-val;
end