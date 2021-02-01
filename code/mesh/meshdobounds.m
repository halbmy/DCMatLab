function Mesh=meshdobounds(Mesh)

if Mesh.dim==3,
    Mesh.nbounds=Mesh.ncells*6;
    Mesh.bound=zeros(Mesh.nbounds,2);
    for i=1:Mesh.ncells,
        Mesh.bound((i-1)*6+(1:6),:)=[Mesh.cell(i,[1 1 1 2 2 3])' Mesh.cell(i,[2 3 4 3 4 4])'];
    end
end
if Mesh.dim==2,
    Mesh.nbounds=Mesh.ncells*3;
    Mesh.bound=zeros(Mesh.nbounds,2);
    for i=1:Mesh.ncells,
        Mesh.bound((i-1)*3+(1:3),:)=[Mesh.cell(i,[1 1 2])' Mesh.cell(i,[2 3 3])'];
    end
    for i=1:Mesh.nbounds,
        bb=Mesh.bound(i,:); % edge nodes
        [fi,jj]=find(Mesh.cell==bb(1));
        for j=2:length(bb),
            [ii,jj]=find(Mesh.cell==bb(j));
            fi=intersect(fi,ii);
        end
        Mesh.boundleft(i)=fi(1);
        if length(fi)>1, Mesh.boundright(i)=fi(2); end
    end
    Mesh.nbounds=size(Mesh.bound,1);
    Mesh.boundnodes=ones(Mesh.nbounds,1)*2;
    Mesh.boundmarker=zeros(Mesh.nbounds,1);
end
