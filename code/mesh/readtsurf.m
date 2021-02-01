filename='bsp_0b.ts';
fid=fopen(filename,'r');
zeile=fgetl(fid);
while isstr(zeile)&&(~strcmp(zeile(1:6),'HEADER')),
    zeile=fgetl(fid);
end
% do something with header
while isstr(zeile)&&(~strcmp(zeile(1),'}')),
    zeile=fgetl(fid);
end
zeile=fgetl(fid);
% read other information about properties
while isstr(zeile)&&((length(zeile)<5)||(~strcmp(zeile(1:5),'TFACE'))),
    zeile=fgetl(fid);
end
zeile=fgetl(fid);
Mesh=[];Mesh.cell=[];
while isstr(zeile)&&(~strcmp(zeile,'END')),
    if length(zeile>3)&&strcmp(zeile(1:4),'VRTX'),
        aa=str2num(zeile(5:end));
        Mesh.node(aa(1),1:length(aa)-1)=aa(2:end);
    end
    if length(zeile>4)&&strcmp(zeile(1:5),'PVRTX'),
        aa=str2num(zeile(6:end));
        Mesh.node(aa(1),1:3)=aa(2:4);
        if length(aa)>4, Mesh.patt(aa(1))=aa(5); end
        if length(aa)>5, Mesh.patt2(aa(1))=aa(6); end
    end
    if length(zeile>3)&&strcmp(zeile(1:4),'TRGL'),
        aa=str2num(zeile(5:end));
        Mesh.cell(end+1,1:length(aa))=aa(1:end);
    end
    zeile=fgetl(fid);
end
fclose(fid);
clf;patch('Vertices',Mesh.node,'Faces',Mesh.cell,'FaceColor',[1 1 1],'EdgeColor',[0 0 0]);
% axis equal tight
Mesh.ncells=size(Mesh.cell,1);
Mesh.nnodes=size(Mesh.node,1);
Mesh.patt(Mesh.patt==-99999)=0;
if isfield(Mesh,'patt')&&(length(Mesh.patt)==Mesh.nnodes),
    Mesh.att=zeros(Mesh.ncells,1);
    for i=1:Mesh.cell, Mesh.att(i)=mean(Mesh.patt(i,:)); end
end