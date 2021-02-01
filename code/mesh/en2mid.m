function [x,y,z]=en2mid(ELE,NODE)

% EN2MID - Element-Node to midpoint conversion
% [x,y,z] = en2mid(Elements,Nodes)
% [x,y,z] = en2mid(Mesh)
% xyz = en2mid(...)

nel=size(ELE,1);
x=zeros(nel,1);y=x;z=x;
for e=1:nel,
    x(e)=mean(NODE(ELE(e,1:4),1));    
    y(e)=mean(NODE(ELE(e,1:4),2));    
    z(e)=mean(NODE(ELE(e,1:4),3));    
end
if nargout==1, %only 1 structure
    x=[x y z];
end