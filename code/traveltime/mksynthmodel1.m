xmax=300; % profile length
zmax=100; % model depth
dg=5; % geophone spreading
dsdg=2; % 2 shots per geophone
z1=10;
Poly=[];
Poly.node=[0 0;xmax 0;xmax z1;xmax zmax;0 zmax;0 z1];
Poly.node(:,3)=0;
Poly.edge=[1 2;2 3;3 4;4 5;5 6;6 1;6 3];
Poly.edge(:,3)=0;
Poly.region=[1 1 1;1 zmax-1 2];
%%
pos=50;len=20;dlen2=5;dep=10;ln=size(Poly.node,1);
Poly.node(end+1:end+4,1:2)=[pos z1;pos+len z1;pos+len-dlen2 z1+dep;pos+dlen2 z1+dep];
Poly.edge(end,2)=ln+1;
Poly.edge(end+1:end+5,1:2)=ln+[1 2;2 3;3 4;4 1;2 0];
Poly.edge(end,2)=3;
Poly.region(end+1,:)=[pos+len/2 z1+1 size(Poly.region,1)+1];
%%
pos=120;len=50;dlen2=10;dep=20;ln=size(Poly.node,1);
Poly.node(end+1:end+4,1:2)=[pos z1;pos+len z1;pos+len-dlen2 z1+dep;pos+dlen2 z1+dep];
Poly.edge(end,2)=ln+1;
Poly.edge(end+1:end+5,1:2)=ln+[1 2;2 3;3 4;4 1;2 0];
Poly.edge(end,2)=3;
Poly.region(end+1,:)=[pos+len/2 z1+1 size(Poly.region,1)+1];
%%
pos=220;len=10;dlen2=2;dep=5;ln=size(Poly.node,1);
Poly.node(end+1:end+4,1:2)=[pos z1;pos+len z1;pos+len-dlen2 z1+dep;pos+dlen2 z1+dep];
Poly.edge(end,2)=ln+1;
Poly.edge(end+1:end+5,1:2)=ln+[1 2;2 3;3 4;4 1;2 0];
Poly.edge(end,2)=3;
Poly.region(end+1,:)=[pos+len/2 z1+1 size(Poly.region,1)+1];
show2dpoly(Poly);
%%
Poly.node(1,3)=-99;
dg=5;gpos=(Poly.node(1,1)+dg:dg:Poly.node(2,1)-dg)';
Poly.node(2,3)=-99;
gpos(:,2)=0;gpos(:,3)=-99;Poly.node=[Poly.node;gpos];
show2dpoly(Poly);
%%
Poly.node(:,2)=-Poly.node(:,2);Poly.region(:,2)=-Poly.region(:,2);
writepoly2d('synth.poly',Poly);
%%
epos=gpos;epos(:,2)=epos(:,2)-1;epos(:,3)=0;
Poly1=Poly;Poly1.node=[Poly.node;epos];
writepoly2d('synth.poly',Poly1);
system('dctriangle -v -q34 synth.poly');
Mesh=loadmesh('synth.bms');
tripatchmod(Mesh)
%%
Shot=[];Shot.pos=Mesh.node(Mesh.nodemarker==-99,1:2);
Shot.ns={};Shot.nx={};Shot.tt={};Shot.t=[];
ng=sum(Mesh.nodemarker==-99);
all=(1:ng)';
i=1;
while i<ng,    
    Shot.ns{end+1}=i;
    Shot.nx{end+1}=all;
    Shot.nx{end}(Shot.nx{end}==i)=[];    
    Shot.tt{end+1}=zeros(size(Shot.nx{end}));
    Shot.t=[Shot.t;Shot.tt{end}];
    i=i+2;
end
%%
veli=[500 1800 1000 1000 1000]';
tripatchmod(Mesh,veli(Mesh.cellattr));
%%
[W,Shot.t]=waymatrix(Mesh,Shot,veli(Mesh.cellattr));
plotshot(Shot,Shot.t);
savesgtfile('synth.sgt',Shot);