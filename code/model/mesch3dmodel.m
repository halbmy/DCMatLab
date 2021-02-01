function [M,x,y,z]=mesch3dmodel(Model)

% MESCHE - Regular mesh from Model structure
% [M,x,y,z] = meshe(Model)

nx=size(Model.M{1},1)*Model.nx(1);
ny=size(Model.M{1},2)*Model.ny(1);
nz=length(Model.z)-1;
x=(0:nx)*Model.dx+Model.x0;
y=(0:ny)*Model.dy+Model.y0;
z=Model.z;
M=ones(nx,ny,nz);
for k=1:nz,
    M(:,:,k)=Model.Bg(k);
end
for k=1:nz,
    nnx=Model.nx(k);
    nny=Model.ny(k);
    rx=floor(mod(nx,nnx)/2);
    rrx=mod(nx,nnx)-rx;
    ry=floor(mod(ny,nny)/2);
    rry=mod(ny,nny)-ry;
    for i=1:nnx,
        for j=1:nny,
            M(i+rx:nnx:end-rrx,j+ry:nny:end-rry,k)=Model.M{k};
        end
    end
end
