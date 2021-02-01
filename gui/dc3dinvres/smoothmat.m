function C=smoothmat(Model)

% SMOOTHMAT - matrix for (horizontal) smoothness constraints
% C = smoothmat(Model)

nnz=0;
for k=1:length(Model.M), nnz=nnz+prod(size(Model.M{k})); end
C=spalloc(nnz,nnz,5*nnz);
old=0;
for k=1:length(Model.M),
    [nx,ny]=size(Model.M{k});
    if nx*ny>1,
        one=ones(nx*ny-1,1);
        Cx=spdiags([-one one],[0 1],nx*ny-1,nx*ny);
        Cx(nx:nx:end,:)=[];
        Cy=spdiags([-one one],[0 nx],nx*(ny-1),nx*ny);
        new=old+nx*ny;
        C(old+1:new,old+1:new)=Cx'*Cx+Cy'*Cy;
        old=new;
    end
end
    