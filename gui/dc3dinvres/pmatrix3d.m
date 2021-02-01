function P = pmatrix3d(Nx,Ny,Divx,Divy,Cov,mincov)

% PMATRIX3D - Form p-matrix to combinate and delete cells
% P = pmatrix3d(Nx,Ny,Divx,Divy[,Cov,mincov])
% Nx,Ny .. number of cells in x/y direction
% Divx/Divy .. Vector of size of new cell in each layer
% optional
% Cov .. Coverage matrix of cells
% mincov .. minium coverage (delete cells below)

if nargin<3, error('3 input arguments are required!'); end
if nargin<4, Divy=Divx; end
if nargin<5, mincov=0; end
Nz=max(length(Divx),length(Divy));
P=speye(Nx*Ny*Nz);
streich=[];
for k = Nz:-1:1
    divx=Divx(k);
    divy=Divy(k);
    if divx*divy>1,
        for j = ceil(Ny/divy):-1:1
            for i = ceil(Nx/divx):-1:1
                nnx=min(divx,Nx-(i-1)*divx);
                nny=min(divy,Ny-(j-1)*divy);
                liste=[];
                start=(k-1)*Nx*Ny+(j-1)*divy*Nx+(i-1)*divx;
                for l = 1:nny,
                    liste=[liste start+(1:nnx)+(l-1)*Nx]; 
                end
                li=sum(P(:,liste),2);
                P(:,liste(1))=li;%/sqrt(sum(li));
                streich=[streich liste(2:end)];
            end
        end
    end
end
P(:,streich)=[];
if nargin>4,
    if nargin<6, mincov=0.4; end
    Voc=P'*Cov(:);
    P(:,find(Voc<mincov))=[];
end
