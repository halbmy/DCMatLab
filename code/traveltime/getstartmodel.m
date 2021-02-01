function velocity=getstartmodel(Mesh,vmin,vmax,zz)

% GETSTARTMODEL - Get Start model
% velocity=getstartmodel(Mesh,v) - uniform
% velocity=getstartmodel(Mesh,vmin,vmax) - gradient model
% velocity=getstartmodel(Mesh,vmin,vmax,zlayer) - 2-layered
% velocity=getstartmodel(Mesh,Shot) - determine from data

if nargin<4, zz=0; end
if nargin<3, vmax=0; end
if nargin<2, error('Must specify mesh and Shot or min/max!'); end
if (~isstruct(Mesh))||(~isfield(Mesh,'node'))||(~isfield(Mesh,'cell')),
    error('1st argument must be a valid Mesh!'); end
if isstruct(vmin), %automatic detection
    error('not yet implemented!');
    %vmin=
    %vmax=
end
% vmin and vmax specified
velocity=ones(Mesh.ncells,1)*vmin;
if vmax>0, % not uniform
    cellmids=zeros(Mesh.ncells,Mesh.dim);
    for i=1:Mesh.ncells, cellmids(i,:)=mean(Mesh.node(Mesh.cell(i,:),:)); end
    A=ones(Mesh.ncells,2);A(:,2)=cellmids(:,1);
    ab=A\cellmids(:,Mesh.dim); % remove linear trend by regression
    cellmids(:,Mesh.dim)=cellmids(:,Mesh.dim)-ab(1)-ab(2)*cellmids(:,1);
    if zz>0, % 2-layered case
        velocity(abs(cellmids(:,Mesh.dim))>zz)=vmax;
    else, % gradient model
%         deprel=abs(cellmids(:,Mesh.dim)-max(Mesh.node(:,Mesh.dim)));
        deprel=abs(cellmids(:,Mesh.dim)-max(cellmids(:,Mesh.dim)));
        deprel=deprel/max(deprel);
        velocity=(vmax/vmin).^deprel*vmin;
        % velocity=vmin+vmax*abs(cellmids(:,2))/max(abs(cellmids(:,2)));
    end
end