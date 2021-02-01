function C = getcmatrix3d(Model,INV)

% GETCMATRIX - Get constraint matrix
% L = getcmatrix3d(Model,INV);
% Model - Model structure (Para or grid)
% INV - structure of options
%       weight = weighting matrix C(for matrix inversion)
%                0 - equal weighting
%                1 - smoothness 1st order
%                2 - smoothness 2nd order
%                3 - smoothness 2nd order (Neumann boundaries)
%                4 - smoothness 2nd order (Mixed boundaries)
%                5 - weighting by coverage
%       rbnorm = normalize distances (1)
%       rbzfak = factor for weighting z-derivatives (1)

global XX YY ZZ
if nargin<1, error('Model input required for input'); end
if nargin<2, INV=struct('default',1); end
%% complement fields
if ~isfield(INV,'weight'), INV.weight=0; end
if ~isfield(INV,'mico'), INV.mico=0; end
if ~isfield(INV,'rbnorm'), INV.rbnorm=1; end
if ~isfield(INV,'rbzfak'), INV.rbzfak=1; end

%% decoupling matrices
Wx=1;Wy=1;Wz=1;
if ~iscell(Model.M),
    simo=size(Model.M);
    if isequal(size(XX),simo-[1 0 0]),
        Wx=spdiags(1-XX(:),0,numel(XX),numel(XX));
    end
    if isequal(size(YY),simo-[0 1 0]),
        Wy=spdiags(1-YY(:),0,numel(YY),numel(YY));
    end
    if isequal(size(ZZ),simo-[0 0 1]),
        Wz=spdiags(1-ZZ(:),0,numel(ZZ),numel(ZZ));
    end
end
%% determine regularization matrices
C=1; % 0th order damping
if INV.weight==1, % first order smoothness
    if iscell(Model.M),
        C=smoothmat(Model);
        C=C+speye(size(C))*0.01;
    else
        [C,Cx,Cy,Cz]=smooth3d1st(Model.x,Model.y,Model.z,INV.rbnorm,INV.rbzfak);
        Cx=Wx*Cx;
        Cy=Wy*Cy;
        Cz=Wz*Cz;
        Cxyz=[Cx;Cy;Cz*sqrt(INV.rbzfak)];
        if isfield(INV,'blocky')&&(INV.blocky>0),
%             MM=log(Model.M(:));
%             sm=abs(Cxyz*MM);su2=sum(sm.^2);sua=sum(sm);
%             wxz=ones(size(sm));fi=find(sm);
%             if ~isempty(fi), wxz(fi)=su2/sua./sm(fi); end
%             wxz(wxz>1)=1;
            wxz=irls(Cxyz*log(Model.M(:)),0,1);
            Cxyz=spdiags(wxz,0,length(wxz),length(wxz))*Cxyz;
        end
        C=Cxyz'*Cxyz;
    end
end
if ismember(INV.weight,2:4), % 2nd order smoothness
    if iscell(Model.M),
        C=smoothmat2(Model);
        C=C+speye(size(C));
    else
        C1=smooth3d2nd(Model.x,Model.y,Model.z,INV.weight-2,INV.rbnorm,INV.rbzfak);
        if isfield(INV,'blocky')&&(INV.blocky>0),
            MM=log(Model.M(:));
            sm=abs(C1*MM);su2=sum(sm.^2);sua=sum(sm);
            wxz=ones(size(sm));fi=find(sm);
            if ~isempty(fi), wxz(fi)=su2/sua./sm(fi); end
            wxz(wxz>1)=1;
            C1=spdiags(wxz,0,length(wxz),length(wxz))*C1;
        end
        C=C1'*C1;
    end
end
if INV.weight==5, % Weighting by coverage
    global Cov
    if isempty(Cov),
        global S
        Cov=ones(size(S,2),1);
        for i=1:length(Cov), Cov(i)=sum(abs(S(:,i))); end
    end
    C=spdiags(1./Cov(:),0,prod(size(Cov)),prod(size(Cov)));
end
if INV.weight==6, % my favourite mix against stripes
    C1=smooth3d2nd(Model.x,Model.y,Model.z,1,INV.rbnorm,INV.rbzfak);
    C=(C1'*C1+smooth3d1st(Model.x,Model.y,Model.z,INV.rbnorm,INV.rbzfak))*0.5+speye(numel(Model.M))*0.01;
end