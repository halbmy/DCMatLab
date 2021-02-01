function [dM,lam] = minvers(S,dR,INV,err,Model,dM0)

% MINVERS - Model inverse step
%           Solve S * dM = dR by way of matrix/tsvd/sirt inverse
% [dM,lambda] = minvers(S,dR[,OPT,err,Model,dM0])
% OPT = structure of possible fields:
%       method = 0(Matrix Inversion) Solution of
%       ((DS)'(DS) + lam C'C) dm = (DS)'D (d-F(m)) 
%                1 Truncated SVD-Inversion (needs lots of memory/time!)
%                2 SIRT
%                3 Truncated least squares (performs lam steps of 0)
%       auto   = automatic regularization choosing
%                0 interactive choosing
%                1 optimize lambda in every step
%                2 fixed lambda 
%       redu   = reduction of free parameters
%                0 no reduction
%                1 delete bad covered
%                2 combine cells
%                3 combine&delete 
%       mico   = minimum coverage for deleting bad cells 
%       weight = weighting matrix C(for matrix inversion)
%                0 - equal weighting
%                1 - smoothness 1st order
%                2 - smoothness 2nd order (Dirichlet boundaries)
%                3 - smoothness 2nd order (Neumann boundaries)
%                4 - smoothness 2nd order (Mixed boundaries)
%                5 - weighting by coverage
%       rbnorm = normalize distances (1)
%       rbzfak = factor for weighting z-derivatives (1)
%       the matrix D is constructed by diag(1./err)

global FIX
if nargin<3, INV=0; end
%     INV=struct('redu',0,'mitschicht',0,'method',2,'auto',0); end
if nargin<4, err=0.01; end
if nargin<6, dM0=0; end

if ~isfield(INV,'lam'), INV.lam=0; end
if ~isfield(INV,'redu'), INV.redu=0; end
if ~isfield(INV,'mitschicht'), INV.mitschicht=0; end
if ~isfield(INV,'method'), INV.method=2; end
if ~isfield(INV,'auto'), INV.auto=0; end
if ~isfield(INV,'weight'), INV.weight=0; end
if ~isfield(INV,'mico'), INV.mico=0; end
if ~isfield(INV,'rbnorm'), INV.rbnorm=0; end
if ~isfield(INV,'rbzfak'), INV.rbzfak=1; end
if (~isfield(INV,'glob'))||(INV.glob==0), dM0(:)=0; end
lam=0;
D=spdiags(1./log(1+err),0,length(err),length(err));
nm=numel(Model.M);
P=speye(nm);
if INV.redu==1, P=pmincov(S,0.4); end
if INV.redu==2,
%     P=1;
    %P=pmatrix3d(length(x)-1,length(y)-1,x); % !!!! to complete
end
%% Fix and compound
if isequal(numel(FIX),nm),
    fi=find(FIX(:)==-1);
    for i=1:length(fi), % fixed parameters
        P(:,find(P(fi(i),:)))=[];
    end
    unf=unique(FIX);unf=unf(unf>0); %compounds do not work
    for u=1:length(unf), % all cell compounds
        fi=find(FIX(:)==unf(u));fi1=zeros(size(fi));
        for i=1:length(fi), fi1(i)=find(P(fi(i),:)); end
        P(:,end+1)=sum(P(:,fi1),2);P(:,fi1)=[];
    end
end
%% compute model updates according to method
if ~isfield(INV,'method')||(INV.method==0),
    C=1;lam=INV.lam;
    maxlam=500;minlam=1;kuh=0.8;
    C=getcmatrix3d(Model,INV);
    if INV.mitschicht, % not yet implemented
    end
    % if INV.const==0, INV.lam=0; end
    if INV.lam>0, minlam=INV.lam; end
    if (INV.auto==2)&(INV.lam==0), INV.auto=2; end
    if INV.auto>1,
        if size(S,2)==length(dR), 
            dM=cglscdpt(S,dR,INV.lam,C,D,P,dM0(:));
        else
            dM=cglscdp(S,dR,INV.lam,C,D,P,dM0(:));
        end
    else
        % [dM,optlam]=invshiftcdp(S,dR,C,D,P,maxlam,minlam,kuh);
        kmax=fix(-log(maxlam/minlam)/log(kuh)+1);
        ak=maxlam*(kuh.^(0:kmax));
        global DM
        [DM,rho,eta,iter] = cglsparcdp(S,dR,C,D,P,ak);    
        save('rhoeta.mat','rho','eta','ak');
        if INV.auto==0,
            global warten
            m_choose;
        else %L-curve
            warten=lkurv(rho,eta,ak); 
        end
        dM=DM(:,warten);
        clear global DM warten
        if isfield(INV,'const')&&(INV.const), lam=ak(warten); end 
    end
elseif INV.method==1,
    dM=invtsvd(D*S,D*dR);        
elseif INV.method==2,
    dM=sirt(D*S,D*dR);
elseif INV.method==3,
    dM=tlscg(S,dR,INV.lam,D,P);
elseif INV.method==4,
    C=getcmatrix3d(Model,INV);
    dM=cglscdp(S,dR,INV.lam,C,D,P,dM0(:),0,INV.lam);
end

function P = pmincov(S,mincov)
% PMINCOV - For p-matrix from sensitivity by minimum coverage
% P = pmincov(S,mincov)
% mincov - minimum coverage (default 0.4)
if nargin<2, mincov=0.4; end
P=speye(size(S,2));
Cov=ones(size(S,2),1);
for i=1:length(Cov), Cov(i)=sum(abs(S(:,i))); end
P(:,find(COV<mincov))=[];