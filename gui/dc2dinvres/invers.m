function [dM,lam] = invers(S,dR,INV,err,dM0)

% INVERSION - Solve S * dM = dR
% [dM,lam] = inversion(S,dR,OPT,err,dModel)
% OPT = structure of possible fields:
%       method = 0 Matrix Inversion
%                1 Truncated SVD-Inversion
%                2 SIRT
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
%       lam    = regularization parameter
%         if auto=0 or 1 lam = minimum lambda to test

global Mod FIX XX ZZ
if nargin<3, 
    INV=struct('redu',0,'mitschicht',0,'method',0,'auto',0); end
if nargin<4, err=0.1; end
if nargin<5, dM0=0; end

if ~isfield(INV,'lam'), INV.lam=0; end
if ~isfield(INV,'redu'), INV.redu=0; end
if ~isfield(INV,'mitschicht'), INV.mitschicht=0; end
if ~isfield(INV,'method'), INV.method=0; end
if ~isfield(INV,'auto'), INV.auto=0; end
if ~isfield(INV,'weight'), INV.weight=0; end
if ~isfield(INV,'mico'), INV.mico=0; end
if ~isfield(INV,'rbnorm'), INV.rbnorm=1; end
if ~isfield(INV,'rbzfak'), INV.rbzfak=1; end
if ~isfield(INV,'silent'), INV.silent=0; end
if ~isfield(INV,'maxiter'), INV.maxiter=1000; end

if (~isfield(INV,'glob'))||(INV.glob==0), dM0(:)=0; end
D=spdiags(1./log(1+err),0,length(err),length(err));
if ~isfield(INV,'method')||(INV.method==0),
    sim=(length(Mod.x)-1)*(length(Mod.z)-1);
    C=1;P=speye(sim);lam=INV.lam;
    maxlam=1000;minlam=1;kuh=0.8;
    if INV.redu==1, P=pmincov(S,INV.mico); end
    if INV.redu==2, P=pmatrix2d(Mod.x,Mod.z); end
    if INV.redu==3, P=pmatrix2d(Mod.x,Mod.z,S,INV.mico); end
    if isequal(numel(FIX),sim),
        unf=unique(FIX);unf=unf(unf>0); %compounds do not work
        if ~isempty(unf), P=speye(sim); end % well with combine cells
        fi=find(FIX(:)==-1);
        for i=1:length(fi), % fixed parameters
           P(:,find(P(fi(i),:)))=[];
        end
        for u=1:length(unf), % all cell compounds
            fi=find(FIX(:)==unf(u));fi1=zeros(size(fi));
            for i=1:length(fi), fi1(i)=find(P(fi(i),:)); end
            P(:,end+1)=sum(P(:,fi1),2);P(:,fi1)=[];
        end
    end
    if INV.weight==1, % Smoothness constraints of first order
        [C,Cx,Cz]=smooth2d1st(Mod.x,Mod.z,INV.rbnorm,INV.rbzfak);
        if isfield(INV,'blocky')&&(INV.blocky>0),
            if 1, % take model as structure info
                Cxz=[Cx;Cz];
                sm=abs(Cxz*dM0);su2=sum(sm.^2);sua=sum(sm);
                wxz=ones(size(sm));fi=find(sm);
                if ~isempty(fi), wxz(fi)=su2/sua./sm(fi); end
                wxz(wxz>1)=1;
                Cxz=spdiags(wxz,0,length(wxz),length(wxz))*Cxz;
                C=Cxz'*Cxz;
            else
                smx=Cx*dM0;xma=max(abs(smx));sux=sum(smx.^2);wx=1;
                smz=Cz*dM0;zma=max(abs(smz));suz=sum(smz.^2);wz=1;
                if sux>0, wx=abs(smx)*sum(abs(smx))/sux;wx(wx<1)=1; end
                if suz>0, wz=abs(smz)*sum(abs(smz))/suz;wz(wz<1)=1; end
                Cx=spdiags(1./wx,0,length(wx),length(wx))*Cx;
                Cz=spdiags(1./wz,0,length(wz),length(wz))*Cz;
                C=Cx'*Cx+Cz'*Cz;
            end
        else
            doneu=0;
            if isequal(size(XX),[length(Mod.x)-1 length(Mod.z)]-1), 
                Cx=spdiags(1-XX(:),0,size(Cx,1),size(Cx,1))*Cx;doneu=1; end
            if isequal(size(ZZ),[length(Mod.x) length(Mod.z)-1]-1), 
                Cz=spdiags(1-ZZ(:),0,size(Cz,1),size(Cz,1))*Cz;doneu=1; end
            if doneu, C=Cx'*Cx+Cz'*Cz; end
        end
        global Mstruc
        if isequal(size(Mstruc),[length(Mod.x) length(Mod.z)]-1),
                Cxz=[Cx;Cz];
                sm=abs(Cxz*log(Mstruc(:)))+1;su2=sum(sm.^2);sua=sum(sm);
                wxz=ones(size(sm));fi=find(sm);
                if ~isempty(fi), wxz(fi)=su2/sua./sm(fi); end
                wxz(wxz>1)=1;
                Cxz=spdiags(wxz,0,length(wxz),length(wxz))*Cxz;
                C=Cxz'*Cxz;
        end
    end
    if ismember(INV.weight,2:4), 
        C1=smooth2d2nd(Mod.x,Mod.z,INV.weight-2,INV.rbnorm,INV.rbzfak);
        C=C1'*C1;
    end
    if INV.weight==5,
        Cov=sum(abs(S));
        C=spdiags(1./Cov(:),0,length(Cov),length(Cov));
    end
    if INV.mitschicht, % not yet implemented
    end
    if INV.lam>0, minlam=INV.lam; end
    if (INV.auto==2)&&(INV.lam==0), INV.auto=2; end
    if INV.auto>1, % fixed regularization
        if length(dR)==size(S,2),
            dM=cglscdpt(S,dR,INV.lam,C,D,P,dM0,dM0*0,INV.maxiter,INV.silent); 
        else
            dM=cglscdp(S,dR,INV.lam,C,D,P,dM0,dM0*0,INV.maxiter,INV.silent);
        end        
        lam=INV.lam;        
    else % optlambda or choose lambda
%         [dM,lam]=invshiftcdp(S,dR,C,D,P,maxlam,minlam,kuh);
        kmax=fix(-log(maxlam/minlam)/log(kuh)+1);
        ak=maxlam*(kuh.^(0:kmax));
        global DM
        [DM,rho,eta,iter] = cglsparcdp(S,dR,C,D,P,ak);    
        save('rhoeta.mat','rho','eta','ak');
        if INV.auto==0, % manual
            global warten
            warten=0;
            choose2dmodel;
            while(warten==0), pause(1.0); end
        else %L-curve
            warten=lkurv(rho,eta,ak);
        end
        message(sprintf('Choosing lambda(%d)=%.1f',warten,ak(warten)));
        dM=DM(:,warten);
        lam=ak(warten);
    end
elseif INV.method==1, % TSVD inversion
%     [dM,lam]=invtsvd(D*S,D*dR);
    [VD,s,VM]=svd(D*S);  % Model & Data Vectors
    s=diag(s);t0=clock;r=fix(size(S,1)/2);
    message(sprintf('ready(%.2fs) min(sv)= %g, max(sv)= %g',...
        etime(clock,t0),min(s),max(s)));
    beta=((D*dR)'*VD)';
    xi=beta(1:r)./s(1:r);
    if INV.auto>1, % fixed
        dM=zeros(size(S,2),1);
        for i=1:round(INV.lam),
            dM=dM+xi(i)*VM(:,i); end
        lam=INV.lam;
    else
        global DM
        DM=zeros(size(S,2),r);
        DM(:,1)=xi(1)*VM(:,1);
        for i=2:r,
            DM(:,i)=DM(:,i-1)+xi(i)*VM(:,i); end
        eta=cumsum(xi.^2);rho=flipud(cumsum(beta(r:-1:1)));
        ak=1:r;save('rhoeta.mat','rho','eta','ak');
        if INV.auto==0, % manual
            global warten
            warten=0;
            choose2dmodel;
            while(warten==0), pause(1.0); end
            lam=warten;
        else
            lam=lkurv(rho,eta,ak);
        end
        message(sprintf('Choosing r=%d',lam));
        dM=DM(:,lam);
    end
elseif INV.method==2,
    dM=sirt(D*S,D*dR);
    lam=0;
elseif INV.method==3,
    if INV.auto>1, % fixed
        dM=regcgls(D*S,D*dR,1,round(INV.lam));
        lam=INV.lam;
    else % choosing
        global DM 
        [dm,rho,eta,j,DM]=regcgls(D*S,D*dR,0.1,1e-6,100);
        ak=1:j;save('rhoeta.mat','rho','eta','ak');
        if INV.auto==0, % manual
            global warten
            warten=0;
            choose2dmodel;
            while(warten==0), pause(1.0); end
            lam=warten;
        else
            lam=lkurv(rho,eta,ak);
        end
        message(sprintf('Choosing r=%d',lam));
        dM=DM(:,lam);
    end
end

function P = pmincov(S,mincov)

% PMINCOV - For p-matrix from sensitivity by minimum coverage
% P = pmincov(S,mincov)
% mincov - minimum coverage (default 0.4)
if nargin<2, mincov=0.4; end
P=speye(size(S,2));
COV=sum(abs(S));
P(:,COV<mincov)=[];