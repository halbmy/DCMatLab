% function Model=invauto3d(N,Model,INV,FOR)
% global S

good=1;
if ~isfield(Model,'R'), Model.R=mfdfwd3d(Model,N,FOR); end
RMS=rms(N.r,Model.R,INV.lolo);
CHIQ=chi2(N.r,Model.R,N.err,INV.lolo);
message(['CHI^2 = ',sprintf('%.2f ',CHIQ) '(RMS=' sprintf('%.2f%%',RMS(end)) ')']);

while good,
    oldMod=Model;oldR=Model.R;
    if ~iscell(Model.M),
        Model.ncells=prod(size(Model.M));
    else
        nn=0;for k=1:length(Model.M), nn=nn+prod(size(Model.M{k})); end
        Model.ncells=nn;
    end

    lbound=0;ubound=0;
    if isfield(INV,'lbound'), lbound=INV.lbound; end
    if isfield(INV,'ubound'), ubound=INV.ubound; end
    if lbound>min(min(N.r),min(Model.R)), lbound=0; end
    if (ubound>0)&&(ubound<max(max(N.r),max(Model.R))), ubound=0; end

    %% Calculating Residuum
    if INV.lolo,
        dR=log(N.r(:)-lbound)-log(Model.R(:)-lbound);
        if ubound>0, dR=dR-log(ubound-N.r(:))+log(ubound-Model.R(:)); end
    else
        dR=N.r(:)-Model.R(:);
    end
    % inverse subproblem
    dBg=[]; % to fix!!!
    if iscell(Model.M),
        MM=[];for k=1:length(Model.M), MM=[MM;Model.M{k}(:)]; end
    else MM=Model.M(:); end
    if INV.lolo,
        dM0=log(MM-lbound);
        if ubound>0, dM0=dM0-log(ubound-Model.M(:)); end
        if isequal(size(MM),size(Mref)),
            dM0=dM0-log(Mref-lbound);
            if ubound>0, dM0=dM0+log(ubound-Mref(:)); end
        end
    else
        if isequal(numel(MM),numel(Mref)), dM0=MM(:)-Mref(:); end
    end
    set(gcf,'Pointer','watch');
    [dM,lam]=minvers(S,dR,INV,N.err,Model,dM0);
    tauopt=1;
    oldR=Model.R;
    if isfield(INV,'linesearch')&&(INV.linesearch>0),
        set(gcf,'Pointer','watch');
        newModel=modelupdate(Model,dM,2-INV.lolo,lbound,ubound);
        %if ~iscell(newModel.M), minmax(newModel.M), end
        [cmin,cmax]=draw3dmodel(newModel,MAL,[]);drawnow;
%         malcbar(handles,cmin,cmax,MAL.clog);
        Model.R=mfdfwd3d(newModel,N,FOR);
        [tauopt,appR]=linesearch(N,oldR,Model.R,INV.lolo);
	newchiq=chi2(N.r,Model.R,N.err,1);
        if (tauopt==0)||(newchiq>CHIQ(end)*0.95), % line search failed
            taug=0.3; % sampling point
            newModel=modelupdate(Model,dM*taug,2-INV.lolo,lbound,ubound);
            Ri=mfdfwd3d(newModel,N,FOR);
            chiq=[CHIQ(end) chi2(N.r,Ri,N.err,1) newchiq];
            taus=[0 taug 1];
            G=[1 0 0;1 taug taug^2;1 1 1];
            xx=G\chiq';
            tauopt=-xx(2)/2/xx(3);
            if tauopt<=0,    
                message('(parabolic) line search failed, stopping iteration');
                dc3dinvres('showmod_Callback',gcbo,[],guidata(gcbo));
                set(gcf,'Pointer','arrow');
                RMS(end+1)=RMS(end);
                CHIQ(end+1)=CHIQ(end);
                dc3dinvres('rmsout_Callback',gcbo,[],guidata(gcbo))
                return;
            else
                if ~isfield(INV,'lolo')||(INV.lolo>0), 
                    appR=oldR.*exp(tauopt*(log(R)-log(oldR))); 
                else appR=oldR+tauopt*(R-oldR); end
                message(sprintf('Parabolic line search parameter %.2f, chi^2=%.1f-%.1f',...
                    tauopt,newchiq,chi2(N.r,appR,N.err,1)));
            end
            if abs(tauopt-taug)<=0.5, % close to sample
                tauopt=taug; % take sampling point and save time
                dM=dM*tauopt;
                Model.R=Ri; 
                tauopt=1; % avoid new calculation
            end
        else
            if tauopt<1, message(sprintf('Line search parameter %.2f, chi^2=%.1f-%.1f',...
                    tauopt,newchiq,chi2(N.r,appR,N.err,1))); end
            dM=dM*tauopt;
        end        
        if tauopt<1, message(sprintf('Line search parameter %.2f, chi^2=%.1f-%.1f',...
                tauopt,newchiq,chi2(N.r,appR,N.err,1))); end
        dM=dM*tauopt;
    end
    set(gcf,'Pointer','arrow');
    % Model update
    if INV.lolo==1, % logarithmic model parameters
        %     ma=find(dM==max(dM));ma=ma(1);mi=find(dM==min(dM));mi=mi(1);
        %     message(sprintf('Max(dM) = %g Min(dM) = %g',exp(dM(ma)),exp(dM(mi))));
        Model=modelupdate(Model,dM,1,lbound,ubound);
        if length(dBg)==length(Model.Bg), Model.Bg(:)=Model.Bg(:).*exp(dBg(:)); end
    else
        %     message(sprintf('Max(dM) = %g Min(dM) = %g',max(dM),min(dM)));
        Model=modelupdate(Model,dM,2);
        if length(dBg)==length(Model.Bg), Model.Bg(:)=Model.Bg(:)+dBg(:); end
    end
    if iscell(Model.M),
        MM=[];for k=1:length(Model.M), MM=[MM;Model.M{k}(:)]; end
        mi=min(MM);ma=max(MM);
    else
        mi=min(Model.M(:));ma=max(Model.M(:));
    end
    message(sprintf('Min(Model) = %.1f Ohmm, Max(Model) = %.1f Ohmm',mi,ma));
    %if ~iscell(Model.M), minmax(Model.M), end
    draw3dmodel(Model,MAL);
    drawnow;
    if FOR.method==2,
        rho=Model.Bg(1);
        if iscell(Model.M), % Para Model
            MM=[];for k=1:length(Model.M), MM=[MM;Model.M{k}(:)]; end
        else MM=Model.M(:); end
        Model.R=rho.*exp(S*(log(MM)-log(rho)));
    else
        set(gcf,'Pointer','watch');
        if ~isfield(INV,'linesearch')||(INV.linesearch==0)||(tauopt<0.95), Model.R=mfdfwd3d(Model,N,FOR); end
    end
    message(sprintf('Min(Model) = %.1f Ohmm, Max(Model) = %.1f Ohmm',mi,ma));
    RMS=[RMS rms(N.r,Model.R,INV.lolo)];
    CHIQ=[CHIQ chi2(N.r,Model.R,N.err,INV.lolo)];
    message(['CHI^2 = ',sprintf('%.2f ',CHIQ) '(RMS=' sprintf('%.2f%%',RMS(end)) ')']);
    if isfield(INV,'robust')&&(INV.robust>0),
        if ~isfield(INV,'lolo')||(INV.lolo==0),
            dR=(N.r-Model.R)./N.err;
        else
            dR=(log(N.r)-log(Model.R))./log(1+N.err);
        end
        w=abs(dR)*sum(abs(dR))/sum(dR.^2);
        w(w<1)=1;
        N.err=N.err.*w;
        message('Iteratively data reweighting (robust inversion)');
    end
    if ~isfield(INV,'broyden')||(INV.broyden>0),
        message('Sensitivity updating by broyden method');
        if size(S,1)==numel(dM),
            if isfield(INV,'lolo')&&(INV.lolo==0),
                ddr=(Model.R-oldR-(dM(:)*S)')'/(dM(:)'*dM(:));
            else
                ddr=(log(Model.R)-log(oldR)-(dM(:)'*S)')'/(dM(:)'*dM(:));
            end
            if issparse(S),
                for i=1:size(S,1),
                    fi=find(S(i,:));
                    S(i,fi)=S(i,fi)+ddr(fi)*dM(i);
                end
            else
                for i=1:size(S,2), S(i,:)=S(i,:)+ddr*dM(i); end
            end
        else
            if isfield(INV,'lolo')&&(INV.lolo==0),
                ddr=(Model.R-oldR-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
            else
                ddr=(log(Model.R)-log(oldR)-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
            end
            if issparse(S),
                for i=1:size(S,2),
                    fi=find(S(:,i));
                    S(fi,i)=S(fi,i)+ddr(fi)*dM(i);
                end
            else
                for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
            end
        end
    end % broyden
    dchiq=1-CHIQ(end)/CHIQ(end-1);
    if dchiq<0, % new model not better
        good=0;
        message('Going back to old model (increasing CHI^2) and stop');
        Model=oldMod;Model.R=oldR;
    end
    if CHIQ(end)<1, message('CHI^2<1, stopping inversion'); end
    good=good&(CHIQ(end)>1)&(dchiq>0.05);    
end
