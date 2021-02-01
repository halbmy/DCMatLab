nM=M;running=1;dM0=0;
lbound=0;ubound=0;
if isfield(INV,'lbound'), lbound=INV.lbound; end
if isfield(INV,'ubound'), ubound=INV.ubound; end
if lbound>min(min(N.r),min(R)), lbound=0; end
if (ubound>0)&&(ubound<max(max(N.r),max(R))), ubound=0; end
RMS=rms(N.r,R,INV.lolo);
CHIQ=chi2(N.r,R,N.err,INV.lolo);
message(['CHI^2 = ',sprintf('%.2f ',CHIQ) '(RMS=' sprintf('%.2f%%',RMS(end)) ')']);
while running,
    % Inverse subproblem with log(Data) and log(Model)
    if lbound>min(min(R),min(N.r)), lbound=0;display('resetting lower bound'); end
    if (ubound>0)&&(ubound<max(max(R),max(N.r))),
        ubound=0;display('resetting upper bound'); end
    if INV.lolo,
        dR=log(N.r(:)-lbound)-log(R(:)-lbound);
        if ubound>0, dR=dR-log(ubound-N.r(:))+log(ubound-R(:)); end
        if isequal(size(M),size(Mref)),
            dM0=log(M(:)-lbound)-log(Mref(:)-lbound); 
            if ubound>0, dM0=dM0-log(ubound-M(:))+log(ubound-Mref(:)); end
        end
    else
        dR=N.r(:)-R(:);
        if isequal(size(M),size(Mref)), dM0=M(:)-Mref(:); end
    end
    fi=(INV.lam==0);
    [dM,lam]=invers(S,dR,INV,N.err,dM0);
    if (INV.auto==1)&&(INV.const==1), 
        INV.auto=2;INV.lam=lam; 
    end
    if fi, INV.first=lam; end
    % Model updating
    nM=M;fak=1; % line search factor
    if isfield(INV,'linesearch')&&(INV.linesearch>0),
        if INV.lolo, 
            if ubound>0,
                nM(:)=(ubound*(M(:)-lbound).*exp(dM)+lbound*(ubound-M(:)))./((M(:)-lbound).*exp(dM)+ubound-M(:));
            else
                nM(:)=(M(:)-lbound).*exp(dM)+lbound;
            end
        else
            nM(:)=M(:)+dM; 
        end
        R1=abs(dcfwd2d(x,z,nM,Lay,N,FOR));
        [fak,appR]=linesearch(N,R,R1,INV.lolo);
        if fak>0.9, fak=1; end
        if fak<0, fak=0.1; end
        if fak<1, dM=dM*fak;
            message(sprintf('Line search factor = %.2f chi^2=%.1f-->%.1f',...
                fak,chi2(N.r,R1,N.err,INV.lolo),chi2(N.r,appR,N.err,INV.lolo))); end
    else
        fak=1.001;
    end
    if INV.lolo,
        if ubound>0,
           nM(:)=(ubound*(M(:)-lbound).*exp(dM)+lbound*(ubound-M(:)))./((M(:)-lbound).*exp(dM)+ubound-M(:));
        else
            nM(:)=(M(:)-lbound).*exp(dM)+lbound;
        end
    else
        nM(:)=M(:)+dM;
    end
    patch2dmodel(x,z,nM,MAL,N);
    % Forward Calculation
    oldR=R;
    if fak==1, R=R1; elseif fak==0, R=oldR;
    else R=abs(dcfwd2d(x,z,nM,Lay,N,FOR)); end
    % stop if rms grows or changes less than 5 percent
    RMS=[RMS rms(N.r,R,INV.lolo)];
    CHIQ=[CHIQ chi2(N.r,R,N.err,INV.lolo)];
    message(['CHI^2 = ',sprintf('%.2f ',CHIQ) '(RMS=' sprintf('%.2f%%',RMS(end)) ')']);
    dchiq=CHIQ(end)-CHIQ(end-1);
    if (dchiq<0)||(length(RMS)<3),
        M=nM;
        running=(abs(dchiq)/CHIQ(end)>0.05);
        if running==0, message('Stopping (dchi^2<5 percent)'); end
    else
        message('Going back to old model and stopping');
        running=0;
        R=oldR;
        patch2dmodel(x,z,M,MAL,N);
        CHIQ(end)=[];RMS(end)=[];
    end
    if CHIQ(end)<1, running=0; end
    if running,
        % S = S + ((log(R)-log(oldR)-S*dM)*(dM'))/(dM'*dM); % BROYDEN
        if (size(S,1)==length(R))&&(size(S,2)==length(dM)),
            ddr=(log(R)-log(oldR)-S*dM(:))/(dM(:)'*dM(:));
            for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
            message('Doing Broyden update of sensitivity matrix!');
        end
        if (size(S,2)==length(R))&&(size(S,1)==length(dM)), %trans+sparse
            ddr=(log(R)-log(oldR)-(dM(:)'*S)')'/(dM(:)'*dM(:));
            for i=1:size(S,1), 
                fi=find(S(i,:));
                S(i,fi)=S(i,fi)+ddr(fi)*dM(i);
            end
            message('Doing Broyden update of sensitivity matrix!');
        end
    end
    if isfield(INV,'robust')&&(INV.robust>0),
        dR=(log(N.r)-log(R))./log(1+N.err);
        w=abs(dR)*sum(abs(dR))/sum(dR.^2);
        w(w<1)=1;
        N.err=N.err.*w;
        message('Iteratively data reweighting (robust inversion)');
    end
end
