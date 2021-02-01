function [R,Rrez,MEA,zus]=fdfwd3d(x,y,z,M,Bg,N,FOR)

% FDFWD3D - 3D DC Forward Calculation with finite differences
% Rhoa = fdfwd3d(x,y,z,M,Bg,N,OPT)
%   x/y/z  - Model Coordinates 
%   M      - Model resistivities (1 shorter)
%   N      - Structure of electrode numbers(a,b,m,n), 
%            k-factors(k) and measurements(r)
%            elec - Electrode Positions
%   OPT    - structure of possible fields
%            method 0-electrodewise (default)
%                   1-as measured
%                   2-by sensitivity
%            acc - (relative) Accuracy for pcg (1e-4)
%            tol - Tolerance for incomplete preconditioner (1e-4)
%            maxit - Maximum Iterations for pcg (50)
%            rand - boundaring cells (4)
%            prolong - prolonging factor (4)
%            zusatz - equdistant cells outside Electrodes (2)
%            direct - equation solver (0=cg, 1=direct, -1=auto)
%
% R is geometrical mean of normal and reciprocal simulation
% [R,Rrez] = ... returns also reciprocal simulation
% [R,Rrez,MEA] = ... returns also potential matrix MEA

if nargin<6, error('At least 6 input arguments!'); end
if nargin<7,
    FOR=struct('method',0,'acc',1e-4,'tol',1e-4,'maxit',50,...
        'rand',4,'prolong',4,'zusatz',2,'direct',-1);
end
if ~isfield(FOR,'method'), FOR.method=0; end
if ~isfield(FOR,'acc'), FOR.acc=1e-4; end
if ~isfield(FOR,'tol'), FOR.tol=1e-4; end
if ~isfield(FOR,'maxit'), FOR.maxit=50; end
if ~isfield(FOR,'zusatz'), FOR.zusatz=2; end
if ~isfield(FOR,'rand'), FOR.rand=4; end
if ~isfield(FOR,'prolong'), FOR.prolong=4; end
if ~isfield(FOR,'refine'), FOR.refine=0; end
if ~isfield(FOR,'fillup'), FOR.fillup=1; end
if ~isfield(FOR,'direct'), FOR.direct=-1; end
if ~isfield(FOR,'directnodes'), FOR.directnodes=100000; end

% something wrong !!!!!
% if isempty(M), %f. compiler
%     message('Model empty, loading from model.mat!');
    M=getfield(load([tempdir 'model.mat']),'M'); 
% end
if isempty(M),
    R=ones(length(N.a),1);
    message('Model still empty, forward calculation failed!');
    return
end
zus=[0 0 0];MEA=[];R=1;Rrez=1;
Rbg=getbg(x,y,M(:,:,1),N.elec); % background resistivities
message(sprintf('Background resistivities of %.1f-%.1f',min(Rbg),max(Rbg)));
%Rbg(:)=median(Rbg); % warum klappt das nicht??!
sq=1/median(Rbg);
if Bg(1)==0, Bg=median(Rbg); end
sq=1/Bg(1);
di=FOR.refine;
if FOR.refine==0, 
    di=round(min(diff(x))/(z(2)-z(1))); 
    if di<1, di=1; end
    dj=round(min(diff(y))/(z(2)-z(1))); 
    if dj<1, dj=1; end
end
dj=di;
%refining?
if di*dj>1,
    xx=zeros(1,(length(x)-1)*di+1);
    xx(1:di:end)=x;
    for ii=2:di,
        xx(ii:di:end-1)=x(1:end-1)+diff(x)*(ii-1)/di;
    end
    x=xx;
    xx=zeros(1,(length(y)-1)*dj+1);
    xx(1:dj:end)=y;
    for jj=2:dj,
        xx(jj:dj:end-1)=y(1:end-1)+diff(y)*(jj-1)/dj;
    end
    y=xx;
end
% Blowing up coordinates
Prolong=1:FOR.zusatz;
pp=1;
for l=1:FOR.rand
    pp=pp*FOR.prolong;
    Prolong=[Prolong Prolong(end)+pp];
end
if FOR.rand>0,
    X=[min(x)-(x(2)-x(1))*flipud(Prolong(:));x(:);max(x)+(x(end)-x(end-1))*Prolong(:)]';
    Y=[min(y)-(y(2)-y(1))*flipud(Prolong(:));y(:);max(y)+(y(end)-y(end-1))*Prolong(:)]';
    Z=[z(:);max(z)+(z(end)-z(end-1))*Prolong(:)]';
    %     dzend=z(end)-z(end-1);
    %     Z=[z max(z)+(1:3)*dzend];
    %     Z=[Z max(Z)+dzend*Prolong];
else 
    X=x(:);
    Y=y(:);
    Z=z(:);
end
horz=max(X(end)-X(1),Y(end)-Y(1));
ra=FOR.rand+FOR.zusatz;
zra=ra;
while Z(end)>horz,
    Z(end)=[];
    zra=zra-1;
end
I=length(X);
J=length(Y);
K=length(Z);
SIGMA=ones(I+1,J+1,K+1)*sq;
SIGMA(:)=sq;
for l=1:length(Bg),
    SIGMA(:,:,l+1)=1/Bg(l);    
end
SIGMA(:,:,l+1:end)=1/Bg(l);
SIGMA(:,:,1)=0.0;
if ~isempty(M),
    %SIGMA(ra+2:end-ra-1,ra+2:end-ra-1,3:end-ra-1)=1./M;
    for ii=0:di-1,
        for jj=0:dj-1,
            SIGMA(ra+2+ii:di:end-ra-1+ii,ra+2+jj:dj:end-ra-1+jj,2:end-zra-1)=1./M;
        end
    end
end
if FOR.fillup,
    for rr=1:ra+1,
        SIGMA(rr,:,:)=SIGMA(ra+2,:,:);
        SIGMA(end-rr+1,:,:)=SIGMA(end-ra-1,:,:);
        SIGMA(:,rr,:)=SIGMA(:,ra+2,:);
        SIGMA(:,end-rr+1,:)=SIGMA(:,end-ra-1,:);
    end
    for rr=1:zra+1,
        SIGMA(:,:,end-rr+1)=SIGMA(:,:,end-zra-1);
    end
end
% imagesc(SIGMA(:,:,2));colorbar;return
% imagesc(squeeze(SIGMA(:,10,1:end))');colorbar;return
% minmax(SIGMA(:,:,2:end)),return
% save('fdfwd.mat','SIGMA','X','Y','Z');
I=length(X);J=length(Y);K=length(Z);
IJK=I*J*K;
if FOR.direct==-1, FOR.direct=(IJK<FOR.directnodes); end
R=zeros(length(N.a),1);
% n. Vorwaertsrechnung(Errechnung von R)
if FOR.method==2, % primitive Version (mit Sensitivity)
    global S M INV
    message('Forward modelling with sensitivity');
    if INV.lolo>0,
        R=exp(S*(log(M(:))-log(Bg(1))))*Bg(1);
    else
        R=S*(M(:)-Bg(1))+Bg(1);   
    end
else              % mit FD
    message(sprintf('Forward modelling %dx%dx%d=%d nodes',I,J,K,IJK));
    message(sprintf('Constructing matrix C for Sigma...(D&M)'));
    C=diskr(X,Y,Z,SIGMA);
    if (FOR.method==0)&&(FOR.direct==1),
        C1=diskr(X,Y,Z,(SIGMA>0));
        PM=potmap(N.elec,X,Y,Z);
        p=symamd(C);t0=clock;
        message('Starting fast direct forward calculation...');
        %save forw X Y Z N C C1 PM Rbg p
        MEA=fd3dmea(X,Y,Z,N.elec,C,C1,PM,Rbg,p);
        %save mea MEA
        message(sprintf('...ready(%.1fs)',etime(clock,t0)));
        [R,Rrez]=collectrhoa(N,MEA);
        rez=(R-Rrez)*2./(R+Rrez);
        message(sprintf('Standard deviation of reciprocity %.2f%% (max %.2f%%)',...
            std(rez)*100,max(abs(rez))*100));
        if nargout<2, R=sqrt(abs(R.*Rrez)); end
        return;
    end
    RA=[Inf Inf 0];RB=[Inf Inf 0];
    OLDRA=[0.11 0.11 0];OLDRB=OLDRA;
    
    t0=clock;
    map=symamd(C);%map=(1:IJK)';
    C=C(map,map);map=map(:);
    if FOR.direct==1,
       message(sprintf('Calculating AMD Cholesky factor...'));
       Chol=chol(C);
       message(sprintf('ready(%.2f seconds) %d nonzero elements(%dMB)',...
         etime(clock,t0),nnz(Chol),round(nnz(Chol)*8/1024/1024)));
    elseif (FOR.direct==2)&&(FOR.method==0),
        %B=zeros(IJK,length(N.a));
        B=zeros(length(N.a),IJK);
    else % iterative solution
        if isfield(FOR,'tol')&&isnumeric(FOR.tol)
            if FOR.tol==999,
                if ~isequal(size(Chol),size(C)),
                    Chol=cholinc(C,1e-4); end
            elseif FOR.tol>=2, %Jacobi
                message(sprintf('Calculating Preconditioner(Jacobi)...'));
                Chol=spdiags(sqrt(diag(C)),0,size(C,1),size(C,1)); 
            elseif FOR.tol>=1, % SSOR
                message(sprintf('Calculating Preconditioner(SSOR)...'));
                L=triu(C);D=spdiags(diag(C),0,size(C,1),size(C,1))/FOR.tol;
                Chol=inv(sqrt(D))*(D+L)/sqrt(2-FOR.tol); 
            elseif FOR.tol>0, 
                    message(sprintf('Calculating Preconditioner(Chol-Inc)...'));
                    Chol=cholinc(C,FOR.tol);
            else     
                message(sprintf('Calculating Preconditioner(Chol-Inc(SP))...'));
                Chol=cholinc(C,'0'); 
            end
        else %standard
            message(sprintf('Calculating Preconditioner(Chol-Inc(SP))...'));
            Chol=cholinc(C,'0');
        end
        message(sprintf('ready(%.2f seconds) %d nonzero elements',etime(clock,t0),nnz(Chol)));
    end % iterative solution
    Cholt=Chol';
    RR=[reshape(repmat(X,1,J*K),1,IJK);reshape(repmat(repmat(Y,1,K),I,1),1,IJK);reshape(repmat(Z,I*J,1),1,IJK)];                                         
    t0=clock;
    data=length(N.a);
    Phis=zeros(I,J,K);
    iter=0;flag=0;
    warning('off','MATLAB:divideByZero');
    if FOR.method==1,
        % ehemals forward;  % Rechnung mit "echten Dipolen"
        lauf=0;itsum=0;
        message(sprintf('Constructing matrix CQ for Sigma_q...(D&M)'));
        DSIGMA=-SIGMA+sq;DSIGMA(:,:,1)=0;
        CQ=diskr(X,Y,Z,DSIGMA);CQ=CQ(map,map);
%         C1=diskr(X,Y,Z,(SIGMA>0));
        wb=waitbar(0,'Forward calculations(as measured)...');
        aller=fix(data/25);mal=aller;
        [XX,YY,ZZ]=ndgrid(X,Y,Z);
        for l = 1:data,
            RB=[Inf Inf 0];
            RA=N.elec(N.a(l),:);
            if N.b(l)>0, RB=N.elec(N.b(l),:); end
            if (norm(RA-OLDRA)+norm(RB-OLDRB))>0,
                lauf=lauf+1;
                if RB(1)==Inf   % single source
                    %   fprintf('A=(%g %g %g)',RA(1),RA(2),RA(3));
                    Phip=calc_phip_pol(X,Y,Z,sq,RR,RA);
                else            % bipole source
                    %   fprintf('A=(%g %g %g) B=(%g %g %g)',RA(1),RA(2),RA(3),RB(1),RB(2),RB(3));                    
%                     Phip=calc_phip(X,Y,Z,sq,RR,RA,RB);
                    rada=sqrt((XX-RA(1)).^2+(YY-RA(2)).^2+(ZZ-RA(3)).^2);
                    radb=sqrt((XX-RB(1)).^2+(YY-RB(2)).^2+(ZZ-RB(3)).^2);
                    radas=sqrt((XX-RA(1)).^2+(YY-RA(2)).^2+(ZZ+RA(3)).^2);
                    radbs=sqrt((XX-RB(1)).^2+(YY-RB(2)).^2+(ZZ+RB(3)).^2);
                    Phip=(1./rada(:)+1./radas(:)-1./radb(:)-1./radbs(:))/(4*pi*sq);
                    fi=isinf(Phip);Phip(fi)=sign(Phip(fi));
                end
                b=CQ*Phip(map);
                if FOR.direct==1,
                    %Phis(map)=Chol\(b'/Chol)'; %Mist, aber wohl am besten?
                    %Phis(map)=((Chol\b)'/Chol)';
                    Phis(map)=Chol\(Cholt\b); % Speicherintensiv aber schnell
                else
                    [Phis(map),flag,err,iter,res]=pcgichol(C,b,FOR.acc,FOR.maxit,Chol,Phis(map));
                    if flag==0,
                        %fprintf('(%d) Time=%g s Iter=%g\n',l,ttt,iter);
                    else
                        message(sprintf('no Convergence (flag=%d)!',flag));
                    end
                end
                itsum=itsum+iter;
                OLDRA=RA;
                OLDRB=RB;
                mal=mal-1;
                if mal==0,
                    waitbar(l/data,wb);
                    mal=aller;
                end
            end
            RM=N.elec(N.m(l),:);
            if N.n(l)>0, 
                RN=N.elec(N.n(l),:);
                R(l)=potxy(RM(1),RM(2),sq,Phis,X,Y,RA,RB)-potxy(RN(1),RN(2),sq,Phis,X,Y,RA,RB);  
            else
                R(l)=potxy(RM(1),RM(2),sq,Phis,X,Y,RA,RB);               
            end
        end
        R=R.*N.k;
        Rrez=R;
        close(wb);
        message(sprintf('Done %d forward calculations, Time=%.2f seconds %d It.(%.1f)',...
            lauf,etime(clock,t0),itsum,itsum/lauf));
    else        
        message(sprintf('Constructing matrix C1 for Sigma=1...(D&M)'));
        C1=diskr(X,Y,Z,(SIGMA>0));C1=C1(map,map);
        anzel=size(N.elec,1);
        % Polweise Errechnung, erlaubt sens_approx und asens
        global PHI
        ia=1+FOR.rand+FOR.zusatz;ie=I-FOR.rand-FOR.zusatz;
        ja=1+FOR.rand+FOR.zusatz;je=J-FOR.rand-FOR.zusatz;
        ka=1;ke=K-zra;%ke=K-FOR.rand-FOR.zusatz;
%         PHI=zeros(length(ia:di:ie)*length(ja:dj:je)*(ke-ka+1),anzel);
%         X(ia:di:ie), Y(ja:dj:je), Z(ka:ke)
        %message(sprintf('Saving potential within x=%g-%g, y=%g-%g, z=%g-%g',...
        %   X(ia),X(ie),Y(ja),Y(je),Z(ka),Z(ke)));
        RB=[Inf Inf 0];
        flag=0;itsum=0;
        MEA=zeros(anzel);
        Phis=zeros(I,J,K);
        aller=fix(anzel/25);
        mal=aller;
        flag=0;iter=0;
        wb=waitbar(0,'Forward Calculation(every electrode)...');
        for l = 1:anzel,
            RA=N.elec(l,:);
            sq=1/Rbg(l);
            Phip=calc_phip_pol(X,Y,Z,sq,RR,RA);
            b=C1*Phip(map)*sq-C*Phip(map);
            if FOR.direct==2,
                B(l,:)=b(map)';
            else % no assembling
                if FOR.direct==1,
                    %Phis(map)=Chol\(Cholt\b);
                    Phis(map)=Chol\(b'/Chol)';
                    %Phis(map)=Cholt\(Chol\b);
		    %if l==1, save('check.mat','C','Chol','b','map'); end
                else % PCG
                    [Phis(map),flag,err,iter,res]=pcg(C,b,FOR.acc,FOR.maxit,Cholt,Chol,Phis(map));
%                    [Phis(map),flag,err,iter,res]=pcgichol(C,b,FOR.acc,FOR.maxit,Chol,Phis(map));
                    itsum=itsum+iter;
                    if flag==0
                        %      fprintf(' t=%g s It=%g E=%g\n',ttt,iter,err);
                    else
                        message(sprintf('no convergence (flag=%d)!',flag));
                    end
                end
                %Phigesamt=Phip+Phis;    % Phis nur ausgeborgt
                %             iphis=Phis(ia:di:ie,ja:dj:je,ka:ke);
                % !!!! noch besser alle, und dann besser integrieren...
                %             PHI(:,l)=iphis(:);
                %%%PHI(:,l)=reshape(,(ie-ia+1)*(je-ja+1)*(ke-ka+1),1);  % Nur für Sens. wichtigen Teil
                %mea=potxy(N.elec(:,1),N.elec(:,2),sq,Phis,X,Y,RA);
                %             phisi=Phis(ia:di:ie,ja:dj:je,ka:ke);PHI(:,l)=phisi(:);
                el=N.elec;el(:,1)=el(:,1)-RA(1);el(:,2)=el(:,2)-RA(2);%el(l,1)=1;
                su=sqrt(sum(el.^2,2));su(l)=1;
                mea=1./su/(2*pi*sq); % Phip, nur für Oberfläche
                mea=mea+interp2(Y,X,Phis(:,:,1),N.elec(:,2),N.elec(:,1));
                mea(l)=0;
                if find(0>mea),
                    message(sprintf('Found %d negative potentials at source %d (%.1f,%.1f)',...
                        length(find(0>mea)),l,N.elec(l,1),N.elec(l,2)));
                end
                MEA(l,:)=mea';  % Potentiale Nr. i bei Einspeisung in j
                mal=mal-1;
                if mal==0,
                    waitbar(l/anzel,wb);
                    mal=aller;
                end
            end
        end
        close(wb)
        if FOR.direct==2,
            %save('numtest.mat','C','B');
            message(sprintf('Calculating AMD Cholesky factor...'));t0=clock;
            Chol=chol(C);
            message(sprintf('ready(%.2f seconds) %d nonzero elements',etime(clock,t0),nnz(Chol)));
            tic;PHIS(map,:)=Chol\(B/Chol)';toc;
            PM=potmap(N.elec,X,Y,Z); % hier fehlt noch was !!!
        end
        message(sprintf('Done %d forward calculations, Time=%.2f sec %d It.(%.1f/El.)',...
            size(N.elec,1),etime(clock,t0),itsum,itsum/anzel));
        [R,Rrez]=collectrhoa(N,MEA);
        rez=(R-Rrez)*2./(R+Rrez);
        message(sprintf('Standard deviation of reciprocity %.2f%% (max %.2f%%)',...
            std(rez)*100,max(abs(rez))*100));
    end % pole-wise
    if nargout<2,
        R=sqrt(abs(R.*Rrez));
    end
end % not by sens
l=length(find(R<0));
if l>0,
    message(sprintf('Attention! %d values of R were negative!',l));
    R=abs(R);
end
l=length(find(R==0));
if l>0, message(sprintf('Attention! %d values of R are zero!',l)); end
