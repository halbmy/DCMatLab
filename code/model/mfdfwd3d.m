function [R,Rrez,MEA]=mfdfwd3d(Mod,N,FOR)

% FDFWD3D - 3D DC Forward Calculation with finite differences
% Rhoa = mfdfwd3d(Model,N,OPT)
% Model    - model structure
%   N      - Structure of electrode numbers(a,b,m,n), 
%            k-factors(k) and measurements(r)
%            elec- Electrode Positions
%   OPT    - structure of possible fields
%            method 0-electrodewise (default)
%                   1-as measured
%                   2-by sensitivity
%            acc - Accuracy of forward step
%            tol - Tolerance for incomplete preconditioner
%            maxit - Maximum Iterations for pcg
%            rand - boundaring cells (3)
%            prolong - prolonging factor (5)
%            zusatz - cells outside Electrodes (2)
%            direct - equation solver (0=cg, 1=direct, -1=auto)

if nargin<2, error('At least 2 input arguments!'); end
if nargin<3,
    FOR=struct('method',0,'acc',1e-3,'tol',1e-4,'maxit',50,...
        'rand',4,'prolong',5,'zusatz',4,'direct',-1); 
end
if ~iscell(Mod.M), % grid model
    if nargout>1,
        if isfield(Mod,'Bg'),
            [R,Rrez,MEA]=fdfwd3d(Mod.x,Mod.y,Mod.z,Mod.M,Mod.Bg,N,FOR);
        else
            [R,Rrez,MEA]=fdfwd3d(Mod.x,Mod.y,Mod.z,Mod.M,0,N,FOR);
        end
    else
        M=Mod.M;
        save([tempdir 'model.mat'],'M');
        if isfield(Mod,'Bg'),
            R=fdfwd3d(Mod.x,Mod.y,Mod.z,Mod.M,Mod.Bg,N,FOR);
        else
            R=fdfwd3d(Mod.x,Mod.y,Mod.z,Mod.M,0,N,FOR);
        end
        %         delete([tempdir 'model.mat']);
    end
    return
else % little hack!
    [M,x,y,z]=mesch3dmodel(Mod);
    R=fdfwd3d(x,y,z,M,Mod.Bg,N,FOR);
    return;
end
if ~isfield(FOR,'method'), FOR.method=0; end
if ~isfield(FOR,'acc'), FOR.acc=1e-3; end
if ~isfield(FOR,'tol'), FOR.tol=1e-4; end
if ~isfield(FOR,'maxit'), FOR.maxit=50; end
if ~isfield(FOR,'zusatz'), FOR.zusatz=1; end
if ~isfield(FOR,'rand'), FOR.rand=3; end
if ~isfield(FOR,'prolong'), FOR.prolong=5; end
if ~isfield(FOR,'refine'), FOR.refine=1; end
if ~isfield(FOR,'direct'), FOR.direct=0; end

nz=length(Mod.z)-1;
nx=size(Mod.M{1},1);ny=size(Mod.M{1},2);
nnx=Mod.nx(1)*nx;nny=Mod.ny(1)*ny;
z=Mod.z;Bg=Mod.Bg;
sq=1/Mod.Bg(1);

x=(0:nnx)*Mod.dx+Mod.x0;
y=(0:nny)*Mod.dy+Mod.y0;
Rbg=getbg(x,y,Mod.M{1},N.elec);

% Blowing up coordinates
Prolong=1:FOR.zusatz;
pp=1;
for l=1:FOR.rand
    pp=pp*FOR.prolong;
    Prolong=[Prolong Prolong(end)+pp];
end
if FOR.rand>0,
    X=[min(x)-(x(2)-x(1))*fliplr(Prolong) x max(x)+(x(end)-x(end-1))*Prolong];
    Y=[min(y)-(y(2)-y(1))*fliplr(Prolong) y max(y)+(y(end)-y(end-1))*Prolong];
    Z=[z max(z)+(z(end)-z(end-1))*Prolong];
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
SIGMA(:,:,1)=0.0;
for l=1:length(Bg),
    SIGMA(:,:,l+1)=1/Bg(l);    
end
SIGMA(:,:,l+1:end)=1/Bg(l);
ra=FOR.rand+FOR.zusatz+1;
for k=1:nz,
    nx=Mod.nx(k);ny=Mod.ny(k);
    modk=Mod.M{k};
    sx=size(modk,1);sy=size(modk,2);
    rx=floor(mod(nnx,sx*Mod.nx(k))/2);
    ry=floor(mod(nny,sy*Mod.ny(k))/2);
    for i=1:nx,
        for j=1:ny,
            SIGMA(ra+i+rx:nx:end-ra-rx+i-nx,ra+j+ry:ny:end-ra-ry+j-ny,k+1)=1./modk;
        end
    end
end
IJK=I*J*K;
if FOR.direct==-1, FOR.direct=(IJK<100000); end
R=zeros(length(N.r),1);
% n. Vorwaertsrechnung(Errechnung von R)
if FOR.method==2, % primitive Version (mit Sensitivity)
    global S INV
    message('Forward modelling with sensitivity');
    M=[];
    for k=1:length(Mod.M), M=[M;Mod.M{k}(:)]; end
    if INV.lolo>0,
        R=exp(S*(log(M(:))-log(Mod.Bg(1))))*Mod.Bg(1);
    else
        R=S*(M(:)-Bg(1))+Bg(1);   
    end
else              % mit FD
    message(sprintf('Forward modelling %dx%dx%d=%d nodes',I,J,K,IJK));
    message(sprintf('Constructing matrix CQ for Sigma_q...(D&M)'));
    CQ=diskr_dm(X,Y,Z,-SIGMA+sq);
    message(sprintf('Constructing matrix C for Sigma...(D&M)'));
    C=diskr_dm(X,Y,Z,SIGMA);
    if FOR.direct==1,
        C1=diskr(X,Y,Z,(SIGMA>0));
        PM=potmap(N.elec,X,Y,Z);
        p=symamd(C);t0=clock;
        message('Starting fast direct forward calculation...');
        MEA=fd3dmea(X,Y,Z,N.elec,C,C1,PM,Rbg,p);
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
    map=symamd(C);%map=(1:IJK)';
    C=C(map,map);map=map(:);
    t0=clock;
    if FOR.direct==1, %direct solver via AMD Cholesky factorization
        message('Calculating AMD Cholesky Factor...');
        Chol=chol(C);
        message(sprintf('ready(%.2f seconds) %d nonzero elements',etime(clock,t0),nnz(Chol)));
    else,
        if isfield(FOR,'tol')&&isnumeric(FOR.tol)&&(FOR.tol>0),
            if FOR.tol>=2, % Jacobi
                message('Calculating Preconditioner(Jacobi)...');
                Chol=spdiags(sqrt(diag(C)),0,size(C,1),size(C,1));
            elseif FOR.tol>=1, % SSOR
                message('Calculating Preconditioner(SSOR)...');
                L=triu(C);D=spdiags(diag(C),0,size(C,1),size(C,1));
                Chol=inv(sqrt(D))*(D+L)/sqrt(2-FOR.tol);
            elseif FOR.tol>0,
                message('Calculating Preconditioner(Chol-Inc)...');
                Chol=cholinc(C,FOR.tol);
            else
                message(sprintf('Calculating Preconditioner(Chol-Inc-SP)...'));
                Chol=cholinc(C,'0');
            end
            Cholt=Chol';
            message(sprintf('ready(%.2f seconds) %d nonzero elements',etime(clock,t0),nnz(Chol)));
        else
            message('Calculating preconditioner(Chol-Inc-SP)...');
            Chol=cholinc(C,'0');
        end
    end
    RR=[reshape(repmat(X,1,J*K),1,IJK);reshape(repmat(repmat(Y,1,K),I,1),1,IJK);reshape(repmat(Z,I*J,1),1,IJK)];                                         
    t0=clock;
    data=length(N.r);
    Phis=zeros(I,J,K);
    if FOR.method==1,
        % ehemals forward;  % Rechnung mit "echten Dipolen"
        lauf=0;itsum=0;
        wb=waitbar(0,'Forward calculations(as measured)...');
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
                    Phip=calc_phip(X,Y,Z,sq,RR,RA,RB);
                end
                b=CQ*Phip(map);
                if FOR.direct==1,
                    Phis(map)=Chol\(b'/Chol)';
                else
                    [Phis(map),flag,err,iter,res]=pcg(C,b,FOR.acc,FOR.maxit,Cholt,Chol,Phis(map));
                end
                if flag~=0,
                    message(sprintf('no Convergence (flag=%d)!',flag));
                end
                itsum=itsum+iter;
                OLDRA=RA;
                OLDRB=RB;
                waitbar(l/data,wb);
            end
            RM=N.elec(N.m(l),:);
            if N.n(l)>0, 
                RN=N.elec(N.n(l),:);
                R(l)=potxy(RM(1),RM(2),sq,Phis,X,Y,RA,RB)-potxy(RN(1),RN(2),sq,Phis,X,Y,RA,RB);  
            else
                R(l)=potxy(RM(1),RM(2),sq,Phis,X,Y,RA,RB);               
            end
        end
        close(wb)
        message(sprintf('Done %d forward calculations, Time=%.2f seconds %d It.(%.1f)',...
            lauf,etime(clock,t0),itsum,itsum/lauf));
    else        
        % Polweise Errechnung, erlaubt sens_approx und asens
        %global PHI
        ia=1+FOR.rand+FOR.zusatz;ie=I-FOR.rand-FOR.zusatz;
        ja=1+FOR.rand+FOR.zusatz;je=J-FOR.rand-FOR.zusatz;
        ka=2;ke=K-FOR.rand-FOR.zusatz;
        di=Mod.nx(1);dj=Mod.ny(1);
        %[X(ia) X(ie) Y(ja) Y(je) Z(ka) Z(ke)]
        anzel=size(N.elec,1);
        RB=[Inf Inf 0];
        flag=0;itsum=0;
        if ~FOR.direct,
            wb=waitbar(0,'Forward Calculation(every electrode)...');
        end
        MEA=zeros(anzel);
        PHI=zeros(length(ia:di:ie)*length(ja:dj:je)*(ke-ka+1),anzel);
        Phis=zeros(I,J,K);
        aller=fix(anzel/25);
        mal=aller;
        flag=0;iter=0;
        if FOR.direct==2, B=zeros(anzel,IJK); end
        C1=diskr_dm(X,Y,Z,SIGMA>0);C1=C1(map,map);
        for l = 1:anzel,
            RA=N.elec(l,:);
            Phip=calc_phip_pol(X,Y,Z,sq,RR,RA);
            b=C1*Phip(map)*sq-C*Phip(map); %(sq*C1-C)*Phip(:);
            if FOR.direct==2,
                B(:,l)=b(map)';
            else
                if FOR.direct==1, % AMD Cholesky
                    Phis(map)=Chol\(b'/Chol)';
                else % PCG
                    [Phis(map),flag,err,iter,res]=pcg(C,CQ*Phip(:),...
                        FOR.acc,FOR.maxit,Cholt,Chol,Phis(map));
                    itsum=itsum+iter;
                    if flag==0
                        %      fprintf(' t=%g s It=%g E=%g\n',ttt,iter,err);
                    else
                        message(sprintf('no convergence (flag=%d)!',flag));
                    end
                end
                %Phigesamt=Phip+Phis;    % Phis nur ausgeborgt
                iphis=Phis(ia:di:ie,ja:dj:je,ka:ke);
                % !!!! noch besser alle, und dann besser integrieren...
                %PHI(:,l)=iphis(:);
                %PHI(:,l)=reshape(Phis(ia:ie,ja:je,ka:ke),(ie-ia+1)*(je-ja+1)*(ke-ka+1),1);  % Nur für Sens. wichtigen Teil
                %mea=potxy(N.elec(:,1),N.elec(:,2),sq,Phis,X,Y,RA);
                el=N.elec;el(:,1)=el(:,1)-RA(1);el(:,2)=el(:,2)-RA(2);el(l,1)=1;
                su=sqrt(sum(el.^2,2));su(l)=1;
                mea=1./su/(2*pi*sq); % Phip, nur für Oberfläche
                mea=mea+interp2(Y,X,Phis(:,:,1),N.elec(:,2),N.elec(:,1));
                mea(l)=0;
                if find(0>mea),
                    [ l N.elec(l,:) length(find(0>mea))]
                    pause(0.1);
                end
                MEA(l,:)=mea';  % Potentiale Nr. i bei Einspeisung in j
                mal=mal-1;
                if mal==0,
                    waitbar(l/anzel,wb);
                    mal=aller;
                end
            end % direct
        end
        if FOR.direct==2,
            %save('numtest.mat','C','B');
            %tic;PHIS=C\B;toc;
            message('Calculating AMD Cholesky factor...');
            Chol=Chol(C);
            PHIS(map,:)=Chol\(B/Chol)';
            PM=potmap(N.elec,X,Y,Z)*PHIS;
        else
            close(wb);
        end
        wb=waitbar(0,'Reordering...');
        aller=fix(data/25);
        mal=aller;
        R=zeros(size(N.r));Rrez=R;
        for l = 1:data,
            R(l)=MEA(N.a(l),N.m(l));
            Rrez(l)=MEA(N.m(l),N.a(l));
            if N.n(l)>0, 
                R(l)=R(l)-MEA(N.a(l),N.n(l)); 
                Rrez(l)=Rrez(l)-MEA(N.n(l),N.a(l)); 
            end
            if N.b(l)>0,
                R(l)=R(l)-MEA(N.b(l),N.m(l));
                Rrez(l)=Rrez(l)-MEA(N.m(l),N.b(l));
                if N.n(l)>0, 
                    R(l)=R(l)+MEA(N.b(l),N.n(l)); 
                    Rrez(l)=Rrez(l)+MEA(N.n(l),N.b(l)); 
                end
            end
            mal=mal-1;
            if mal==0,
                mal=aller;
                waitbar(l/data,wb);
            end
        end
        close(wb)
        message(sprintf('Done %d forward calculations, Time=%.2f sec %d It.(%.1f/El.)',...
            size(N.elec,1),etime(clock,t0),itsum,itsum/anzel));
    end
    R(:)=R(:).*N.k(:);
    Rrez(:)=Rrez(:).*N.k(:);
    if FOR.direct,
        R=R+1/sq;
        Rrez=Rrez+1/sq;
    end
    rez=(R-Rrez)*2./(R+Rrez);
    message(sprintf('Standard deviation of reciprocity %.2f%%',std(rez)*100));
    if nargout<2,
        R=sqrt(abs(R.*Rrez));
    end
end
