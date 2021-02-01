function [R,Rrez,MEA]=fdfwd3d(x,y,z,M,Bg,N,FOR)

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
%            acc - Accuracy of forward step
%            tol - Tolerance for incomplete preconditioner
%            maxit - Maximum Iterations for pcg
%            rand - boundaring cells (3)
%            prolong - prolonging factor (5)
%            zusatz - cells outside Electrodes (2)
% [R,Rrez,MEA] = ... gives also reciproce measurements and MEA matrix

if nargin<6, error('At least 6 input arguments!'); end
if nargin<7,
    FOR=struct('method',0,'acc',1e-3,'tol',1e-4,'maxit',50,...
        'rand',3,'prolong',5,'zusatz',1);
end
if ~isfield(FOR,'method'), FOR.method=0; end
if ~isfield(FOR,'acc'), FOR.acc=1e-3; end
if ~isfield(FOR,'tol'), FOR.tol=1e-4; end
if ~isfield(FOR,'maxit'), FOR.maxit=50; end
if ~isfield(FOR,'zusatz'), FOR.zusatz=1; end
if ~isfield(FOR,'rand'), FOR.rand=3; end
if ~isfield(FOR,'prolong'), FOR.prolong=5; end
if ~isfield(FOR,'refine'), FOR.refine=0; end

% something wrong !!!!!
if isempty(M), load('model.mat'); end
if isempty(M),
    R=ones(length(N.r,1),1);
    message('Model empty, forward calculation failed!');
    return
end
%end

Rbg=getbg(x,y,M(:,:,1),N.elec); % background resistivities
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
    X=[min(x)-(x(2)-x(1))*fliplr(Prolong) x max(x)+(x(end)-x(end-1))*Prolong];
    Y=[min(y)-(y(2)-y(1))*fliplr(Prolong) y max(y)+(y(end)-y(end-1))*Prolong];
    Z=[z max(z)+(z(end)-z(end-1))*Prolong];
    %     dzend=z(end)-z(end-1);
    %     Z=[z max(z)+(1:3)*dzend];
    %     Z=[Z max(Z)+dzend*Prolong];
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
if ~isempty(M),
    %SIGMA(ra+2:end-ra-1,ra+2:end-ra-1,3:end-ra-1)=1./M;
    for ii=0:di-1,
        for jj=0:dj-1,
            SIGMA(ra+2+ii:di:end-ra-1+ii,ra+2+jj:dj:end-ra-1+jj,2:end-zra-1)=1./M;
        end
    end
end
% save('fdfwdo.mat','SIGMA','X','Y','Z');
I=length(X);J=length(Y);K=length(Z);
IJK=I*J*K;
R=zeros(length(N.r),1);
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
    RA=[Inf Inf 0];RB=[Inf Inf 0];
    OLDRA=[0.11 0.11 0];OLDRB=OLDRA;
    
    message(sprintf('Calculating Preconditioner(Chol-Inc)...'));
    t0=clock;
    %[L,U]=luinc(C,FOR.tol);
    %     Chol=ichol(C);
    if isfield(FOR,'tol')&&isnumeric(FOR.tol)&&(FOR.tol>0),
        Chol=cholinc(C,FOR.tol);
    else
        Chol=cholinc(C,'0');
    end
    Cholt=Chol';
    message(sprintf('ready(%.2f seconds)',etime(clock,t0)));
    
    RR=[reshape(repmat(X,1,J*K),1,IJK);reshape(repmat(repmat(Y,1,K),I,1),1,IJK);reshape(repmat(Z,I*J,1),1,IJK)];                                         
    t0=clock;
    data=length(N.r);
    Phis=zeros(I,J,K);
    if FOR.method==1,
        % ehemals forward;  % Rechnung mit "echten Dipolen"
        lauf=0;itsum=0;
        message(sprintf('Constructing matrix CQ for Sigma_q...(D&M)'));
        CQ=diskr(X,Y,Z,-SIGMA+sq);
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
                q=CQ*Phip(:);
                %tic;[Phis,flag,err,iter,res]=pcg(C,q,FOR.acc,FOR.maxit,L,U,Phis(:));ttt=toc;
                tic;[Phis,flag,err,iter,res]=pcg(C,q,FOR.acc,FOR.maxit,Cholt,Chol,Phis(:));ttt=toc;
                Phis=reshape(Phis,I,J,K);
                if flag==0,
                    %fprintf('(%d) Time=%g s Iter=%g\n',l,ttt,iter);
                else
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
        message(sprintf('Constructing matrix C1 for Sigma=1...(D&M)'));
        C1=diskr(X,Y,Z,(SIGMA>0));
        anzel=size(N.elec,1);
        % Polweise Errechnung, erlaubt sens_approx und asens
        global PHI
        ia=1+FOR.rand+FOR.zusatz;ie=I-FOR.rand-FOR.zusatz;
        ja=1+FOR.rand+FOR.zusatz;je=J-FOR.rand-FOR.zusatz;
        ka=1;ke=K-FOR.rand-FOR.zusatz;
        PHI=zeros(length(ia:di:ie)*length(ja:dj:je)*(ke-ka+1),anzel);
%         X(ia:di:ie), Y(ja:dj:je), Z(ka:ke)
        %[X(ia) X(ie) Y(ja) Y(je) Z(ka) Z(ke)]
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
            [Phis(:),flag,err,iter,res]=pcg(C,C1*Phip(:)*sq-C*Phip(:),FOR.acc,FOR.maxit,Cholt,Chol,Phis(:));
            itsum=itsum+iter;
            if flag==0
                %      fprintf(' t=%g s It=%g E=%g\n',ttt,iter,err);
            else
                message(sprintf('no convergence (flag=%d)!',flag));
            end
            %Phigesamt=Phip+Phis;    % Phis nur ausgeborgt
            %             iphis=Phis(ia:di:ie,ja:dj:je,ka:ke);
            % !!!! noch besser alle, und dann besser integrieren...
            %             PHI(:,l)=iphis(:);
            phisi=Phis(ia:ie,ja:je,ka:ke);PHI(:,l)=phisi(:);
            %PHI(:,l)=reshape(,(ie-ia+1)*(je-ja+1)*(ke-ka+1),1);  % Nur für Sens. wichtigen Teil
            %mea=potxy(N.elec(:,1),N.elec(:,2),sq,Phis,X,Y,RA);
            el=N.elec;el(:,1)=el(:,1)-RA(1);el(:,2)=el(:,2)-RA(2);%el(l,1)=1;
            su=sqrt(sum(el.^2,2));su(l)=1;
            mea=1./su/(2*pi*sq); % Phip, nur für Oberfläche
            mea=mea+interp2(Y,X,Phis(:,:,1),N.elec(:,2),N.elec(:,1));
            mea(l)=0;
            if find(0>mea),
                [ l N.elec(l,:) length(find(0>mea))]
%                 pause(0.1);
            end
            MEA(l,:)=mea';  % Potentiale Nr. i bei Einspeisung in j
            mal=mal-1;
            if mal==0,
                waitbar(l/anzel,wb);
                mal=aller;
            end
        end
        close(wb)
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
        rez=(R-Rrez)*2./(R+Rrez);
        message(sprintf('Standard deviation of reciprocity %.2f%%',std(rez)*100));
        Rrez=Rrez(:).*N.k(:);
    end % pole-wise
    R=R(:).*N.k(:);
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
