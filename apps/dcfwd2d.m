function [R,Rrez,MEA]=dcfwd2d(x,z,M,Lay,N,FOR)

% DCFWD2D - 2D DC Forward Calculation (polwise)
% R = dcfwd2d(x,z,M,Lay,N,FOR)
% R = dcfwd2d(Mod,N,FOR)
% [R,Rrez,MEA]=dcfwd2d(x,z,M,Lay,N,FOR)
% x,z..Vectors of nodes
% M..Cell Resistivities[length(x)*length(z)]
% Lay..Layered Model (or background res.)
% N..struct of electrode numbers(a,b,m,n)
%    resistivities(r), konfiguration factors(k)
%     and electrode positions(elec)
% FOR..structure of possible fields
%      rand..border cells[3]
%      prolong..prolonging factor[5]
%      refine..refining factor for x[0..automatic]
%      acc..Solving accuration[1e-4]
%      tol..accuracy of preconditioner[1e-4]
%      maxit..maximum number of iterations[200]

if isstruct(x),
    if nargin<2, error('Too less input arguments!'); end
    N=z;
    if nargin>2, FOR=M; else FOR=struct('rand',4,'prolong',4,'refine',0); end
    if isfield(x,'Lay'), Lay=x.Lay; else Lay=0; end
    if isfield(x,'M'), M=x.M; else error('specify model!'); end
    if isfield(x,'z'), z=x.z; else error('specify z!'); end
    if isfield(x,'x'), x=x.x; else error('specify x!'); end
else
    if nargin<5, error('Too less input arguments!'); end
    if nargin<6, FOR=struct('rand',4,'prolong',4,'refine',0); end
end
if ~isfield(FOR,'rand'), FOR.rand=4; end
if ~isfield(FOR,'prolong'), FOR.prolong=4; end
if ~isfield(FOR,'acc'), FOR.acc=1e-4; end
if ~isfield(FOR,'tol'), FOR.tol=1e-4; end
if ~isfield(FOR,'maxit'), FOR.maxit=200; end
if ~isfield(FOR,'refine'), FOR.refine=0; end
if ~isfield(FOR,'zusatz'), FOR.zusatz=2; end
if ~isfield(FOR,'silent'), FOR.silent=0; end
zref=1;
if isfield(FOR,'zref'), zref=FOR.zref; end
if isfield(FOR,'fillup')&&(FOR.fillup),
    Lay=zeros(length(z),1);
end

if zref>1,
  Z=ones((length(z)-1)*zref+1,1);
  for k=1:zref, Z(k:zref:end-1)=z(1:end-1)+diff(z)*(k-1)/zref; end
  Z(end)=z(end);
else
    Z=z(:);
end
if FOR.refine==0,
    ref=round(min(diff(x))/min(diff(Z)));
    if ref<1, ref=1; end
else
    ref=FOR.refine;
end
X=ones((length(x)-1)*ref+1,1);
for k=1:ref, X(k:ref:end-1)=x(1:end-1)+(k-1).*diff(x)/ref; end
X(end)=x(end);

if FOR.rand>0
    if FOR.zusatz>0, Prolong=(1:FOR.zusatz)'; else Prolong=[]; end
    zProlong=Prolong;
    if max(N.elec(:,2))>0, zProlong=Prolong; end % !!! noch zu verbessern
    pp=1;
    for ll=1:FOR.rand
        pp=pp*FOR.prolong;
        if ~isempty(Prolong), Prolong=[Prolong;Prolong(end)+pp];
        else Prolong=pp; end
        if ~isempty(zProlong), zProlong=[zProlong;zProlong(end)+pp];
        else zProlong=[zProlong;pp]; end
    end
    X=[min(X)-(X(2)-X(1))*flipud(Prolong);X;max(X)+(X(end)-X(end-1))*Prolong];
    Z=[Z(:);max(Z)+(Z(end)-Z(end-1))*zProlong];
    zrand=length(zProlong);
    while(Z(end)/sqrt(FOR.prolong)>X(end)-X(1)), 
        Z(end)=[];zrand=zrand-1;
    end
    if find([(diff(X)<=0);(diff(Z)<=0)]), error('Something is wrong!'); end
    if(~FOR.silent),
        message(sprintf('Adding boundaries x=%.1f-%.1f  z=%.1f-%.1f',...
        min(X),max(X),min(Z),max(Z)));
    end
end

if Lay(1)==0, sig=0.01; else sig=1/Lay(1); end
Sigmac=ones(length(Z)-1,length(X)-1)*sig;
if length(Lay)>1,
    for k=1:length(Lay)-1, 
        if Lay(k)>0, % background
            sig=1./Lay(k);
            for kk=1:zref, Sigmac((k-1)*zref+kk,:)=sig; end
        else % boundary
            sig=1/M(1,k);    
            sig2=1./M(end,k);
            for kk=1:zref,
                Sigmac((k-1)*zref+kk,:)=sig;    
                Sigmac((k-1)*zref+kk,end-FOR.zusatz-FOR.rand+1:end)=sig2;
            end
        end
    end
end
for k=1:ref,
    for i=1:zref,
        Sigmac(i:zref:end-zrand,...
            FOR.zusatz+FOR.rand+k:ref:end-FOR.rand-FOR.zusatz)=1./M';
    end
end
% lLay=length(Lay)-1;
if Lay(end)>0, % background
    sig=1/Lay(end);
    Sigmac(size(M,2)*zref+1:end,:)=sig;
else % boundary
    for k=size(M,2)*zref+1:size(Sigmac,1),
        Sigmac(k,:)=Sigmac(size(M,2)*zref,:);
    end
end
% figure(10);imagesc(log10(1./Sigmac));colorbar horiz

if(~FOR.silent),
    message(sprintf('Forward Calculation with %dx%d=%d cells (refine=%d,%d)',...
        size(Sigmac,1),size(Sigmac,2),numel(Sigmac),ref,zref));
end
if find(N.elec(:,2)~=0), Sigmac(end,end)=Sigmac(end,end)*1.000001; end
t0=clock;
%MEA=dc2dmea(X,Z,Sigmac,N.elec(:,1),N.elec(:,2));
if isfield(FOR,'method')&&(FOR.method==1),%as measured
    % eigentlich gesondertes geo2d mit 2 Quellen
    un=unique([N.a;N.b]);
    MEA=geo2d(X,Z,Sigmac,N.elec(un,1),N.elec(un,2),N.elec(:,1),N.elec(:,2));
    map=ones(size(N.elec,1));
    map(un)=1:length(un);
elseif isfield(FOR,'method')&&(FOR.method==2),%reverse
    % eigentlich gesondertes geo2d mit 2 Quellen
    un=unique([N.m;N.n]);
    MEA=geo2d(X,Z,Sigmac,N.elec(un,1),N.elec(un,2),N.elec(:,1),N.elec(:,2));
    map=ones(size(N.elec,1));
    map(un)=1:length(un);
else
    MEA=geo2d(X,Z,Sigmac,N.elec(:,1),N.elec(:,2));
    for l=1:size(MEA,1), MEA(l,l)=0; end
    map=1:size(N.elec,1);
end
fi=find(MEA<0);
if (length(fi)>0)&(~FOR.silent),
    message(sprintf('Found %d=%.f%% values below zero',length(fi),length(fi)/numel(MEA)*100));
    %MEA(fi)=0;
end
%maxmea=max(max((MEA-MEA')./MEA));
%MEA=(MEA+MEA')/2;
R=zeros(length(N.a),1);Rrez=R;
for l = 1:length(N.a),
    rho=0;rez=0;
    if N.m(l)>0,
        if N.a(l)>0,
            rho=rho+MEA(N.m(l),map(N.a(l)));
            rez=rez+MEA(N.a(l),map(N.m(l)));
        end
        if N.b(l)>0,
            rho=rho-MEA(N.m(l),map(N.b(l)));
            rez=rez-MEA(N.b(l),map(N.m(l)));
        end
    end
    if N.n(l)>0,
        if N.a(l)>0,
            rho=rho-MEA(N.n(l),map(N.a(l)));
            rez=rez-MEA(N.a(l),map(N.n(l)));
        end
        if N.b(l)>0,
            rho=rho+MEA(N.n(l),map(N.b(l)));
            rez=rez+MEA(N.b(l),map(N.n(l)));
        end
    end
    R(l)=rho*N.k(l);
    Rrez(l)=rez*N.k(l);
end
% [R,Rrez]=collectrhoa(N,MEA);
if isfield(FOR,'method')&&(FOR.method==1),
    if ~(FOR.silent), message(sprintf('ready(%.1fs)',etime(clock,t0))); end
elseif isfield(FOR,'method')&&(FOR.method==2),
    if ~(FOR.silent), message(sprintf('ready(%.1fs)',etime(clock,t0))); end
    R=Rrez;
else
    rez=2*(R-Rrez)./(R+Rrez);
    maxrez=max(abs(rez));
    stdrez=std(rez);
    if nargout<2,
        RR=sqrt(abs(R.*Rrez));
        %fi=find(Rrez<=0);RR(fi)=R(fi);
        %fi=find(R<=0);RR(fi)=Rrez(fi);
        fi=find((N.elec(N.m,2)>0)&(N.b==0));RR(fi)=R(fi);
        fi=find((N.elec(N.a,2)>0)&(N.n==0));RR(fi)=Rrez(fi);
        R=RR;
    end
    if ~(FOR.silent),
        message(sprintf('ready(%.1fs), reciprocity error: sd=%.1f%%, max %.1f%%',...
            etime(clock,t0),stdrez*100,maxrez*100));
    end
end
