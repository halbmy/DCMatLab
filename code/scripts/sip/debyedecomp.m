function [m,tau,phiout]=debyedecomp(f,phi,dphi,tau,lam,damp)

% DEBYEDECOMP - Debye relaxation decomposition
%               using positivity constraints
% [m,tau] = debyedecomp(f,phi)
% [m,tau] = debyedecomp(f,phi,dphi)
% [m,tau] = debyedecomp(f,phi,dphi,lam)
% m = debyedecomp(f,phi,dphi,tau,lam)
% m = debyedecomp(f,phi,dphi,tau)
% m..spectral chargeability for relaxation times..tau
% phi..phases for frequencies..f with phase errors..dphi
% lam..degree of smoothness (1-1000, default 100)

if nargin<6, damp=0; end
if nargin<5, lam=100; end
if (nargin==4)&&(length(tau)==1), lam=tau; end
if (nargin<4)||(length(tau)<2),
    taumax=1/min(f)/2/pi*4;
    taumin=1/max(f)/2/pi/8;
    ntau=100;
    tau=logspace(log10(taumin),log10(taumax),ntau)'; % tau discretisation
else
    ntau=length(tau);
end
if (nargin<3)||isempty(dphi), dphi=ones(size(phi))*1e-4; end
G=ones(length(f),length(tau)); % kernel function
for i=1:length(tau),
    wt=f*2*pi*tau(i);
    g=wt./(wt.^2+1);
    G(:,i)=g;
end
r2r=sin(phi(:)); % relative imaginary conductivity
dr2r=dphi(:).*cos(phi(:)); % error in r2r from error in phase
D=diag(1./dr2r); % error weighting matrix
DG=D*G; % weighted
mmin=1e-5;dampval=0.005;
if damp>0,
    one=ones(ntau,1);
    C=spdiags(one,0,ntau,ntau);
else
    one=ones(ntau-1,1);
    C=spdiags([-one one],[0 1],ntau-1,ntau); % 1D smoothness matrix
    if damp<0,
       C(ntau,1)=dampval;
       C(ntau+1,ntau)=dampval;
    end
end
%%
mstart=mean(r2r./sum(G,2));
m=ones(length(tau),1)*mstart;
d1=zeros(100,1);
%%
for i=1:5,
    DGT=DG*diag(m);
    dd=r2r-G*m;
%     dm=(DGT'*DGT+(C'*C)*lam)\(DGT'*(D*dd)-C'*(C*log(m))*lam*(1-damp)); % spectral relaxation
    rhs=C*log(m);
    if damp>0, rhs=zeros(ntau,1); end
    if damp<0, rhs(end-1:end)=rhs(end-1:end)-log(mmin); end 
    dm=(DGT'*DGT+(C'*C)*lam)\(DGT'*(D*dd)-C'*(rhs)*lam);
    for t=1:100, % line search
        m1=m.*exp(dm*t/100);
        d1(t)=norm(D*(r2r-G*m1));
    end
    if damp>0, lam=lam*0.8; end
    [dopt,topt]=min(d1);
    m=m.*exp(dm*topt/100);
end
phiout=asin(G*m);
if (nargout==2)&&(nargin>=4), %tau given
    tau=phiout;
end
if nargout<1,
    subplot(211);semilogx(f,phi*1000,'bx-',f,asin(G*m)*1000,'ro-');
    grid on;xlabel('f in Hz');ylabel('\phi in mrad');
    subplot(212);semilogx(tau,m,'.-');grid on;
    xlabel('\tau in s');ylabel('spectral charcheability');
end