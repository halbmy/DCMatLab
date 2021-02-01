function [phi1,tau,a,a1]=removeemcoupling(f,phi,dphi,opt)

% REMOVEEMCOUPLING - remove EM coupling by weighted Debye fit
% phi1 = removeemcoupling(f,phi[,opt])
% [phi1,tau,spect_charg] = ...
% f..frequency vector,phi..phase vector
% opt..option structure with fields
%      fdamp - damping frequency (2000Hz)
%      lam - regularization/smoothness parameter (1)
%      cexp - coupling exponent (2)
%      taumin - minimum relaxation time in s (1e-6)
%      taumax - maximum relaxation time in s (1 / 2 pi fmax)
%      ntau - number of relaxation times (50)

if nargin<4, opt=[]; end
if nargin<3, dphi=ones(size(phi))*1e-4; end
taumin=1e-6; % typical for EM
taumax=1/min(f)/2/pi;
ntau=50;
f1=2e3;
lam=1;
cexp=2;
if isfield(opt,'taumin'), taumin=opt.taumin; end
if isfield(opt,'taumax'), taumax=opt.taumax; end
if isfield(opt,'ntau'), ntau=opt.ntau; end
if isfield(opt,'lam'), lam=opt.lam; end
if isfield(opt,'fdamp'), f1=opt.fdamp; end
if isfield(opt,'cexp'), cexp=opt.cexp; end
tau=logspace(log10(taumin),log10(taumax),ntau)'; % tau discretisation
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
one=ones(ntau-1,1);
C=spdiags([-one one],[0 1],ntau-1,ntau); % 1D smoothness matrix
a=(DG'*DG+(C'*C)*lam)\(DG'*(D*r2r)); % spectral relaxation
w1t=f1*2*pi*tau; % damping frequency
a1=a.*w1t.^cexp./(w1t.^cexp+1);  % damped spectral relaxation
phi1=asin(G*a1); % resulting phase spectrum
