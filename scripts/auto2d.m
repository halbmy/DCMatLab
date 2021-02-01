%function auto2d(datfile,varargin)
load setup
Err=struct('proz',3,'umin',100e-6,'curr',0.1);

datfile='ggd/comb.dat'; zmax=30;nlay=11;xtra=5;
INV.lam=3;%FOR.zref=2;

datfile='wnr_0813.dat';INV.lam=3;Err.proz=10;Err.umin=100e-6;zmax=12;

global x z M N S
N=read2dfile(datfile);
if ~isfield(N,'i'), N.i=Err.curr; end
N.err=estimateerror(N,Err.proz,Err.umin,N.i); % 3% plus 100myV/100mA
fprintf('Error estimate: %.1f%% + %dmyV/%dmA: min=%.1f%% max=%.1f%%\n',...
Err.proz,round(Err.umin),round(Err.curr),min(N.err)*100,max(N.err)*100);
%[x,z,M]=modelfromdata2d(N);
z=zparam(N,nlay,zmax);
del=diff(N.elec(:,1));dx=min(del(find(del)))/2;
x=min(N.elec(:,1))-xtra*dx:dx:max(N.elec(:,1))+xtra*dx;
rhostart=median(N.r);
M=ones(length(x)-1,length(z)-1)*rhostart;
if ~exist('Lay'),
  Lay=M(1,1);
else
  M(:)=Lay(end);rho=Lay(1);
  for k=1:min(length(Lay),size(M,2)), M(:,k)=Lay(k); end
end
S=getsens2d(datfile,x,z,N);Cov=M;Cov(:)=sum(abs(S));
R=ones(size(N.r))*Lay(1);
C=smooth2d2nd(x,z);C=C'*C;
% Cov=sum(abs(S));C=spdiags(1./sqrt(Cov(:)),0,size(S,2),size(S,2));
P=pmatrix2d(x,z);P=1;
D=spdiags(1./log(1+N.err),0,length(N.err),length(N.err));
lam=10;if isfield(INV,'lam'), lam=INV.lam; end
M0=M;
RMSerr=rms(N.r,R,INV.lolo);
CHIQ=chi2(N.r,R,N.err);
tmf=norm(D*(log(N.r)-log(R)))^2+lam*norm(C*(log(M(:)./M0(:))))^2;
message(sprintf('RMS = %.1f Chi^2 = %.1f, TMF = %.1f',RMSerr,CHIQ,tmf));
if isfield(INV,'eind')&&(INV.eind),
  message('Estimating 1D-Starting Model...'); 
  rho=Lay(1);FOR.fillup=0;
  [Lay,M]=onedinv(N,z,M);
end
if length(unique(Lay))>1,
  message(['Layer resistivities: ' sprintf('%g ',Lay)]);
  dM=log(M(:))-log(rho);
  if max(N.elec(:,2))>0,
      R=dcfwd2d(x,z,M,Lay,N,FOR);
  else
      message('Calculating forward response of 1D model...');
      R=fwd2d1d(N,Lay,z);
  end 
  S = S + ((log(R)-log(rho)-S*dM)*(dM'))/(dM'*dM); % broyden 
  RMSerr=rms(N.r,R,INV.lolo);
  CHIQ=chi2(N.r,R,N.err);
  tmf=norm(D*(log(N.r)-log(R)))^2+lam*norm(C*(log(M(:)./M0(:))))^2;
  message(sprintf('RMS = %.1f Chi^2 = %.1f, TMF = %.1f',RMSerr,CHIQ,tmf));
end
Mref=M;nM=M;
running=1;
while running,
    dR=log(N.r(:)-INV.lbound)-log(R(:)-INV.lbound);
    dM0=log(M(:)-INV.lbound)-log(Mref(:)-INV.lbound);
    [dM,INV.lam]=invers(S,dR,INV,N.err,dM0);
    if (INV.auto==1)&&(INV.const==1), INV.auto=2; end
    nM(:)=(M(:)-INV.lbound).*exp(dM)+INV.lbound;
    oldR=R;
    R=dcfwd2d(x,z,nM,Lay,N,FOR);
    if isfield(INV,'linesearch')&&(INV.linesearch>0),
        tau=linesearch(N,oldR,R);
        if tau==0, break; end
        if tau<1, 
            message(sprintf('line search parameter %.2f  ',tau));
            nM(:)=(M(:)-INV.lbound).*exp(dM*tau)+INV.lbound;
            R=dcfwd2d(x,z,nM,Lay,N,FOR); 
        end
    end
    draw2dmodel(x,z,nM,MAL);drawnow;    
    RMSerr=[RMSerr rms(N.r,R,INV.lolo)];
    CHIQ=[CHIQ chi2(N.r,R,N.err)];
    tmf=[tmf norm(D*(log(N.r)-log(R)))^2+lam*norm(C*(log(M(:)./M0(:))))^2];
    message(sprintf('RMS = %.1f Chi^2 = %.1f, TMF = %.1f',RMSerr(end),CHIQ(end),tmf(end)));
    dchiq=CHIQ(end)-CHIQ(end-1);
    S = S + ((log(R)-log(oldR)-S*dM)*(dM'))/(dM'*dM); % BROYDEN
    message('Doing Broyden Update of Sensitivity');
    if (dchiq<0)|(length(CHIQ)<3),
        M=nM;
        running=(abs(dchiq)/CHIQ(end)>0.05);
        if running==0, message('Stopping (dCHI^2<5 percent)'); end
    else
        message('Going back to old model and stopping');
        running=0;
    end
    if CHIQ(end)<0.9, running=0; end
end
draw2dmodel(x,z,M,MAL);
[path,name,ext]=fileparts(datfile); 
outfile=strrep(datfile,ext,'.mod'); 
% modelexport2d(outfile,M,x,z);
