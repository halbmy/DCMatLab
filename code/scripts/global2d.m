lam=10;
datfile='rothsch\neu\dipol-dipolLcirc1-15.dat';
load setup
if exist(datfile)~=2, error(['Data file ' datfile ' does not exist!']); end
[fpath,fname,fext]=fileparts(datfile);
N=read2dfile(datfile);data=length(N.r);N.err=[];
if ~isfield(N,'err')||(length(N.err)~=data),
    N.err=estimateerror(N,2); end    
[S,x,z]=getsens2d(datfile);
if isempty(S), 
    [x,z]=modelfromdata2d(N,8); 
    S=getsens2d(datfile,x,z,N);
end
rho0=median(N.r);
M0=ones(length(x)-1,length(z)-1)*rho0;
R0=ones(size(N.r))*rho0;
M=M0;nM=M;
D=spdiags(1./log(1+N.err(:)),0,length(N.err),length(N.err));
L=smooth2d1st(x,z,1,0.5);
% C=smooth2d2nd(x,z,0,1,0.5);L=C'*C;
MAL=struct('log',1,'cauto',1,'cmin',80,'cmax',1000,'cbar',0);
mm=zeros(size(S,2),1);MM=mm'*L*mm;
M=M0;R=R0;
CHIQ=chi2(N.r,R,N.err,1);
TMF=CHIQ*data+lam*MM;
it=0;running=1;
while running,
    it=it+1;
    dM=cglscdp(S,log(N.r)-log(R),lam,L,D,1,mm);
    nM(:)=M(:).*exp(dM);
    oldR=R;R=dcfwd2d(x,z,nM,100,N,FOR);
    % line search procedure
    for i=1:20,
        taus(i)=i*0.05;appR=oldR.*exp(taus(i)*(log(R)-log(oldR)));
        nM(:)=M(:).*exp(taus(i)*dM);mm=log(nM(:));
        tmf(i)=chi2(N.r,appR,N.err,1)*data+lam*mm'*L*mm;
    end
    [xx,nn]=min(tmf);tau=taus(nn);
    M(:)=M(:).*exp(tau*dM);
    if tau<1,
        fprintf('tau=%.2f\n',tau);
        R=dcfwd2d(x,z,M,100,N,FOR);
    end
    S = S + ((log(R)-log(oldR)-S*dM)*(dM'))/(dM'*dM); % broyden 
    mm=log(M(:))-log(M0(:));
    CHIQ(end+1)=chi2(N.r,R,N.err,1);
    fprintf('CHQ = ');fprintf('%.1f ',CHIQ);fprintf('\r\n');
    MM=mm'*L*mm;
    TMF(end+1)=CHIQ(end)*data+lam*MM;
    fprintf('TMF = ');fprintf('%.1f ',TMF);fprintf('\r\n');
    running=(CHIQ(end)/CHIQ(end-1)<0.95)|(it<5)|(CHIQ(end)>10);
    draw2dmodel(x,z,M,MAL);drawnow;
    text(0,z(end)*1.1,sprintf('lam=%.1f Iteration %d \\chi^2=%.1f modfun=%.1f tmf=%.1f',...
        lam,it,CHIQ(end),MM(end),TMF(end)));
end
modelexport2d(strrep(datfile,fext,'-g.mod'),M,x,z);
return
MAL.cbar=1;
Cov=M;Cov(:)=sum(abs(S));
draw2dmodel(x,z,M,MAL,Cov);drawnow;
text(0,z(end)*1.1,sprintf('lam=%.1f Iteration %d \\chi^2=%.1f modfun=%.1f tmf=%.1f',lam,it,CHIQ(end),MM(end),TMF(end)));
set(gcf,'PaperpositionMode','auto');
print(gcf,'-depsc2',strrep(datfile,fext,'-g.eps'));
dos(['epstopdf ' strrep(datfile,fext,'-g.eps')]);
