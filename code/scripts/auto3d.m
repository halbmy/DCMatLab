global Model N S
[path,name,ext]=fileparts(datfile); 
if strcmp(lower(ext),'.pro'),
    N=readpro(datfile);
else
  if checkfile(datfile)==2,
      N=readinv3dfile(datfile);
  else
      N=readres3dinvfile(datfile);
  end
  N=getpseudos(N);
end
if (~isfield(N,'err'))&&(length(N.err)~=length(N.a)),
  if ~isfield(N,'i'), N.i=Err.curr*1e-3; end
  N.err=estimateerror(N,Err.proz,Err.umin*1e-6,N.i); % 3% plus 100myV/100mA
  fprintf('Error estimate: %.1f%% + %dmyV/%dmA: min=%.1f%% max=%.1f%%\n',...
  Err.proz,round(Err.umin),round(Err.curr),min(N.err)*100,max(N.err)*100);
end
if exist('DX'),
  [Model.x,Model.y,Model.z,Model.M,Model.Bg]=modelfromdata3d(N,DX,DX);    
else
  [Model.x,Model.y,Model.z,Model.M,Model.Bg]=modelfromdata3d(N);    
end
R=ones(size(N.r))*Model.Bg(1);
RMSerr=rms(N.r,R,INV.lolo);
CHIQ=chi2(N.r,R,N.err);
message(sprintf('RMS = %.1f Chi^2 = %.1f',RMSerr,CHIQ));
if exist('z')==1,
  Model.z=z;
  Model.M=ones(length(Model.x)-1,length(Model.y)-1,length(Model.z)-1)*Model.M(1,1,1);
  Model.Bg=ones(size(z))*Model.Bg(1);
end
if iscell(Model.M),
  S=modsens(Model);
else
  S=calcsens3d(Model.x,Model.y,Model.z,N);
end
if exist('Bg')==1,
  Model.Bg=Bg;oldR=R;INV.eind=0;
  if length(unique(Model.Bg))==1,
    R(:)=Model.Bg(1);
    RMSerr=rms(N.r,R,INV.lolo);
    CHIQ=chi2(N.r,R,N.err);
    message(sprintf('RMS = %.1f Chi^2 = %.1f',RMSerr,CHIQ));
  else
    for k=1:size(Model.M,3), Model.M(:,:,k)=Bg(min(k,length(Bg))); end
    dM=log(Model.M(:))-log(Bg(1));
    message('Calulating forward response of 1d model...');
    R=fwd3d1d(N,Model.Bg,Model.z); 
    RMSerr=rms(N.r,R,INV.lolo);
    CHIQ=chi2(N.r,R,N.err);
    message(sprintf('RMS = %.1f Chi^2 = %.1f',RMSerr,CHIQ));
    fprintf('Sensitivity updating by broyden method...');
    ddr=(log(R)-log(oldR)-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
    for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
    %S = S + ((log(R)-log(oldR)-S*dM(:))*(dM(:)'))/(dM(:)'*dM(:)); 
    fprintf('ready\n');
  end
end
Cov=Model.M;Cov(:)=sum(abs(S));
C=smooth3d2nd(Model.x,Model.y,Model.z);C=C'*C;
% Cov=sum(abs(S));C=spdiags(1./sqrt(Cov(:)),0,size(S,2),size(S,2));
P=1;
D=spdiags(1./log(1+N.err),0,length(N.err),length(N.err));
lam=10;
M0=Model.M;
if isfield(INV,'eind')&&(INV.eind),
  [Model,dM]=onedinv(N,Model); 
  message(['Layer resistivities: ' sprintf('%g ',Model.Bg)]);
  oldR=R;
  R=fwd3d1d(N,Model.Bg,Model.z); 
  RMSerr=rms(N.r,R,INV.lolo);
  CHIQ=chi2(N.r,R,N.err);
  message(sprintf('RMS = %.1f Chi^2 = %.1f',RMSerr,CHIQ));
  fprintf('Sensitivity updating by broyden method...');
  ddr=(log(R)-log(oldR)-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
  for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
  %S = S + ((log(R)-log(oldR)-S*dM(:))*(dM(:)'))/(dM(:)'*dM(:)); 
  fprintf('ready\n');
end
nM=Model;
running=1;lb=1;ub=10000;
while running,
    dR=log(N.r(:)-INV.lbound)-log(R(:)-INV.lbound);
    %dM=cglscdp(S,dR,lam,C,D,P);%,log(M(:)./M0(:)));
    [dM,INV.lam]=minvers(S,dR,INV,N.err,Model);
    if (INV.auto==1)&&(INV.const==1), INV.auto=2; end
    linesearch=1;fak=1;
    while linesearch,
      nM.M(:)=(Model.M(:)-INV.lbound).*exp(dM*fak)+INV.lbound;
      linesearch=(min(nM.M(:))<lb)|(max(nM.M(:))>ub);
      if linesearch, fak=fak*0.9; end
    end
    if fak<1, message(sprintf('Trusted region line search factor %.2f',fak)); end
    oldR=R;
    R=mfdfwd3d(nM,N,FOR);
    RMSerr=[RMSerr rms(N.r,R,INV.lolo)];
    CHIQ=[CHIQ chi2(N.r,R,N.err)];
    message(sprintf('RMS = %.1f Chi^2 = %.1f',RMSerr(end),CHIQ(end)));
    dchiq=CHIQ(end)-CHIQ(end-1);
    fprintf('Doing Broyden Update of Sensitivity...');
    ddr=(log(R)-log(oldR)-S*dM(:))/(dM(:)'*dM(:));ddr=ddr(:);
    for i=1:size(S,2), S(:,i)=S(:,i)+ddr*dM(i); end
    fprintf('ready\n');
    %S = S + ((log(R)-log(oldR)-S*dM)*(dM'))/(dM'*dM); % BROYDEN
    if (dchiq<0)|(length(CHIQ)<3),
        Model=nM;
        running=(abs(dchiq)/CHIQ(end)>0.05);
        if running==0, message('Stopping (dCHI^2<5 percent)'); end
    else
        message('Going back to old model and stopping');
        running=0;
    end
    if CHIQ(end)<0.9, running=0; end
end
outfile=strrep(datfile,ext,'.mod'); 
fid=fopen(outfile,'w'); 
ss='%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f\r\n'; 
if ~iscell(Model.M),
    for k = 1:size(Model.M,3),
        z1=Model.z(k);z2=Model.z(k+1);
        for j = 1:size(Model.M,2),
            y1=Model.y(j);y2=Model.y(j+1);
            for i = 1:size(Model.M,1),
                x1=Model.x(i);x2=Model.x(i+1);
                fprintf(fid,ss,x1,x2,y1,y2,z1,z2,Model.M(i,j,k));
            end
        end
    end      
end
fclose(fid);
figure;close;
h=figure('visible','off');
MAL.cbar=0;MAL.alpha=1;
[cmin,cmax]=draw3dmodel(Model,MAL);
print(h,'-depsc2',strrep(datfile,ext,'.eps'));
%clf;cbar(cmin,cmax,MAL.log);
%print(h,'-depsc2',strrep(datfile,ext,'-cbar.eps'));
dR=(1-R(:)./N.r(:))*100;
mal=MAL;mal.cauto=0;
mal.cmax=max([N.r(:);R(:)]);
mal.cmin=min([N.r(:);R(:)]);
plotprofiles(N,[],mal);
print(h,'-depsc2',strrep(datfile,ext,'-data.eps'));
plotprofiles(N,R,mal);
print(h,'-depsc2',strrep(datfile,ext,'-response.eps'));
