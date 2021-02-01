%% load data and make model
A=load('ab2mn2rhoaphi-all.txt');
ab2=A(:,1);mn2=A(:,2);rhoa=A(:,3);
f=[0.73 1.46 2.93 5.86 11.72 23.44 46.88];% 93.75];
tau=logspace(-3.5,0,15);
nf=length(f);nt=length(tau);
thk=[0.5 1 1.5 2 3 4 4.5 5]; %lt. Bohrung ENG03
z=[0 cumsum(thk) sum(thk)+3];zm=[z(1:end-1)+[thk/2 1]];
rho=[1500 800 400 400 400 400 20 200 10]
[rho,R]=invfile1drho(A(:,1:3),thk,0.00001,rho);
P=A(:,4:end-1);
dd=18:19;P(dd,:)=[];ab2(dd)=[];mn2(dd)=[];rhoa(dd)=[];R(dd)=[];
nr=length(rhoa);nc=length(rho);
%% compute unit sensitivity matrix
S1=ones(nr,nc);fak=0.1;
for i=1:nc,
    rho1=rho; % get old
    rho1(i)=rho(i)*(1+fak); % change by 5%
    R1=fwdschlum(rho1,thk,ab2,mn2); % response of changed model
    S1(:,i)=(R1-R)/(rho(i)*fak); % change in response by change in model
end
%% computer full Jacobian matrix
S1rho=S1*diag(rho);
WT1=f(:)*tau*pi*2;
W1=WT1./(WT1.^2+1);
S=kron(W1,S1rho);
%% transform data and estimate error
RI=diag(R)*P/1000;
% P([18:19],:)=0;
P(P<0)=0;P(16:end,6)=0;P(14:end,7)=0;
errP=ones(size(P))*0.3;
errP(:,end-1:end)=errP(:,end-1:end)*[2 0;0 3];
errP(P==0)=1e6;
errRI=diag(R)*errP/1000;
D=diag(1./errRI(:));
SC=diag(1./R)*1000;%
%% actual inversion
[L,C1,C2]=smooth2d1st((0:nc),(0:nt),0,1);
MY=ones(nc,nt)*0.002;
fRI=RI;fRI(:)=S*MY(:);
cax=[0 12];
subplot(2,3,1);imagesc(SC*RI);alpha(double(P>0));caxis(cax);;colorbar
set(gca,'YTick',1:nr,'YTickLabel',num2strcell(ab2),'XTick',1:nf,'XTickLabel',num2strcell(rndig(f)));
title('\phi_{obs} in mrad');xlabel('f in Hz');ylabel('AB/2 in m');
subplot(2,3,4:5);imagesc(MY);colorbar
set(gca,'XTick',1:nt,'XTickLabel',num2strcell(rndig(tau,2)),'YTick',0.5:nc+0.5,'YTickLabel',num2strcell(z));
set(gca,'YTick',1:nr,'YTickLabel',num2strcell(ab2),'XTick',1:nf,'XTickLabel',num2strcell(f));
xlabel('depth in m');ylabel('\tau in ms');
subplot(2,3,2);imagesc(SC*fRI);alpha(double(P>0));caxis(cax);colorbar
set(gca,'YTick',1:nr,'YTickLabel',num2strcell(ab2),'XTick',1:nf,'XTickLabel',num2strcell(f));
title('\phi_{pred} in mrad');xlabel('f in Hz');ylabel('AB/2 in m');
%
lam=0.1;
for i=1:8,
    DST=D*S*diag(MY(:));
    dM=(DST'*DST+L*lam)\(DST'*(D*(RI(:)-fRI(:)))-L*log(MY(:))*lam);
    MY(:)=MY(:).*exp(dM*i/(i+1));
    fRI(:)=S*MY(:);
%     subplot(2,3,4:5);imagesc(MY);caxis([0 max(MY(:))]);colorbar;
    subplot(2,3,4:5);imagesc(log(MY));colorbar;
    set(gca,'XTick',1:2:nt,'XTickLabel',num2strcell(rndig(tau(1:2:end)*1000,2)),'YTick',0.5:nc+0.5,'YTickLabel',num2strcell(z));
    ylabel('depth in m');xlabel('\tau in ms');title('spectral chargeability');
    subplot(2,3,2);imagesc(SC*fRI);alpha(double(P>0));caxis(cax);colorbar
    set(gca,'YTick',1:nr,'YTickLabel',num2strcell(ab2),'XTick',1:nf,'XTickLabel',num2strcell(rndig(f)));
    title('\phi_{pred} in mrad');xlabel('f in Hz');ylabel('AB/2 in m');
    deltaPhi=SC*(fRI-RI);
    subplot(2,3,3);imagesc(deltaPhi);alpha(double(P>0));caxis(minmax(deltaPhi(P>0)));colorbar
    set(gca,'YTick',1:nr,'YTickLabel',num2strcell(ab2),'XTick',1:nf,'XTickLabel',num2strcell(rndig(f)));
    title('\Delta\phi in mrad');xlabel('f in Hz');ylabel('AB/2 in m');
    fprintf('Iter %d, chi^2=%.2f\n',i,mean((D*(RI(:)-fRI(:))).^2));
    pause(1.0);
end
%%
subplot(2,3,6);
maxMY=max(MY,[],2);nMY=diag(1./maxMY)*MY;
ind=[3 5 7 8 9];plot(log(tau),nMY(ind,:),'x-');legend(num2strcell(zm(ind)));
axis tight;grid on;set(gca,'XTickLabel',num2strcell(rndig(exp(get(gca,'XTick'))*1000)));
% set(gcf,'Renderer','OpenGL')
%% resolution analysis
GI=inv(DST'*DST+L*lam)*(DST');
MDM=GI'*GI;
DMYrel=MY;
DMYrel(:)=1./sqrt(diag(MCM));
showftaumodel(DMYrel,z,tau)
%%
RM=inv(DST'*DST+L*lam)*(DST'*DST);
RMd=MY;RMd(:)=RM(:,6*9+5);
showftaumodel(RMd,z,tau)
imagesc(RMd);colorbar