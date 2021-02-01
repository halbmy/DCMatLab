%% synthetic model
rho=[1000 200 40];nm=length(rho);
taui=[0.1 0.1 0.03];mi=[0 0.01 0.01];
thk=[1 3];
% ab2=[1.5 2 2.5 3 4 5 6 8 10 12 15 20 25];
ab2=[1 2 5 10 20 50];
mn2=ones(size(ab2))*0.3;
rhoa=fwdschlum(rho,thk,ab2,mn2);nd=length(rhoa);
figure(1);subplot(1,2,1);loglog(rhoa,ab2,'x-');
axis ij tight;grid on
%% normal sensitivity
S1=ones(nd,nm);
for i=1:nm,
    rho1=rho; % get old
    rho1(i)=rho(i)*1.05; % change by 5%
    R1=fwdschlum(rho1,thk,ab2,mn2); % response of changed model
    S1(:,i)=(R1-rhoa)/(rho(i)*0.05); % change in response by change in model
end
%%
% f=[1 3 10];tau=[0.03 0.1]; % very small
f=[0.5 1 2 5 10 20];tau=[0.01 0.03 0.1 0.3];
% f=[0.3 0.7 1.4 3 6 12 24];tau=logspace(-3,0.5,20);
nf=length(f);nt=length(tau);
S1rho=S1*diag(rho);
WT1=f(:)*tau*pi*2;
W1=WT1./(WT1.^2+1);
S=kron(W1,S1rho);
%% spectral synthetic model (not used)
M=zeros(nm,nt);
M(2,2)=mi(2);
M(3,1)=mi(3);
RI=zeros(nd,nf);
RI(:)=S*M(:); %calculate
SC=diag(1./rhoa)*1000;
imagesc(SC*RI1);colorbar
%%
[L,C1,C2]=smooth2d1st(0:nt,(0:nm)*10,0);
errP=ones(nd,nf)*0.0001;
errRI=diag(rhoa)*errP;
D=diag(1./errRI(:));
RIn=RI+randn(nd,nf).*errRI;
%%
MI=ones(nm,nt)*0.001;
fRI=RI;fRI(:)=S*MI(:);
figure(1);
subplot(1,2,1);imagesc(M);caxis([0 max(mi)*1.2]);colorbar
subplot(1,2,2);imagesc(SC*RIn);cax=caxis;colorbar
figure(2);
subplot(1,2,1);imagesc(MI);caxis([0 max(mi)*1.2]);colorbar
subplot(1,2,2);imagesc(SC*fRI);caxis(cax);colorbar
%%
lam=1;
for i=1:8,
    DST=D*S*diag(MI(:));
    dM=(DST'*DST+L*lam)\(DST'*(D*(RIn(:)-fRI(:)))-L*log(MI(:))*lam);
    MI(:)=MI(:).*exp(dM*0.5);
    fRI(:)=S*MI(:);
    subplot(1,2,1);imagesc(MI);caxis([0 max(mi)*1.2]);colorbar
    set(gca
    subplot(1,2,2);imagesc(SC*fRI);caxis(cax);colorbar
    mean((D*(RIn(:)-fRI(:))).^2)
    pause(1.0);
end
%%
fprintf('tau2=%g(%g)\ttau3=%g(%g)\n',tau(MI(2,:)==max(MI(2,:))),...
    taui(2),tau(MI(3,:)==max(MI(3,:))),taui(3));
