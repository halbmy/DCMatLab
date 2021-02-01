%% synthetic model
rho=[1000 200 40];nm=length(rho);
taui=[0.1 0.08 0.02];mi=[0 0.01 0.01];
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
% f=[0.5 1 2 5 10 20];tau=[0.01 0.03 0.1 0.3];
f=[0.3 0.7 1.4 3 6 12 24];tau=logspace(-3,0.5,20);
nf=length(f);nt=length(tau);
A=repmat(f,nd,1);fvec=A(:);
A=repmat(tau,nm,1);tvec=A(:);
WT=fvec*tvec'*pi*2;
rhovec=repmat(rho,1,nt);
rhoavec=repmat(rhoa,1,nf);
S=diag(1./rhoavec)*repmat(S1,nf,nt).*WT./(WT.^2+1)*diag(rhovec);
subplot(1,2,2);imagesc(S);
%% spectral synthetic model (not used)
M=zeros(nm,nt);
M(2,3)=mi(2);
M(3,2)=mi(3);
P1=zeros(nd,nf);
P1(:)=S*M(:); %calculate
imagesc(P1);colorbar
%% check whether it is the same by summation
P=zeros(nd,nf);
for i=1:nf,
    for j=1:nm,
        wt=f(i)*taui(j)*2*pi;
        rhoi=rho(j)*mi(j)*wt/(wt^2+1);
        P(:,i)=P(:,i)+S1(:,j)*rhoi./rhoa(:);        
    end
end
subplot(1,2,1);imagesc(P);colorbar
%%
[L,C1,C2]=smooth2d1st(0:nt,(0:nm)*10,0);
errP=ones(nd,nf)*0.0005;
D=diag(1./errP(:));
Pn=P+randn(nd,nf).*errP;
MI=zeros(nm,nt);
MI(:)=(S'*S+L*0.05)\(S'*Pn(:));
[tau;MI]
imagesc(MI);colorbar
%%
MI=ones(nm,nt)*mean(Pn(:));
fP=Pn;fP(:)=S*MI(:);
figure(1);
subplot(1,2,1);imagesc(M);caxis([0 max(mi)*1.2]);colorbar
subplot(1,2,2);imagesc(Pn);cax=caxis;colorbar
figure(2);
subplot(1,2,1);imagesc(MI);caxis([0 max(mi)*1.2]);colorbar
subplot(1,2,2);imagesc(fP);caxis(cax);colorbar
%%
lam=1;
for i=1:8,
    DST=D*S*diag(MI(:));
    dM=(DST'*DST+L*lam)\(DST'*(D*(Pn(:)-fP(:)))-L*log(MI(:))*lam);
    MI(:)=MI(:).*exp(dM*0.5);
    fP(:)=S*MI(:);
    subplot(1,2,1);imagesc(MI);caxis([0 max(mi)*1.2]);colorbar
    subplot(1,2,2);imagesc(fP);caxis(cax);colorbar
    mean((D*(Pn(:)-fP(:))).^2)
    pause(1.0);
end
%%
fprintf('tau2=%g(%g)\ttau3=%g(%g)\n',tau(MI(2,:)==max(MI(2,:))),...
    taui(2),tau(MI(3,:)==max(MI(3,:))),taui(3));
%%




return
