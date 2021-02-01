% k=1/39.8;
clc
fname='Eng03\ENG03-17.83-Kl-2.dat';
[f,R,phi,dR,dphi]=read1resfile(fname);
% dphi(4)=1e8;
subplot(211);errorbar(log10(f),phi,phi-dphi,phi+dphi);grid on;
%%
semilogx(f,phi,'bx-');
for i=1:length(f), line([1 1]*f(i),[-1 1]*dphi(i)+phi(i),'Color','blue'); end
%%
tau=logspace(-6,2,100)';
nt=length(tau);
G=ones(length(f),length(tau));
subplot(212);
for i=1:length(tau),
    wt=f*2*pi*tau(i);
    g=wt./(wt.^2+1);
    G(:,i)=g;
    semilogx(f,g);
    hold on;
end
hold off
%
r2r=sin(phi)';
dr2r=dphi'.*cos(phi');
D=diag(1./dr2r);
DG=D*G;
one=ones(nt-1,1);
C=spdiags([-one one],[0 1],nt-1,nt);
lam=1;
a=(DG'*DG+(C'*C)*lam)\(DG'*(D*r2r));
subplot(212);semilogx(tau,a,'bx-');
subplot(211);semilogx(f,phi*1000,'bx-');
%%
w1=1e4;tw1=tau*w1;
c=2;fk=tw1.^c./(tw1.^c+1);
a1=a.*fk;
phi1=asin(G*a1);
subplot(212);semilogx(tau,a,'bx-',tau,a1,'ro-');grid on;
xlabel('\tau in s');ylabel('spectral chargeability');
subplot(211);semilogx(f,phi*1000,'bx-',f,phi1*1000,'ro-');grid on;
xlabel('f in Hz');ylabel('\phi in mrad');
xlim(minmax(f));ylim([0 45]);grid on;
% ylim(minmax(phi1)*1000);
