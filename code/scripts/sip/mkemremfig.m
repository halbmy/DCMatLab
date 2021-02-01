fname='Eng03\ENG03-17.83-Kl-2.dat';
[f,R,phi,dR,dphi]=read1resfile(fname);
subplot(211);set(gca,'FontSize',18);
semilogx(f,phi*1000,'bx-');
axis tight;grid on;
xlabel('f in Hz');ylabel('\phi in mrad');
for i=1:length(f), line([1 1]*f(i),[-1 1]*dphi(i)*1000+phi(i)*1000,'Color','blue'); end
tau=logspace(-6,2,100)';
my=debyedecomp(f,phi,dphi,tau,10);
subplot(212);set(gca,'FontSize',18);
semilogx(tau,a,'bx-');
axis tight;grid on;
xlabel('\tau in s');ylabel('\my');
