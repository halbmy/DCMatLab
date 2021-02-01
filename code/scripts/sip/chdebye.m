% fname='ENG03_17m29\ENG03-17.29-Kl-1.dat';
% fname='ENG03_17m83\ENG03-17.83-Kl-8.dat';
% fname='Eng03_4m15\ENG03-4.15-Kl-1.dat';
% fname='Eng03_4m55 (Torf)\ENG03-4.55-Kl-1.dat';
fname='Eng03_16m15\ENG03-16.15-Kl-1.dat';
% fname='Eng03_2m35\ENG03-2.35-3.dat';
[f,R,phi,dR,dphi]=read1resfile(fname);
% subplot(211);errorbar(log10(f),phi,phi-dphi,phi+dphi);grid on;
%%
phi1=removeemcoupling(f,phi,dphi,struct('fdamp',1e3));
subplot(211);semilogx(f,phi1);grid on;
debyedecomp(f,phi1,dphi,10);