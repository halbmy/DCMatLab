function kl = getkellipse2d(N,l,d)

% GETKELLIPSE - Get geometric factor for half-ellipsoidal electrodes
% k = getkellipse(Data,l,d)

%program kfaktor_linie_dc2dinvres.m
%Berechnet k-Kaktor für linienhafte Elektroden unter Beruecksichtigung des
%Einflusses der Stecktiefe der Elektroden bei 2d Anordnung (d.h. z=konst, 
%x,y beliebig
%Elektrodenkonfigurationen muessen im dc2dinvres Format vorliegen
%kfaktoren werden in variable konf geschrieben, und in dc2dinvres
%uebernommen und mit abgespeichert
%damit diese dann beruecksichtigt werden, muessen die rhoa entfernt werden,
%neue rhoa werden dann aus u i k berechnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3, d=0.0045; end %input('Elektrodendurchmesser [m]: ');
if nargin<2, l=0.03; end %input('Elektrodenstecktiefe [m]: '); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if l==0,
    kl=getkonf(N);
else
    l=max(l,d*0.501);
    am=sqrt(sum((N.elec(N.a,:)-N.elec(N.m,:)).^2,2));
    fb=find(N.b);
    bm=sqrt(sum((N.elec(N.b(fb),:)-N.elec(N.m(fb),:)).^2,2));
    fn=find(N.n);
    an=sqrt(sum((N.elec(N.a(fn),:)-N.elec(N.n(fn),:)).^2,2));
    fbn=find(N.b.*N.n);
    bn=sqrt(sum((N.elec(N.b(fbn),:)-N.elec(N.n(fbn),:)).^2,2));
    kl=kellipse(l,d,am);
    if ~isempty(fb), kl(fb)=kl(fb)-kellipse(l,d,bm); end
    if ~isempty(fn), kl(fn)=kl(fn)-kellipse(l,d,an); end
    if ~isempty(fbn), kl(fbn)=kl(fbn)+kellipse(l,d,bn); end
    kl=1./kl;
end


function k = kellipse(l,d,r1);
% function[kl] = kellipse(l,d,r1,r2,r3,r4);
%FKLINIE berechnet k-Faktor fuer Stabelektroden der Länge l, Durchmesser d,
%nach Formel Militzer, Sommerfeld

nn=32;
z=[0:(l/nn):l]; %Integration des Potentials in 16 Schritten ueber Sondenlaenge
z=(0.5:0.5:nn-0.5)/nn*l; %Integration des Potentials in 16 Schritten ueber Sondenlaenge
e=sqrt(l.^2-d^2/4); %(alternativ: e=l) ...Brennpunktabstand s. Militzer

k=zeros(size(r1));
for i=1:length(z),
    k=k+log(abs((z(i)+e+sqrt(r1.^2+(z(i)+e).^2))./(z(i)-e+sqrt(r1.^2+(z(i)-e).^2))));
end
k=k/(4*pi*e)/length(z);
% k=(4*pi*e).*log(abs((z+e+sqrt(r1^2+(z+e).^2))./(z-e+sqrt(r1^2+(z-e).^2))));
% k=(4*pi*e).*(...
%     log(abs( (z+e+sqrt(r1^2+(z+e).^2))./(z-e+sqrt(r1^2+(z-e).^2)) ))-...
%     log(abs( (z+e+sqrt(r2^2+(z+e).^2))./(z-e+sqrt(r2^2+(z-e).^2)) ))-...
%     log(abs( (z+e+sqrt(r3^2+(z+e).^2))./(z-e+sqrt(r3^2+(z-e).^2)) ))+...
%     log(abs( (z+e+sqrt(r4^2+(z+e).^2))./(z-e+sqrt(r4^2+(z-e).^2))
%     ))).^(-1);