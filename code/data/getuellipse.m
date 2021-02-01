function u = getkellipse2d(pos,l,d)

% GETUELLIPSE - Get analytical potential for half-ellipsoidal electrodes
% u = getuellipse(Data,l,d)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3, d=0.0045; end %input('Elektrodendurchmesser [m]: ');
if nargin<2, l=0.03; end %input('Elektrodenstecktiefe [m]: '); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e=sqrt(l.^2-d^2/4); %(alternativ: e=l) ...Brennpunktabstand s. Militzer
z=pos(:,end);
x=sqrt(sum(pos(:,1:end-1).^2,2));
u=log(abs((z+e+sqrt(x.^2+(z+e).^2))./(z-e+sqrt(x.^2+(z-e).^2))))/(4*pi*e);
