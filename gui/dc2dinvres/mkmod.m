global Mod
%Mod=modelimport2d
mal=struct('cauto',0,'clog',1,'cmin',1,'cmax',100,'high',7,'canot','Ohmm');
Mod1=[];Mod1.z=Mod.z;
nx=500;
Mod1.x=Mod.x(1:nx+1);
Mod1.M=Mod.M(1:nx,:);
patch2dmodel(Mod1.x,Mod1.z,Mod1.M,mal);
hold on;
plot(N.elec(:,1),N.elec(:,2),'w.');
plot(xd(:,1),xd(:,2),'k-');
hold off
epsprint(3,'line3-6000-roll1',1)