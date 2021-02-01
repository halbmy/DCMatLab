global Mod N FIX ZZ datfile
nmz=size(Mod.M,2);
FIX=zeros(size(Mod.M));
ZZ=zeros(size(Mod.M)-[0 1]);
xm=Mod.x(1:end-1)+diff(Mod.x)/2;
xd=load(strrep(datfile,'.dat','.xd'));
zz=interp1(xd(:,1),xd(:,2),xm,'nearest','extrap');

for i=1:length(zz),
   k=max(find(Mod.z<=zz(i)))-1;
   if k>0, 
       FIX(i,1:min(k,nmz))=-1; 
%        ZZ(i,min(k,nmz-1))=1;
%        ZZ(i,min(k,nmz-2)+(0:1))=1;
       dz=Mod.z(k+2)-Mod.z(k+1);
       if (zz(i)-Mod.z(k+1))/dz>0.4, ZZ(i,min(k+1,nmz-1))=1; end
       if (Mod.z(k+2)-zz(i))/dz>0.4, ZZ(i,min(k,nmz-1))=1; end
   end
end

return
xd=load(strrep(datfile,'.ohm','.xd'));
plot(xd(:,1),xd(:,2),'r-',N.elec(:,1),N.elec(:,2),'b-');axis ij
