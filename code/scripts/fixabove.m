global Mod N FIX ZZ
nmz=size(Mod.M,2);
if ~isequal(size(Mod.M),size(FIX)), FIX=zeros(size(Mod.M)); end
if ~isequal(size(Mod.M)-[0 1],size(ZZ)), ZZ=zeros(size(Mod.M)-[0 1]); end
xm=Mod.x(1:end-1)+diff(Mod.x)/2;
zz=interp1(N.elec(:,1),N.elec(:,2),xm,'nearest','extrap');
for i=1:length(zz),
   k=max(find(Mod.z<=zz(i)))-1;
   if k>0, 
       FIX(i,1:min(k,nmz))=-1; 
%        ZZ(i,min(k,nmz-1))=1;
%        ZZ(i,min(k,nmz-2)+(0:1))=1;
       dz=Mod.z(k+1)-Mod.z(k);
       if (zz(i)-Mod.z(k+1))/dz>0.1, ZZ(i,min(k+1,nmz-1))=1; end
       if (Mod.z(k+2)-zz(i))/dz>0.1, ZZ(i,min(k,nmz-1))=1; end
   end
end