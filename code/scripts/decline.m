global Mod ZZ
nmz=size(Mod.M,2);
if ~isequal(size(Mod.M)-[0 1],size(ZZ)), ZZ=zeros(size(M)-[0 1]); end
xm=Mod.x(1:end-1)+diff(Mod.x)/2;
zz=interp1(A(:,1),A(:,2),xm,'nearest','extrap');
for i=1:length(zz),
   k=max(find(Mod.z<=zz(i)))-1;
   if k>0, 
       dz=Mod.z(k+1)-Mod.z(k);
       if (zz(i)-Mod.z(k+1))/dz>0.4, ZZ(i,min(k+1,nmz-1))=1; end
       if (Mod.z(k+2)-zz(i))/dz>0.4, ZZ(i,min(k,nmz-1))=1; end
   end
end
