Topo=load(topofile);
xzr=load(xzrfile);
xhrfile=[strrep(xzrfile,'.xzr','') '.xhr'];
xm=(xzr(:,1)+xzr(:,2))/2;
zm=(xzr(:,3)+xzr(:,4))/2;
xz=interp1(Topo(:,1),Topo(:,2),xm);
xz(find(xm<Topo(1,1)))=Topo(1,2);
xz(find(xm>Topo(end,1)))=Topo(end,2);
xhr=xm;
xhr(:,2)=xz-zm;
xhr(:,3)=xzr(:,5);
fid=fopen(xhrfile,'w');
fprintf(fid,'%.2f\t%.2f\t%.4f\n',xhr');
fclose(fid);