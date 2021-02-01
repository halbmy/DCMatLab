datfile='D:\Guenther.T\2d\ginderes\e1000n\e1000n_pd1_topo.dat';
N=readinv2dfile(datfile);
plot(N.topo(:,1),N.topo(:,2),'.-');
to=N.topo(1,1)+[0;cumsum(sqrt(sum(diff(N.topo).^2,2)))]; %topo auf in mbm
tnr=105;enr=300; % Eichpunkt
[N.topo(tnr,:) N.elec(enr,:)] % stand vorher, gleiche dx versch. höhen
dx=N.topo(tnr,1)-to(tnr); %korrektur (x)
nelec=N.elec;nelec(:,1)=nelec(:,1)-dx;
elec=topowickel(nelec,N.topo);
plot(N.topo(:,1),N.topo(:,2),'.-');axis equal tight;%ylim([200 235]);
hold on;plot(elec(:,1),elec(:,2),'r.','MarkerSize',2);hold off
[N.elec(enr,:) elec(enr,:)]
del=sqrt(sum(diff(elec).^2,2));
return
N.elec=elec;N=rmfield(N,'topo');
saveinv2dfile(strrep(datfile,'.dat','1.dat'),N,0);
N.r=N.r./N.k;
saveinv2dfile(strrep(datfile,'.dat','1.ohm'),N,0);
