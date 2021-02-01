global N datfile
[fp,fn,fe]=fileparts(datfile);
T=rmfield(N,'topo');
T.elec=mbm2xz(N.elec(:,1),N.topo,1,3);
T.r=N.r./N.k;
saveinv2dfile(strrep(datfile,fe,'.ohm'),T);