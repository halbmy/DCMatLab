function bor2xz(Bor,filename,N)

% BOR2XZ - Convert Bore hole information as from reading *.bor file
%          into *.xz file for constraining the dcfemlib model
% bor2xz(Bor,filename[,N])
% Bor..borehole structure with positions (Bor.pos) and layers(Bor.lay)
% N.topo is used (if given) to add topography

if nargin<2, filename='data.xz'; end
if nargin<3,
    posz=zeros(size(Bor.pos));
else
    if isfield(N,'topo'),
        posz=interp1(N.topo(:,1),N.topo(:,2),Bor.pos);
    else
        posz=interp1(N.elec(:,1),N.elec(:,2),Bor.pos);
    end
end
% filename='wob-all.xz';
% fid=1;
fak=0.5;
% clf;plot(Bor.pos,posz);
fid=fopen(filename,'w');
for i=1:length(Bor.pos),
    for j=1:length(Bor.lay{i}),
        if i*j>1, fprintf(fid,'\n'); end % empty row for new line
        xx=Bor.pos(i)+[-1 1]*Bor.lay{i}(j)*fak;
        zz=[1 1]*(posz(i)-Bor.lay{i}(j));
%         line(xx,zz);
        fprintf(fid,'%g\t%g\n',[xx;zz]);
    end
end
fclose(fid);