function exportres2dinv(N,basename)

% EXPORTRES2DINV - Export res2dinv file(s)
% exportres2dinv(N,basename)
% file(s) are named basename_dd, _we etc.

if nargin<2, basename='aaa'; end
if nargin<1, error('Specify data structure!'); end
[fpath,fname,fext]=fileparts(basename);
basename=strrep(basename,fext,'');
[mids,seps,ii,kk]=midkonf2d(N);
del=min(diff(unique(N.elec(:,1))));
type=fix(seps/10000)+1;
dlen=fix(mod(seps(kk),10000)/100)+1;
nn=mod(seps(kk),100);
addname={'pp','pd','dp','we','sl','dd'};
loketype=[2 6 6 1 7 3 3];
for tt=unique(type),
    fid=fopen([basename '_' addname{tt} '.dat'],'w');
    fprintf(fid,'Test\n');
    fprintf(fid,'%.1f\n',del);
    fprintf(fid,'%d\n',loketype(tt));
    fprintf(fid,'%d\n',length(N.r));
    ss='%g\t%g\t';
    if tt==2,
        fprintf(fid,'%d\n',0);
        aa=N.elec(N.a,1);
    else
        fprintf(fid,'%d\n',1);
        aa=mids(ii);aa=aa(:);
    end
    fprintf(fid,'%d\n%d\n',1,0);
    if ~ismember(tt,[1 4]), 
        ss=[ss '%g\t']; 
        aa=[aa dlen(:)];
    end
    ss=[ss '%.2f\n'];
    aa=[aa nn(:)];
    aa=[aa N.r(:)];
    fprintf(fid,ss,aa');
    for i=1:4, fprintf(fid,'%d\n',0); end
    fclose(fid);
end