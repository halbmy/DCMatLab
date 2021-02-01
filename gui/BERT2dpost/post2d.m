function [Mesh,N]=post2d()

if ~exist('filename','var'), filename='*.zip;*.cfg'; end
if exist('post2d.last'), 
    fid=fopen('post2d.last','r');
    filename=fgetl(fid);
    fclose(fid);
end
[fname,pname]=uigetfile({'*.zip;*.cfg';'*.zip';'*.cfg'},'Choose result',filename);
if ~isstr(fname), return; end
filename=fullfile(pname,fname)
for i=1:4, if ishandle(i), close(i); end; end
[Mesh,N]=postmodel2d(filename);
out=sprintf('%d iterations: RMS=%.1f%%, Chi^2=%.1f\n',Mesh.iter,rms(N.r,N.response),chi2(N.r,N.response,N.err,1));
fprintf(out);
basename=strrep(strrep(strrep(fname,'.zip',''),'.cfg',''),'.','_');
if isequal(questdlg(out,'Export figures?'),'Yes'),
    if ishandle(1), epsprint(1,[pname filesep basename '-model'],1); end
    if ishandle(2), epsprint(2,[pname filesep basename '-data'],1); end
    if ishandle(3), epsprint(3,[pname filesep basename '-topoeff'],1); end
    if ishandle(4), epsprint(4,[pname filesep basename '-ipmodel'],1); end
end
fid=fopen('post2d.last','w');
fprintf(fid,'%s\n',filename);
fclose(fid);
