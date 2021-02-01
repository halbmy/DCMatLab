function N=readdatafile(datfile)

% READDATAFILE - tries to read any data file
if exist(datfile)~=2, error('File does not exist!'); end
[fpath,name,ext]=fileparts(datfile);
lext=lower(ext);
% check for res3dinv or unified data in 3d
is3d=false;
fid=fopen(datfile);
for i=1:4, zeile=fgetl(fid); end
if length(sscanf(destrip(zeile),'%f'))>2, is3d=true; end
fclose(fid);
if ismember(lext,{'.pro','s3d'})||is3d,
    N=read3dfile(datfile);
else
    N=read2dfile(datfile);
end