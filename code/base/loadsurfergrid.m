function [A,x,y]=loadsurfergrid(fname);

% LOADSURFERGRID - Load surfer grid (*.grd) file
% [A,x,y] = loadsurfergrid(filename);

fid=fopen(fname,'r');
if fid<0, display('Could not open file'); return; end
% Header section
dsrb=fread(fid,4,'char');
nbytes=fread(fid,1,'long');
vers=fread(fid,1,'long');
% Grid section
dirg=fread(fid,4,'char');
nbytes=fread(fid,1,'long');
nrow=fread(fid,1,'long');
ncol=fread(fid,1,'long');
xll=fread(fid,1,'double');
yll=fread(fid,1,'double');
xsize=fread(fid,1,'double');
ysize=fread(fid,1,'double');
zmin=fread(fid,1,'double');
zmax=fread(fid,1,'double');
rotation=fread(fid,1,'double');
blankvalue=fread(fid,1,'double');
% Data section
atad=fread(fid,4,'char');
eins=fread(fid,1,'long')
A=fread(fid,[ncol nrow],'double')';
fclose(fid);
x=(0:ncol-1)*xsize+xll;
y=(0:nrow-1)*ysize+yll;
if nargout<1,
    imagesc(x,y,A);axis equal tight;colorbar
end
