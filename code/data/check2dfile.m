function erg=check2dfile(fname)

% CHECK2DFILE - Check 2d data file for file type
% type=check2dfile(filename)
% 1 - for res2dinv data files
% 2 - for inv2d data
% 3 - for raw data file
% 4 - resecs data file

fid=fopen(fname,'r');
erg=0;
if fid<0, return; end
zeile='';
while isempty(zeile),
    zeile=destrip(fgetl(fid)); end
if strfind(zeile,'C1(x)'), erg=4;return; end
zeile='';
while isempty(zeile),
    zeile=destrip(fgetl(fid)); end
% while double(zeile(1))>64, zeile(1)=''; end
erg=length(sscanf(zeile,'%f '));
if erg>3,
    erg=3;
else
    if erg>=2, erg=2; end
end
fclose(fid);

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'#');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end
