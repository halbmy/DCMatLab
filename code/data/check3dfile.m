function erg=check3dfile(fname)

% CHECK3DFILE - Check 3d data file for file type
%   type=check3dfile(filename)
%     1 - for res3dinv data files
%     2 - for inv3d data
%     3 - for raw data file

fid=fopen(fname,'r');
erg=0;
if fid<0, return; end
zeile='';
while isempty(zeile),
    zeile=destrip(fgetl(fid)); end
if strfind(zeile,'C1(x)'), erg=4;return; end
zeile='';
while isempty(zeile),
    zeile=fgetl(fid);
    if (zeile(1)=='#'), erg=5; return; end
    zeile=destrip(zeile);
end
erg=length(sscanf(zeile,'%f '));
if erg>3,
    erg=3;
else
    if erg>=2, erg=2; end
end

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'#');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end
