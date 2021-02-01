function A = mytextscan(fid,formstr,ndata,iscomment,commstr)

% MYTEXTSCAN - Version of textscan for Matlab<R14
% A = mytextscan(fid,formstr,ndata)
% fid - file identifier of an opened file
% formstr - format string for each line (before # sign)
% ndata - number of lines to read, may be neglected

if nargin<4, iscomment=0; end
if nargin<3, ndata=200000; end

cols=sort([strfind(formstr,'%d') strfind(formstr,'%f')]);
scols=strfind(formstr,'%s');
ncols=length(cols);
for i=1:ncols, A{i}=zeros(ndata,1); end
% if iscomment, zeile=fgetl(fid); end
for i=1:ndata,
    zeile=fgetl(fid);
    if isequal(zeile,-1), %premature end
        for j=1:length(A), A{j}(i:end)=[]; end
        return;
    end
    zeile=strrep(destrip(zeile),',','.');
    zeile=strrep(zeile,'na','NaN');
    while isempty(zeile), zeile=destrip(fgetl(fid)); end
    if ~ischar(zeile), 
        for j=1:length(A), A{j}(i:end)=[]; end
        break; 
    end    
    mess=sscanf(zeile,formstr);
    if length(mess)<ncols, ncols=length(mess); end
    for j=1:ncols, A{j}(i)=mess(j); end    
end

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'#');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end