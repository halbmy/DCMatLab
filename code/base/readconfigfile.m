function BB = readconfigfile(filename)

% READCONFIGFILE = reads config file into struct variable
% STRUC = readconfigfile(filename)
% Example: file looks like
% a=1
% b=23.1
% Result is struct('a',1,'b',23.1)
% STRUC =
%   a: 1
%   b: 23.1

BB=[];
fid=fopen(filename,'r');
if fid==-1, return; end
zeile=destrip(fgetl(fid));
while isstr(zeile),
    [T,R]=strtok(zeile,'=');
    if ~isempty(R),
        R(1)=[];
        D=str2num(R);
        if ~isempty(D), BB=setfield(BB,T,D);
        else BB=setfield(BB,T,R); end
    end
    zeile=destrip(fgetl(fid));
end
fclose(fid);