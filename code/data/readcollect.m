function [MEA,ELPOS]=readcollect(filename,withoutnr)

% READCOLLECT - read collect file
% MEA = readcollect(filename)
% [MEA,elec] = ...
% MEA .. multielectrodename
% elec .. electrode positions

%filename='ggfem3d.collect';
if nargin<2, withoutnr=-1; end
if withoutnr==-1, %autodetect
    fid=fopen(filename,'r');
    zeile=fgetl(fid);
    si=str2num(zeile);
    nel=si(1);
    ELPOS=fscanf(fid,'%f',[3 nel])';
    zeile=fgetl(fid);
    while isempty(zeile), zeile=fgetl(fid); end
    withoutnr=double(length(str2num(zeile))==nel);
    fclose(fid);
end
fid=fopen(filename,'r');
zeile=fgetl(fid);
si=str2num(zeile);
nel=si(1);
ELPOS=fscanf(fid,'%f',[3 nel])';
MEA=fscanf(fid,'%f',[nel+1-withoutnr nel])';
fclose(fid);
if withoutnr==0, MEA(:,1)=[]; end
