function [t,dep]=readmon(filename,isxtra)
if nargin<2, isxtra=0; end
fid=fopen(filename);
zeile=fgetl(fid);i=1;
while (length(zeile)<6)||(~strcmp(lower(zeile(1:6)),'[data]')), 
    zeile=fgetl(fid);i=i+1; 
end
ndata=str2num(fgetl(fid));i=i+1;
% A=mytextscan(fid,'%*s%s%f%f',ndata);
fclose(fid);
formstr='%*s%s%f%f';
if isxtra, formstr=[formstr '%*f']; end
[Time,dep,cond]=textread(filename,formstr,ndata,'headerlines',i);
t=datenum(Time);
dep=dep/100;
