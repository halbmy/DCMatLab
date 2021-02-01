function S = loadsens(sensname,inttype)

% LOADSENS - Load sensitivity from binary file
% S = loadsens(sensname)

if nargin<1, sensname=['tmp' filesep 'sensMat.mat']; end
if nargin<2, inttype='long'; end
fid=fopen(sensname,'r');
ndata=fread(fid,1,inttype);
nmodel=fread(fid,1,inttype);
S=zeros(ndata,nmodel);
for i=1:ndata, S(i,:)=fread(fid,nmodel,'double')'; end
fclose(fid);