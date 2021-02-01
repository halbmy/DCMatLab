function S = loadsens(sensname)

% LOADSENS - Load sensitivity from binary file
% S = loadsens(sensname)

if nargin<1, sensname=['tmp' filesep 'sensMat.mat']; end
fid=fopen(sensname,'r');
ndata=fread(fid,1,'long')
nmodel=fread(fid,1,'long')
S=zeros(ndata,nmodel);
for i=1:ndata, S(i,:)=fread(fid,nmodel,'double')'; end
fclose(fid);