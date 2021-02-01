function savesensmat(S,sensname)

% SAVESENSMAT - Save sensitivity matrix to binary file
% savesensmat(S,sensname)

if nargin<1, sensname=['tmp' filesep '\sensMat.mat']; end
fid=fopen(sensname,'w');
ndata=size(S,1);nmodel=size(S,2);
fwrite(fid,ndata,'long');
fwrite(fid,nmodel,'long');
for i=1:ndata, fwrite(fid,S(i,:),'double'); end
fclose(fid);