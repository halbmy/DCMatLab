function P = loadprimpot(sensname)

if nargin<1, sensname='tmp\secPot.mat'; end
fid=fopen(sensname,'r');
npots=fread(fid,1,'long');
nnodes=fread(fid,1,'long');
P=zeros(npots,nnodes);
for i=1:npots, P(i,:)=fread(fid,nnodes,'double')'; end
fclose(fid);