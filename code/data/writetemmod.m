function writetemmod(temfile,filename,nlay)

% WRITETEMMOD - Write TEM mod-file for inversion with EM1dInv
% writetemmod(temfile,filename,nlay)

if nargin<3, nlay=4; end
if nargin<2, filename='test.mod'; end
if nargin<1, temfile='test.tem'; end
res=ones(nlay,1)*100;
thk=ones(nlay-1,1)*10;
dep=[0;cumsum(thk)];
fid=fopen(filename,'w');
fprintf(fid,'Comment\n');
fprintf(fid,'1 0\n');
fprintf(fid,'%d %d %s\n',1,1,temfile);
niter=50;
fprintf(fid,'%d\n',niter);
fprintf(fid,'%d\n',nlay);
lam=0.01;
for i=1:nlay, fprintf(fid,'%g\t%g\t%g\n',res(i),-1,lam); end
for i=1:nlay-1, fprintf(fid,'%g\t%g\t%g\n',thk(i),-1,lam); end
for i=1:nlay-1, fprintf(fid,'%g\t%g\n',dep(i),-1); end
fclose(fid);
