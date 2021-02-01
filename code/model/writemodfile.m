function writemodfile(sond,filename,res,thk,lam)

% WRITEMODFILE - Write mod file for EM1dInv use
% writemodfile(sounding_numbers,filename,res,thk)

% filename='test.mod';sond=1;
% filename='test3.mod';sond=1:3;
if nargin<1, sond=1; end
if nargin<2, filename='test.mod'; end
if nargin<3, res=[0.3 1.0]; end
if nargin<4, thk=[5 1 2 3 4]; end
if nargin<5, lam=0.1; end
mode1=1;
comment='test';
fixedlayers=1;
fixfirst=1;
d=cumsum(thk);
cmode=1+(length(sond)>1);cmode=0;%!!!
niter=10;
while length(res)<=length(thk), res(end+1)=res(end); end
nlay=length(res);
cres=ones(size(res))*(-1);
cthk=ones(size(thk))*(-1);
if ~fixedlayers, cres(:)=0.1; end
if fixfirst, cres(1)=1e-3;cthk(1)=1e-3; end
% [res(:),[thk(:);0],[d(:);0],cres(:),[cthk(:);0]]
if length(lam)<length(res), lam=ones(size(res))*lam(1); end
fid=fopen(filename,'w');
fprintf(fid,'%s\r\n',comment);
fprintf(fid,'%d %d\t!# of data&constraint mode\r\n',length(sond),cmode);
for i=1:length(sond), 
    ii=sond(i);
    if mode1, ii=1; end
    fprintf(fid,'%d %d sond%d.dcp\r\n',ii,1,sond(i)); 
end
fprintf(fid,'%d\t!# of iterations\r\n',niter);
maxi=length(sond);if mode1, maxi=1; end
for i=1:maxi,
    fprintf(fid,'%d\t!# of layers\r\n',nlay);
    fprintf(fid,'%g\t%g\t%g\r\n',[res;cres;lam]);
    if fixedlayers, fprintf(fid,'%g\t1.0e-3\t9e9\r\n',thk);
    else fprintf(fid,'%g\t%g\t9e9\r\n',[thk;cthk]); end
    fprintf(fid,'%g\t-1\r\n',d);
end
fclose(fid);
