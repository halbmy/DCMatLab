function saveresecsfile(filename,N)

% SAVERESECSFILE - Save aata in Resecs export format
% saveresecsfile(filename,N)
% N..data struct with elec,a,b,m,n,...

if nargin<2, error('Specify filename and Data Struct!'); end
if size(N.elec,2)<3, N.elec(1,3)=0; end
ndata=length(N.a);
% Nel=size(N.elec,1);
A=zeros(ndata,12);
A(:,1:3)=N.elec(N.a,1:3);
A(:,4:6)=N.elec(N.b,1:3);
A(:,7:9)=N.elec(N.m,1:3);
A(:,10:12)=N.elec(N.n,1:3);
headr='C1(x)\tC1(y)\tC1(z)\tC2(x)\tC2(y)\tC2(z)\tP1(x)\tP1(y)\tP1(z)\tP2(x)\tP2(y)\tP2(z)';
formstr='%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g';
if isfield(N,'u')&&(length(N.u)==length(N.a)),
    A(:,end+1)=N.u(:)*1000;
    headr=[headr '\tU'];
    formstr=[formstr '\t%g'];
end
if isfield(N,'i')&&(length(N.i)==length(N.a)),
    A(:,end+1)=N.i(:)*1000;
    headr=[headr '\tI'];
    formstr=[formstr '\t%g'];
end
if isfield(N,'r')&&(length(N.r)==length(N.a)),
    A(:,end+1)=N.r(:);
    headr=[headr '\tRho'];
    formstr=[formstr '\t%g'];
end
if isfield(N,'err')&&(length(N.err)==length(N.a)),
    A(:,end+1)=N.err(:)*100;
    headr=[headr '\tD'];
    formstr=[formstr '\t%g'];
end
if isfield(N,'ip')&&(length(N.ip)==length(N.a)),
    A(:,end+1)=N.ip(:);
    headr=[headr '\tM'];
    formstr=[formstr '\t%g'];
end
fid=fopen(filename,'w');
if fid<0, error('Could not write file!'); end
fprintf(fid,[headr '\r\n']);
fprintf(fid,[formstr '\r\n'],A');
fclose(fid);