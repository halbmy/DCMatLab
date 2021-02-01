function N = read3drawfile(datafile)

% READ3DRAWFILE - Read 3D Raw Data File
% N = read3drawfile(filename);
% filename must consist of the columns
% x_A y_A x_B y_B x_M y_M x_N y_N rho_a

fid=fopen(datafile,'r');
if fid<0, error(['File ' datafile ' not found!']); end
nr='';headerlines=-1;
while isempty(nr),
    zeile=fgetl(fid);
    nr=sscanf(zeile,'%f');
    headerlines=headerlines+1;
end
fclose(fid);
ncols=length(str2num(zeile));
fid=fopen(datafile,'r');
for i=1:headerlines, zeile=fgetl(fid); end
DATA=fscanf(fid,'%f',[ncols Inf])';
fclose(fid);
nxa=1;nxb=3;nxm=5;nxn=7;nrho=9;nip=0;nerr=0;
if ncols>12, % apparently xyz for each
   nxb=4;nxm=7;nxn=10;nrho=13; 
else
    if ncols>9, nerr=10; end
    if ncols>10, nip=11; end
end
AA=DATA(:,nxa:nxa+1);BB=[];NN=[];
if nxb>0, BB=DATA(:,nxb:nxb+1); end
MM=DATA(:,nxm:nxm+1);
if nxn>0, NN=DATA(:,nxn:nxn+1); end
N.r=DATA(:,nrho);
if nip, N.ip=DATA(:,nip); end
if nerr, N.err=DATA(:,nerr); end
N.elec=unique(sortrows([AA;BB;MM;NN]),'rows');
data=size(DATA,1);
N.b=zeros(data,1);N.n=zeros(data,1);
[tf,N.a]=ismember(AA,N.elec,'rows');
if ~isempty(BB), [tf,N.b]=ismember(BB,N.elec,'rows'); end
[tf,N.m]=ismember(MM,N.elec,'rows');
if ~isempty(NN), [tf,N.n]=ismember(NN,N.elec,'rows'); end
N.elec(:,3)=0;
N.a=N.a(:);N.b=N.b(:);N.m=N.m(:);N.n=N.n(:);N.r=N.r(:);
N.k=getkonf(N);