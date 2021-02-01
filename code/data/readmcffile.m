function N=readmcffile(mcffile,txtfile,nx,ny,dx,dy,x0,y0,ismeand)

% mcffile='d:\Guenther.T\3d\jan\d16reip.mcf';
% mcffile='d:\Guenther.T\3d\jan\d16aeq.mcf';
% nx=16;ny=16;dx=0.1;dy=0.1;x0=6.74;y0=6.8;

if nargin<3, nx=16; end
if nargin<4, ny=16; end
if nargin<5, dx=0.1; end
if nargin<6, dy=dx; end
if nargin<7, x0=6.74; end
if nargin<8, y0=6.8; end
if nargin<9, ismeand=0; end

if ~exist(mcffile,'file'), display('Could not find mcf file!');return; end
% txtfile=strrep(mcffile,'.mcf','.txt');
if ~exist(txtfile,'file'), 
%     txtfile=strrep(txtfile,'.txt','ASCII\.txt');%check in ASCII dir
    display('Could not find txt file!');return; 
end
N=[];
x=x0+(0:nx-1)'*dx;
y=y0+(0:ny-1)'*dy;
N.elec=zeros(nx*ny,3);
for i=1:length(x),
    idx=(i-1)*ny+(1:ny)';
    N.elec(idx,1)=x(i);
    N.elec(idx,2)=y;
end
N.a=zeros(20000,1);N.b=N.a;N.m=N.a;N.n=N.a;
fid=fopen(mcffile,'r');
zeile=fgetl(fid);
while zeile(1)=='A', zeile=fgetl(fid); end
l=0;
while isstr(zeile),
    if zeile(1)=='M',
        aa=str2num(zeile(3:end));
        nch=fix((length(aa)-2)/2);
        for i=1:nch,
            l=l+1;
            N.a(l)=aa(1)+1;
            N.b(l)=aa(2)+1;
            N.m(l)=aa(nch*2+1)+1;
            N.n(l)=aa(nch*2+2)+1;
        end
    end
    zeile=fgetl(fid);
end
N.a=N.a(1:l);N.b=N.b(1:l);N.m=N.m(1:l);N.n=N.n(1:l);
fclose(fid);
fid=fopen(txtfile,'r');
zeile=fgetl(fid);
iu=0;ii=0;id=0;ip=0;i=0;ival=0;
while ~isempty(zeile),
    i=i+1;
    [tok,zeile]=strtok(zeile);
    if strcmp(tok,'D'), id=i-1;ival=ival+1; end
    if strcmp(tok,'P'), ip=i-1;ival=ival+1; end
    if strcmp(tok,'U'), iu=i-1;ival=ival+1; end
    if strcmp(tok,'I'), ii=i-1;ival=ival+1; end
end
formstr='%*s';
for lauf=1:i-1, formstr=[formstr '%f']; end
A=mytextscan(fid,formstr);
fclose(fid);
if iu, N.u=A{iu}/1000; end
if ii, N.i=A{ii}/1000; end
if id, N.err=A{id}/100; end
if ip, N.ip=A{ip}; end
lu=length(N.u);
if lu<length(N.a), 
    N.a(lu+1:end)=[];
    N.b(lu+1:end)=[];
    N.m(lu+1:end)=[];
    N.n(lu+1:end)=[];
end
N.k=getkonf(N);
if length(N.u)==length(N.k),
    N.r=N.u./N.i.*N.k;
else
    display('lengths mismatch!');
end
