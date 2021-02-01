function N=readampfile(filename)

% READAMPFILE - Read ABEM multi purpose file (*.amp)
% N = readampfile(filename)

N=[];
fid=fopen(filename,'r');
if isequal(fid,-1), error(['Could not open amp file ' filename '!']);return; end
zeile='bla';first='';ndata=0;
while ~isempty(zeile),
    zeile=fgetl(fid);
    dp=findstr(zeile,':');first='';last='';
    if dp, first=zeile(1:dp(1)-1);last=zeile(dp(1)+1:end); end
    if length(first)>3, first=first(1:4); end
    switch first,
        case 'Rows', 
            nn=sscanf(last,'%d',3); %header/data/topo
            nheader=nn(1);ndata=nn(2);ntopo=nn(3);
        case 'Acqu', acq=sscanf(last,'%d'); %1=SP,2=res,3=IP
        case 'Smal', del=sscanf(last,'%f',1); %electrode distance
        case 'Elec', % Electrode layout, layout string, Abbrev.
            typ=sscanf(last,'%d\t%*s\t%*s');
            bla=sscanf(last,'%*d\t%s\t%*s');
            abb=sscanf(last,'%*d\t%*s\t%s');
    end
end
header=fgetl(fid);
%%
lA=3;lB=4;lM=5;lN=6;lI=7;lU=8;lR=9;lE=10;
lTx=0;lRx=0;lDx=0;
i=0;
%%
while ~isempty(header),
   [tok,header]=strtok(header); 
   i=i+1;
   ltok=lower(tok);
   ltok(6:end)='';
   if isequal(ltok,'tx'), lTx=i; end
   if isequal(ltok,'rx'), lRx=i; end
   if isequal(ltok,'dx'), lDx=i; end
   if isequal(ltok,'i(ma)'), lI=i; end
   if isequal(ltok,'volta'), lU=i; end
   if isequal(ltok,'app.r'), lR=i; end
   if isequal(ltok,'error'), lE=i; end
end
%%
if lTx>0, lA=lTx;lB=0; end
if lRx>0, lM=lRx;lB=0; end
% header auswerten!!!
AA=zeros(ndata,1);BB=AA;MM=AA;NN=AA;
N.r=AA;
for i=1:ndata,
   zeile=fgetl(fid);
   if isequal(zeile,'-1'), error('not enough data present!');return; end
   zeile=destrip(zeile);
   aa=str2num(zeile);
   if lA, AA(i)=aa(lA); end
   if lTx, AA(i)=aa(lTx); end
   if lB, BB(i)=aa(lB); end
   if lM, MM(i)=aa(lM); end
   if lRx, MM(i)=aa(lRx)+aa(lTx); end
   if lN, NN(i)=aa(lN); end
   if lR, N.r(i)=aa(lR); end
   if lU, N.u(i)=aa(lU); end
   if lI, N.i(i)=aa(lI)/1000; end
   if lE, N.err(i)=aa(lE)/100; end
   if (lTx>0)&&(lDx>0), BB(i)=AA(i)+aa(lDx); end
   if (lRx>0)&&(lDx>0), NN(i)=MM(i)+aa(lDx); end
end
fclose(fid);
AA=AA*del;BB=BB*del;MM=MM*del;NN=NN*del;    
% fi=find(isinf(MM)|isinf(NN));
% if ~isempty(fi),
%     AA(fi)=[];BB(fi)=[];MM(fi)=[];NN(fi)=[];
%     if lR, N.r(fi)=[]; end
%     if lU, N.u(fi)=[]; end
%     if lI, N.i(fi)=[]; end
%     if lE, N.err(fi)=[]; end
% end
[N.elec,SI,SJ]=unique([AA;BB;MM;NN],'rows');
[TF,N.a]=ismember(AA,N.elec);
[TF,N.b]=ismember(BB,N.elec);
[TF,N.m]=ismember(MM,N.elec);
[TF,N.n]=ismember(NN,N.elec);
if isinf(N.elec(end,1)),
    le=size(N.elec,1);
    N.a(N.a==le)=0;
    N.b(N.b==le)=0;
    N.m(N.m==le)=0;
    N.n(N.n==le)=0;
    N.elec(end,:)=[];
end
if size(N,2)<2, N.elec(:,2)=0; end
if isfield(N,'err'), N.err=N.err(:); end
if isfield(N,'u'), N.u=N.u(:); end
if isfield(N,'i'), N.i=N.i(:); end
if isfield(N,'ip'), N.ip=N.ip(:); end
N.k=getkonf(N);
mi=min(N.err);
N.err=N.err-mi+0.02; % noch genauer zu prüfen!

function zeile=destrip(zeile)
% strip string from comments (with # character)
aa=strfind(zeile,'//');
if ~isempty(aa), zeile=zeile(1:aa(1)-1); end