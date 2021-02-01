function N = readres3dinvfile(datafile,raw)

% READRES3DINVFILE - Read 3D Loke&Barker File
% N = readres3dinvfile('filename.dat');
% N.....structure of arrays: a,b,m,n = electrode numbers(elec)
%             r(+ip) = measurements  k = konfiguration factor 
%       elec = Electrode position ( x,y,z )

if nargin<2, raw=0; end
input=fopen(datafile,'r');

if input<0, 
    message(sprintf('Could not open datafile: %s',datafile));
    return; 
end
message(['Reading file ' datafile]);
nn=9; % Anzahl Zahlen/Zeile Standard-> 4*2 (x,y) + Rho_a
data=inf;  % Maximal möglich
ippresent=0;
erstezeile=[];
xgrid=0;ygrid=0;xunit=0;yunit=0;
startrow=1;
if ~raw,    %% RES3DINV-File!
    titel=fgetl(input);
    xgrid=fscanf(input,'%d\n',1);
    ygrid=fscanf(input,'%d\n',1);
    xunit=fscanf(input,'%f\n',1);
    yunit=fscanf(input,'%f\n',1);
    typ=fscanf(input,'%d\n',1);
    data=fscanf(input,'%d\n',1);
    if typ==1,
        message(sprintf('Wenner Array %d Datum Points',data));
        nn=9;
    end
    if typ==2,
        message(sprintf('Pole-Pole Array %d Datum Points',data)); 
        nn=5;
    end
    if typ==3,
        message(sprintf('Dipole-Dipole Array %d Datum Points',data)); 
        nn=9;
    end
    if typ==6,
        message(sprintf('Pole-Dipole Array %d Datum Points',data)); 
        nn=7;
    end
    if typ==7,
        message(sprintf('Schlumberger Array %d Datum Points',data));
        nn=9;
    end
    if typ==8,
        message(sprintf('Equatorial Dipole-Dipole Array %d Datum Points',data)); 
        nn=9;
    end
    isresistance=0;
    if typ==11, %general array
        if data==0,
            zeile=fgetl(input);
            isresistance=fscanf(input,'%d\n',1);
            data=fscanf(input,'%d\n',1);
        end
        message(sprintf('General Array %d Datum Points',data)); 
        nn=10;startrow=2;
    end
    zeile=fgetl(input);
    if findstr(zeile,'IP')
      ippresent=1;
    else
      erstezeile=str2num(zeile);
    end
    if ippresent==1 % 4 Zeilen zusätzlich
        for iii=1:3, ttt=fgetl(input); end
    end
else
    message(sprintf('Variable array %d datum points',data));
end
ss='%f';
for nnn=2:nn+ippresent
    ss=strcat(ss,' %f');
end
ss=strcat(ss,'\n');
[DATA,data]=fscanf(input,ss,[nn+ippresent,data-1+ippresent]);
if ~isempty(erstezeile)
  DATA=[erstezeile' DATA];
  data=data+nn+ippresent;
end
%DATA=sortrows(round(DATA'*10)/10,[1 2])';
DATA=round(DATA*100)/100;
zeile=fgetl(input);
if isequal(lower(zeile(1:min(4,length(zeile)))),'topo'),
    num=fscanf(input,'%d\n',1);
    N.topom=zeros(xgrid,ygrid);
    for i=1:xgrid,
        for j=1:ygrid,
            N.topom(i,j)=fscanf(input,'%f',1);
        end
    end
end
fclose(input);
data=fix(data/(nn+ippresent));
%message(sprintf('Read %d datum points from data file',data));

BB=[];NN=[];
nnn=startrow;
AA=DATA(nnn:nnn+1,:)';AA(1,3)=0; % 1. Spalte
nnn=nnn+2;
if nn>8, % 4Punkt
    BB=DATA(nnn:nnn+1,:)';
    BB(1,3)=0;
    nnn=nnn+2; 
end
MM=DATA(nnn:nnn+1,:)';
MM(1,3)=0;
nnn=nnn+2;
if nn>6, % 
    NN=DATA(nnn:nnn+1,:)';
    NN(1,3)=0; 
end
if isresistance,
    N.rho=DATA(nn,:)';
else
    N.r=DATA(nn,:)'; %%Data
end
if ippresent,
    N.ip=DATA(nn+1,:)';
end
N.elec=unique(sortrows([AA;BB;MM;NN]),'rows');
anzel=length(N.elec);
if length(N.elec(1,:))<3, N.elec(:,3)=0; end
if 1,
    N.b=zeros(data,1);N.n=zeros(data,1);
    [tf,N.a]=ismember(AA,N.elec,'rows');
    if ~isempty(BB), [tf,N.b]=ismember(BB,N.elec,'rows'); end
    [tf,N.m]=ismember(MM,N.elec,'rows');
    if ~isempty(NN), [tf,N.n]=ismember(NN,N.elec,'rows'); end
else
    N.a=zeros(data,1);
    N.b=N.a;
    N.m=N.a;
    N.n=N.a;
    aller=fix(data/25);
    mal=aller;
    wb=waitbar(0,'Reordering Datum points...');
    for l = 1:data,
        RRR=DATA(:,l);lr=length(RRR)+1-startrow;
        RRR(lr-ippresent:lr)=[];
        lr=length(RRR);
        RRR=reshape(RRR,2,lr/2)';
        RRR(:,3)=0;
        [RRR,ind]=sortrows(RRR,[1 2]);
        %RRR=[AA(l,:);BB(l,:);MM(l,:);NN(l,:)];
        [nn,n1,n2]=intersect(N.elec,RRR,'rows');
        ln=length(n1);
        nn=ind(ind(n2));
        if ln>3,
            N.a(l)=n1(nn(1));N.b(l)=n1(nn(2));N.m(l)=n1(nn(3));N.n(l)=n1(nn(4));
        elseif ln>2,
            N.a(l)=n1(nn(1));N.m(l)=n1(nn(2));N.n(l)=n1(nn(3));
        elseif ln>1,
            N.a(l)=n1(nn(1));N.m(l)=n1(nn(2));
        else
            strcat(sprintf('Error(l=%d): n1=',l),sprintf('%d ',n1));
        end
        mal=mal-1;
        if mal==0,
            waitbar(l/data,wb);
            mal=aller;
        end
    end
    close(wb);
end

if(xgrid*ygrid*xunit*yunit>0),  %% Grid present
%   message(sprintf('Found %dx%d Grid dx=%.1f dy=%.1f',xgrid,ygrid,xunit,yunit));
  minx=min(AA(:,1));
  miny=min(AA(:,2));
  D=min(xunit,yunit);
  %X=(0:xgrid-1)*xunit+minx;
  %Y=(0:ygrid-1)*yunit+miny;
  Z=[];
end

% Konfigurationsfaktoren
N.k=getkonf(N);
if isresistance, N.r=N.rho.*N.k; end