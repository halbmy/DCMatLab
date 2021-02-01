function N = readres2dinvfile(datafile)
% READRES2DINVFILE - Read 2D Loke&Barker File
% N = readres2dinvfile('filename.dat');
% N.....structure of * arrays a,b,m,n = electrode numbers(elec)
%                                   k = konfiguration factor
%                    * elec..Electrode positions ( x,z )

N.elec=[];
if nargin>1,
    input=fopen('data/rothdi.dat','r');
else
    input=fopen(datafile,'r');
end
if input<0, 
    message(sprintf('Could not open datafile: %s',datafile));
    return; 
else
    mess1=sprintf('Opening datafile: %s ',datafile);
end

nn=9; % Anzahl Zahlen/Zeile Standard-> 4*2 (x,y) + Rho_a
data=inf;  % Maximal möglich
ippresent=0;
erstezeile=[];
xgrid=0;ygrid=0;xunit=0;yunit=0;
zeile='';
while isempty(destrip(zeile,'//')), zeile=fgetl(input); end
titel=zeile;
zeile=';';while zeile(1)==';', zeile=fgetl(input); end
xgrid=sscanf(zeile,'%f',1);
zeile=';';while zeile(1)==';', zeile=fgetl(input); end
typ=sscanf(zeile,'%d',1);
if ~isnumeric(typ), fclose(input);return; end
isR=0;
if typ==11, %independent electrode positions
    subtype=fscanf(input,'%d\n',1);
    zeile=fgetl(input);
    isR=fscanf(input,'%d',1);
end
data=fscanf(input,'%d',1);
N=struct('a',0,'b',0,'m',0,'n',0,'k',0,'r',0);
location=fscanf(input,'%d',1);
ippresent=fscanf(input,'%d\n',1);
if ippresent==1, % read 3 lines if ip present
    ss=fgetl(input);
    ss=fgetl(input);
    ss=fgetl(input);
end
% if typ==8, % DIPOLE EQUATORIAL
%     mess2=sprintf('Dipol-Dipol-Equatorial NOT yet supported!');
%     return;
% end
if (typ==12)||(typ==13), % CROSS BORHOLE
    %mess2=sprintf('Cross-Borehole NOT yet supported!');
    %return;
    surfel=fgetl(input);
    n0=fscanf(input,'%d\n',1);
    N.elec=[];
    for n=1:n0,
        zeile=fgetl(input);
        el=sscanf(strrep(zeile,',',' '),'%f %f');
        N.elec=[N.elec;el(:)'];
    end
    borel=fgetl(input);
    if strcmp(borel(1:4),'Bore'),
        nbor=2;
    else
        nbor=fscanf(input,'%d\n',1);
    end
    for k=1:nbor,
        nb=[];
        while isempty(nb),
            zeile=fgetl(input);            
            nb=sscanf(zeile,'%d\n',1);
        end
        for n=1:nb,
            zeile=fgetl(input);
            el=sscanf(strrep(zeile,',',' '),'%f %f');
            N.elec=[N.elec;el(:)'];
        end
        borel=fgetl(input);
    end
end
if (typ==11)||(typ==12)||(typ==13), % mixed array
    mess2=sprintf('Mixed Array %d Datum Points',data); 
    AA=[];BB=AA;MM=AA;NN=AA;
    N.r=zeros(data,1);
    if ippresent, N.ip=N.r; end
    for l=1:data,
       zeile='';
       while isempty(zeile), zeile=str2num(fgetl(input)); end
       num=zeile(1);
       N.r(l)=zeile(end-ippresent);
       if ippresent, N.ip(l)=zeile(end); end
       AA=[AA;zeile(2:3)];
       bb=[9999 0];nn=bb;
       if num==2,
           mm=zeile(4:5);
       end
       if num==3,
           mm=zeile(4:5);
           nn=zeile(6:7);
       end
       if num==4,
           bb=zeile(4:5);
           mm=zeile(6:7);
           nn=zeile(8:9);
       end
       BB=[BB;bb];MM=[MM;mm];NN=[NN;nn];
    end
    istopo=fscanf(input,'%d',1);
    if istopo,
        nrtopo=fscanf(input,'%d\n',1);
        for nt=1:nrtopo,
            zeile=strrep(fgetl(input),',',' ');
            N.topo(nt,1:2)=sscanf(zeile,'%f',2);
        end
        %        N.topo=fscanf(input,'%f',[2 nrtopo])';
    end
else
    nntyp=[3 3 4 3 3 4 4 3 0 0 8 10]; %No. of columns
    nn=nntyp(typ);
    ss='%f';
    for nnn=2:nn+ippresent
        ss=strcat(ss,' %f');
    end
    ss=strcat(ss,'\n');
    %alles=fscanf(input,'%s\n',data);
    %[DATA,data]=sscanf(alles,ss,[nn+ippresent,data]);
    %[DATA,data]=fscanf(input,ss,[nn+ippresent,data]);
    DATA=zeros(nn+ippresent,data);
    iserr=0;
    for l=1:data,
        zeile=fgetl(input);
        if (l==1)&&any(strfind(zeile,'Error')),
            iserr=1;
            DATA(end+1,1)=0;
            ss=['%f ' ss];
            zeile=fgetl(input); % description
            errtype=str2num(fgetl(input)); % 
            zeile=fgetl(input); % first data row
        end
        zeile(zeile==',')=' ';
        zeile(zeile==';')='';
        if l==1,
            vals=str2num(zeile);
            lvals=length(vals);
            DATA(lvals,1)=0;
        end
        vals=sscanf(zeile,'%f',[lvals,1]);
        if length(vals)==size(DATA,1), DATA(:,l)=vals; end
%         DATA(:,l)=sscanf(zeile,ss,[nn+ippresent+iserr,1]);
    end
    istopo=fscanf(input,'%d',1);
    if istopo,
       nrtopo=fscanf(input,'%d\n',1);
       for nt=1:nrtopo,
          zeile=strrep(fgetl(input),',',' ');
          N.topo(nt,1:2)=sscanf(zeile,'%f',2);
       end
%        N.topo=fscanf(input,'%f',[2 nrtopo])';
    end
    %data=fix(data/(nn+ippresent));
    XX=DATA(1,:)';
    EL=DATA(2,:)';
    N.r=DATA(nn,:)';
    sida=size(DATA,1);
    simu=nn+ippresent+iserr;
    if sida>simu, N.i=DATA(simu+1,:)'/1000; end
    if sida>simu+1, N.u=DATA(simu+2,:)'/1000; end
    if ippresent, N.ip=DATA(nn+1,:)'; end
    if nn==4, % 4-Point-Measurement
        SP=round(DATA(3,:)*1200)'/1200; % thirds etc.
    end
%     if isip, N.ip=DATA(nn+ippresent,:)'; end
    if iserr, 
        N.err=DATA(nn+ippresent+1,:)';
        if errtype==0, N.err=N.err./N.r; end
    end
    BB=[];NN=[];
    if typ==1, % WENNER
        mess2=sprintf('Wenner Array %d Datum Points',data);
        AA=XX-location*EL*1.5;
        MM=AA+EL;
        NN=MM+EL;
        BB=NN+EL;
    end
    if typ==2, % POLE-POLE
        mess2=sprintf('Pole-Pole Array %d Datum Points',data); 
        AA=XX-location*EL*0.5;
        MM=AA+EL;
    end
    if typ==3, % DIPOLE-DIPOLE
        mess2=sprintf('Dipole-Dipole Array %d Datum Points',data); 
        AA=XX-location*EL.*(SP/2+1);
        BB=AA+EL;
        MM=BB+SP.*EL;
        NN=MM+EL;
    end
    if typ==4, % WENNER-BETA
        mess2=sprintf('Wenner-Beta Array %d Datum Points',data);
        AA=XX-location*EL*1.5;
        BB=AA+EL;
        MM=BB+EL;
        NN=MM+EL;
    end
    if typ==5, % WENNER-GAMMA
        mess2=sprintf('Wenner-Gamma Array %d Datum Points',data);
        AA=XX-location*EL*1.5;
        MM=AA+EL;
        BB=MM+EL;
        NN=BB+EL;
    end
    if typ==6, % POLE-DIPOLE
        mess2=sprintf('Pole-Dipole Array %d Datum Points',data); 
        AA=XX-location*SP.*EL-(SP<0).*(SP-1).*EL;
        MM=AA+SP.*EL;
        NN=MM+sign(SP).*EL;
%         fi=find(SP<0);du=MM(fi);MM(fi)=NN(fi);NN(fi)=du;
    end
    if typ==7, % SCHLUMBERGER
        mess2=sprintf('Schlumberger Array %d Datum Points',data);
        AA=XX-location*EL.*(SP+0.5);
        MM=AA+SP.*EL;
        NN=MM+EL;
        BB=NN+SP.*EL;
    end
    if typ==8, % HALFWENNER
        
    end
end
fclose(input);
fak=200;
AA=round(AA*fak)/fak;BB=round(BB*fak)/fak;
MM=round(MM*fak)/fak;NN=round(NN*fak)/fak;
% message(strcat(mess1,' --> ',titel));
% message(sprintf('%s min=%.1f,max=%.1f',mess2,min(N.r),max(N.r)));
N.elec=unique(sortrows([AA;BB;MM;NN]),'rows');
anzel=size(N.elec,1);
if size(N.elec,2)<2, N.elec(:,2)=0; end
N.a=zeros(data,1);
N.b=zeros(data,1);
N.m=zeros(data,1);
N.n=zeros(data,1);
N.k=zeros(data,1);

fac=111111;
qel=N.elec(:,1)+N.elec(:,2)*fac;
qa=AA(:,1);qm=MM(:,1);
if size(AA,2)>1, qa=qa+AA(:,2)*fac; end
if size(MM,2)>1, qm=qm+MM(:,2)*fac; end
[jj,ii]=ismember(qa,qel);
% [aa,bb]=meshgrid(qa,qel);
% [ii,jj]=find((aa-bb)==0);
N.a(jj)=ii;
if ~isempty(BB),
    qb=BB(:,1);
    if size(BB,2)>1, qb=qb+BB(:,2)*fac; end
    [jj,ii]=ismember(qb,qel);
%     [aa,bb]=meshgrid(qb,qel);
%     [ii,jj]=find((aa-bb)==0);
    N.b(jj)=ii;
end
[jj,ii]=ismember(qm,qel);
% [aa,bb]=meshgrid(qm,qel);
% [ii,jj]=find((aa-bb)==0);
N.m(jj)=ii;
if ~isempty(NN),
    qn=NN(:,1);
    if size(NN,2)>1, qn=qn+NN(:,2)*fac; end
    [jj,ii]=ismember(qn,qel);
%     [aa,bb]=meshgrid(qn,qel);
%     [ii,jj]=find((aa-bb)==0);
    N.n(jj)=ii;
end

if (typ==11)||(typ==12), % Korrektur der unendlichen Elektroden von Typ 11
    num=find(N.elec(:,1)==9999);
    if ~isempty(num),
        N.b(N.b==num)=0;
        fb=find(N.b>num);
        N.b(fb)=N.b(fb)-1;
        N.n(N.n==num)=0;
        fn=find(N.n>num);
        N.n(fn)=N.n(fn)-1;
        N.elec(num,:)=[];
    end
end
N=delmeasurement(N,N.r==0);
N.k=getkonf(N);
if isR|(typ==13),
    N.rho=N.r;
    N.r=N.rho.*N.k;
end
