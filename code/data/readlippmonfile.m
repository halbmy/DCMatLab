function [Data,zeit]=readlippmonfile(filename)
% filename='100210_1305.txt';
% filename='091222_1315.txt';
fid=fopen(filename);
try,
    for i=1:11, zeile=fgetl(fid); end
    del=str2num(zeile);
    zeile=fgetl(fid);
    fel=str2num(zeile);
    zeile=fgetl(fid);
    aa=str2num(zeile);
    nel=aa(2);
    for i=1:3, zeile=fgetl(fid); end
    ndata=str2num(zeile);
    abmn=fscanf(fid,'%d',[4 ndata])';
    for i=1:2, zeile=fgetl(fid); end
    zeit=datenum(zeile,'dd.mm.yyyy HH:MM:SS');
    for i=1:2, zeile=fgetl(fid); end
    A=fscanf(fid,'%f',[6 ndata])';
catch,
    
end
fclose(fid);
%%
Data.elec=zeros(nel,2);
Data.elec(:,2)=-((0:nel-1)'*del+fel);
Data.x=Data.elec(:,1);
Data.z=Data.elec(:,2);
Data.d=-Data.z;
Data.a=abmn(:,1);Data.b=abmn(:,2);Data.m=abmn(:,3);Data.n=abmn(:,4);
Data.i=A(:,3)/1000;
Data.u=A(:,1)/1000;
Data.ip=atan(-A(:,2)./A(:,1));
Data=delmeasurement(Data,(Data.a==255)|(Data.b==255)|(Data.m==255)|(Data.n==255));
Data.k=getkonf(Data);
Data.r=Data.u./Data.i.*Data.k;

