function All = readgeos1dsounding(file)

% READGEOS1DSOUNDING - Read DC 1D sounding in GEOS format
% Data = readgeos1dsounding(filename)

if nargin<1, file='D:\Guenther.T\pro\hydro\borkum\dc1d\GGA-GGS-BORKUM-1991-001.txt'; end
All=[];
% name=textread(file,'%*s:%s',1,'headerlines',5);
% All.name=name{1};
% All.xy=textread(file,'%*s:%f%*s%*s',2,'headerlines',16)';
% All.rw=All.xy(1);
% All.hw=All.xy(2);
% sau=textread(file,'%*s:%s',1,'headerlines',23);
% All.date=textread(file,'%*s:%s',1,'headerlines',9);
nhl=0;
fid=fopen(file);
zeile=fgetl(fid);
while zeile(1)=='#',
    nhl=nhl+1;    
    if strfind(zeile,'Messung'), All.name=sscanf(zeile,'# %*s : %s'); end
    if strfind(zeile,'Rechtswert'), All.rw=sscanf(zeile,'# %*s : %f'); end
    if strfind(zeile,'Hochwert'), All.hw=sscanf(zeile,'# %*s : %f'); end
    if strfind(zeile,'Startdatum'), 
        All.datestr=sscanf(zeile,'# %*s : %s'); 
        All.date=datenum(All.datestr,'dd.mm.yyyy'); 
    end
    
    zeile=fgetl(fid);
end
fclose(fid);
A=textread(file,'','commentstyle','shell');%'headerlines',nhl);
All.ab2=A(:,1);
All.mn2=A(:,2);
All.rhoa=A(:,3);
if (size(A,2)>3)&&(min(A(:,4))>0), A.err=A(:,4); end
