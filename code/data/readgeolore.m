function Geolore = readgeolore(geolorefile,gpslogfile,gpsfile)
if nargin<3, gpsfile=''; end

A=load(geolorefile);
Geolore=[];
Geolore.time=A(:,2)/24;
Geolore.data=A(:,3:end);
% Geolore.depth=(A(:,3)-0.456)/0.043;
if exist(gpsfile,'file')
    ds=textread(gpsfile,'%*s%s%*s',1,'headerlines',2);
    ts=textread(gpsfile,'%*s%s%*s',1,'headerlines',3);
%     tt=textread(gpsfile,'%*s%f%*s',1,'headerlines',3);
%     nulltime=floor(tt/10000)/24+floor(mod(tt,10000)/100)/24/60+mod(tt,100)/86400;
    nulltime=datenum([ds{1} ts{1}],'ddmmyyHHMMSS');    
    Geolore.time=Geolore.time+nulltime;
    Geolore.gpstime=nulltime;
    [ew,Geolore.lon]=textread(gpsfile,'%*s%s%f','headerlines',4);
    [ns,Geolore.lat]=textread(gpsfile,'%*s%s%f','headerlines',5);
end
if exist(gpslogfile,'file')
    [S1,S2,datum,zeit,Geolore.lat,Geolore.lon]=textread(gpslogfile,...
        '%d/%d%d%d%*s%f%*s%f');
    Geolore.gpstime=floor(zeit/10000)/24+...
        floor(mod(zeit,10000)/100)/24/60+mod(zeit,100)/86400;
    dl=sqrt(diff(Geolore.lat).^2+diff(Geolore.lon).^2);
    Geolore.x=cumsum([0;dl])*1000;
    Geolore.v=dl./(diff(Geolore.gpstime)*24);
end

